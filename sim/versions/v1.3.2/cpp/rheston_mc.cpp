// rheston_mc.cpp — rough Heston Monte-Carlo, multithreaded, for the parameter
// tables of the ADS2014-style document.
//
// Model (El Euch--Rosenbaum), alpha = H + 1/2 in (1/2, 1]:
//
//     dS_t   = r S_t dt + S_t sqrt(V_t) dW_t
//     V_t    = V_0 + int_0^t K(t-s) [ kappa (theta - V_s) ds
//                                     + eta sqrt(V_s) dB_s ]
//     K(u)   = u^{alpha-1} / Gamma(alpha),     d<W,B>_t = rho dt.
//
// At H = 1/2 the kernel is the constant 1/Gamma(1) = 1 and the equation is the
// classical Heston CIR variance, which is the anchor of the whole document: there
// the Fourier reference reduces to Heston's closed form.
//
// Discretisation.  Euler--Volterra with the LEFT-ENDPOINT kernel whose smallest
// lag is dt, i.e.
//
//     V_{k+1} = V_0 + sum_{j=0}^{k} K((k-j+1) dt) g_j,
//     g_j     = kappa (theta - V_j^+) dt + eta sqrt(V_j^+) dB_j,
//
// the same convention as the project's "exact convolution": no evaluation of K
// at 0, where it is infinite for H < 1/2.  The variance is truncated at zero and
// the number of truncations is reported, because that truncation is the known
// source of the scheme's upward bias.
//
// Two independent passes per parameter set.
//
//   European  no trajectory storage, antithetic pairs, put and call from the
//             same paths.  The reported standard error is that of the
//             antithetic-paired estimator (the pair is one sample).
//
//   American  Longstaff--Schwartz with a REGRESSION sample and a disjoint
//             VALUATION sample.  The policy is fitted on the first and applied
//             on the second, so the reported price carries no in-sample
//             look-ahead bias: it is a genuine lower bound on the Bermudan value
//             for the chosen exercise grid, with an honest standard error.  The
//             in-sample figure is reported too, as the two bracket the estimate.
//
// The American CALL is the control: with r >= 0 and no dividends early exercise
// of a call is never optimal, so its true value is the European call, which the
// Fourier reference gives independently.  Any gap is the cost of the LSM policy,
// measured rather than assumed.
//
// Build:  clang++ -O3 -march=native -std=c++17 -pthread
// Usage:  rheston_mc <config.csv >out.csv     (see run_rheston_tables.py)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

// ------------------------------------------------------------------------- rng
// xoshiro256++ seeded by splitmix64.  One instance per thread per block, seeded
// from (config seed, block index), so the whole run is reproducible and does not
// depend on how the blocks are scheduled across threads.
struct Rng {
    uint64_t s[4];
    explicit Rng(uint64_t seed) {
        for (int i = 0; i < 4; ++i) {
            seed += 0x9E3779B97F4A7C15ULL;
            uint64_t z = seed;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            s[i] = z ^ (z >> 31);
        }
    }
    static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    inline uint64_t next() {
        uint64_t r = rotl(s[0] + s[3], 23) + s[0];
        uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t;    s[3] = rotl(s[3], 45);
        return r;
    }
    // uniform on (0,1), never 0 so log() is safe
    inline double unif() {
        return (double)((next() >> 11) + 1) * (1.0 / 9007199254740992.0);
    }
    inline void norm2(double& a, double& b) {
        double u1 = unif(), u2 = unif();
        double rr = std::sqrt(-2.0 * std::log(u1));
        double th = 6.28318530717958647692 * u2;
        a = rr * std::cos(th);
        b = rr * std::sin(th);
    }
};

// ---------------------------------------------------------------------- config
struct Cfg {
    std::string id;
    double H, V0, theta, kappa, eta, rho, r, T, S0, K;
    int steps;          // Euler--Volterra steps
    long paths_eu;      // European: antithetic PAIRS = paths_eu/2
    long paths_reg;     // American: regression sample
    long paths_val;     // American: valuation sample
    int ex_stride;      // exercise every ex_stride steps
    uint64_t seed;
};

struct EuOut {
    double put = 0, put_se = 0, call = 0, call_se = 0;
    double mean_ST = 0, mean_VT = 0;
    long long neg = 0;
    double secs = 0;
};

struct AmOut {
    double put = 0, put_se = 0, call = 0, call_se = 0;
    double put_insample = 0, call_insample = 0;
    double secs = 0;
    int dates = 0;
};

static int n_threads() {
    unsigned h = std::thread::hardware_concurrency();
    if (h == 0) h = 4;
    if (h > 2) h -= 1;                 // leave one core for the driver
    return (int)h;
}

// ------------------------------------------------------- the path kernel, once
// Simulates ONE path with the supplied standard normals (sign = +1 or -1 for the
// antithetic twin).  If store != nullptr, records (S, V, A) at the exercise
// dates.  Returns S_T; writes the terminal variance and the truncation count.
struct Store {
    float* S; float* V; float* A;      // date-major: [d * npaths + p]
    long npaths; long p; int ndates; int stride;
};

static inline double one_path(const Cfg& c, const double* ker, const double* z1,
                             const double* z2, double sgn, double* g,
                             double& vT, long long& neg, Store* st) {
    const double dt   = c.T / c.steps;
    const double sqdt = std::sqrt(dt);
    const double rc   = std::sqrt(std::max(0.0, 1.0 - c.rho * c.rho));
    double v    = c.V0;
    double logS = std::log(c.S0);
    double A    = 0.0;
    for (int k = 0; k < c.steps; ++k) {
        if (v < 0.0) ++neg;
        const double vp = v > 0.0 ? v : 0.0;
        const double sv = std::sqrt(vp);
        A += vp * dt;
        const double dB = sgn * z1[k] * sqdt;
        const double dW = c.rho * dB + rc * sgn * z2[k] * sqdt;
        logS += (c.r - 0.5 * vp) * dt + sv * dW;
        g[k] = c.kappa * (c.theta - vp) * dt + c.eta * sv * dB;
        if (k + 1 < c.steps) {
            double acc = 0.0;
            for (int j = 0; j <= k; ++j) acc += ker[k - j] * g[j];
            v = c.V0 + acc;
        }
        if (st && ((k + 1) % st->stride == 0)) {
            const int d = (k + 1) / st->stride - 1;
            if (d < st->ndates) {
                const long o = (long)d * st->npaths + st->p;
                st->S[o] = (float)std::exp(logS);
                st->V[o] = (float)(v > 0.0 ? v : 0.0);
                st->A[o] = (float)(A / ((k + 1) * dt));
            }
        }
    }
    vT = v > 0.0 ? v : 0.0;
    return std::exp(logS);
}

// -------------------------------------------------------------- European pass
static EuOut run_european(const Cfg& c) {
    const auto t0 = std::chrono::steady_clock::now();
    const double dt = c.T / c.steps;
    const double alpha = c.H + 0.5;
    std::vector<double> ker(c.steps);
    const double gA = std::tgamma(alpha);
    for (int i = 0; i < c.steps; ++i) ker[i] = std::pow((i + 1) * dt, alpha - 1.0) / gA;

    const long npairs = c.paths_eu / 2;
    const int nt = n_threads();
    std::vector<double> sp(nt, 0.0), sp2(nt, 0.0), sc(nt, 0.0), sc2(nt, 0.0);
    std::vector<double> sst(nt, 0.0), svt(nt, 0.0);
    std::vector<long long> sneg(nt, 0);

    auto work = [&](int tid) {
        std::vector<double> z1(c.steps), z2(c.steps), g(c.steps);
        double lp = 0, lp2 = 0, lc = 0, lc2 = 0, lst = 0, lvt = 0;
        long long lneg = 0;
        for (long i = tid; i < npairs; i += nt) {
            Rng rng(c.seed * 0x9E3779B97F4A7C15ULL + (uint64_t)i);
            for (int k = 0; k < c.steps; ++k) rng.norm2(z1[k], z2[k]);
            double vT = 0;
            const double s_p = one_path(c, ker.data(), z1.data(), z2.data(), +1.0,
                                        g.data(), vT, lneg, nullptr);
            lvt += vT;
            const double s_m = one_path(c, ker.data(), z1.data(), z2.data(), -1.0,
                                        g.data(), vT, lneg, nullptr);
            lvt += vT;
            const double put  = 0.5 * (std::max(c.K - s_p, 0.0) + std::max(c.K - s_m, 0.0));
            const double call = 0.5 * (std::max(s_p - c.K, 0.0) + std::max(s_m - c.K, 0.0));
            lp += put;  lp2 += put * put;
            lc += call; lc2 += call * call;
            lst += 0.5 * (s_p + s_m);
        }
        sp[tid] = lp; sp2[tid] = lp2; sc[tid] = lc; sc2[tid] = lc2;
        sst[tid] = lst; svt[tid] = lvt; sneg[tid] = lneg;
    };
    std::vector<std::thread> th;
    for (int t = 0; t < nt; ++t) th.emplace_back(work, t);
    for (auto& t : th) t.join();

    double P = 0, P2 = 0, C = 0, C2 = 0, ST = 0, VT = 0;
    long long neg = 0;
    for (int t = 0; t < nt; ++t) {
        P += sp[t]; P2 += sp2[t]; C += sc[t]; C2 += sc2[t];
        ST += sst[t]; VT += svt[t]; neg += sneg[t];
    }
    const double n = (double)npairs;
    const double disc = std::exp(-c.r * c.T);
    EuOut o;
    const double mp = P / n, mc = C / n;
    o.put  = disc * mp;
    o.call = disc * mc;
    o.put_se  = disc * std::sqrt(std::max(0.0, P2 / n - mp * mp) / n);
    o.call_se = disc * std::sqrt(std::max(0.0, C2 / n - mc * mc) / n);
    o.mean_ST = ST / n;
    o.mean_VT = VT / (2.0 * n);
    o.neg = neg;
    o.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return o;
}

// ----------------------------------------------------- least-squares machinery
// Solves (X'X + lam I) b = X'y by Cholesky.  nb is small (8), so this is trivial
// work next to the simulation; it is written out to avoid a BLAS dependency.
static bool solve_sym(std::vector<double>& M, std::vector<double>& b, int nb) {
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j <= i; ++j) {
            double s = M[i * nb + j];
            for (int k = 0; k < j; ++k) s -= M[i * nb + k] * M[j * nb + k];
            if (i == j) {
                if (s <= 1e-14) return false;
                M[i * nb + j] = std::sqrt(s);
            } else {
                M[i * nb + j] = s / M[j * nb + j];
            }
        }
    }
    for (int i = 0; i < nb; ++i) {
        double s = b[i];
        for (int k = 0; k < i; ++k) s -= M[i * nb + k] * b[k];
        b[i] = s / M[i * nb + i];
    }
    for (int i = nb - 1; i >= 0; --i) {
        double s = b[i];
        for (int k = i + 1; k < nb; ++k) s -= M[k * nb + i] * b[k];
        b[i] = s / M[i * nb + i];
    }
    return true;
}

static const int NB = 8;
static inline void basis(double S, double V, double A, double K, double theta,
                         double* f) {
    const double x = S / K;
    const double v = V / theta;
    const double a = A / theta;
    f[0] = 1.0; f[1] = x; f[2] = x * x; f[3] = x * x * x;
    f[4] = v;   f[5] = v * v; f[6] = x * v; f[7] = a;
}

// --------------------------------------------------------------- American pass
static AmOut run_american(const Cfg& c) {
    const auto t0 = std::chrono::steady_clock::now();
    const double dt = c.T / c.steps;
    const double alpha = c.H + 0.5;
    std::vector<double> ker(c.steps);
    const double gA = std::tgamma(alpha);
    for (int i = 0; i < c.steps; ++i) ker[i] = std::pow((i + 1) * dt, alpha - 1.0) / gA;

    const int ndates = c.steps / c.ex_stride;
    const long nreg = c.paths_reg, nval = c.paths_val, np = nreg + nval;
    std::vector<float> S((size_t)ndates * np), V((size_t)ndates * np), A((size_t)ndates * np);

    const int nt = n_threads();
    std::vector<long long> sneg(nt, 0);
    auto work = [&](int tid) {
        std::vector<double> z1(c.steps), z2(c.steps), g(c.steps);
        long long lneg = 0;
        // antithetic pairs, laid out so that path 2i and 2i+1 are a pair and both
        // land in the same block (nreg and nval are even)
        for (long i = tid; i < np / 2; i += nt) {
            Rng rng(c.seed * 0xD1B54A32D192ED03ULL + (uint64_t)i + 1);
            for (int k = 0; k < c.steps; ++k) rng.norm2(z1[k], z2[k]);
            double vT = 0;
            Store st{S.data(), V.data(), A.data(), np, 2 * i, ndates, c.ex_stride};
            one_path(c, ker.data(), z1.data(), z2.data(), +1.0, g.data(), vT, lneg, &st);
            st.p = 2 * i + 1;
            one_path(c, ker.data(), z1.data(), z2.data(), -1.0, g.data(), vT, lneg, &st);
        }
        sneg[tid] = lneg;
    };
    {
        std::vector<std::thread> th;
        for (int t = 0; t < nt; ++t) th.emplace_back(work, t);
        for (auto& t : th) t.join();
    }

    const double dtau = c.ex_stride * dt;
    const double dfac = std::exp(-c.r * dtau);

    AmOut o;
    o.dates = ndates;

    // one backward induction per payoff sign: +1 = put, -1 = call
    for (int which = 0; which < 2; ++which) {
        const bool is_put = (which == 0);
        std::vector<double> cash(np);
        const size_t last = (size_t)(ndates - 1) * np;
        for (long p = 0; p < np; ++p) {
            const double s = S[last + p];
            cash[p] = is_put ? std::max(c.K - s, 0.0) : std::max(s - c.K, 0.0);
        }
        for (int d = ndates - 2; d >= 0; --d) {
            const size_t off = (size_t)d * np;
            for (long p = 0; p < np; ++p) cash[p] *= dfac;
            // fit on the regression block, in-the-money paths only
            std::vector<double> M((size_t)NB * NB, 0.0), b(NB, 0.0), f(NB);
            long nin = 0;
            for (long p = 0; p < nreg; ++p) {
                const double s = S[off + p];
                const double ex = is_put ? c.K - s : s - c.K;
                if (ex <= 0.0) continue;
                ++nin;
                basis(s, V[off + p], A[off + p], c.K, c.theta, f.data());
                for (int i = 0; i < NB; ++i) {
                    b[i] += f[i] * cash[p];
                    for (int j = 0; j <= i; ++j) M[i * NB + j] += f[i] * f[j];
                }
            }
            if (nin < 4 * NB) continue;                 // too few: never exercise here
            for (int i = 0; i < NB; ++i) M[i * NB + i] += 1e-8 * M[0] ;
            std::vector<double> Mc = M, bc = b;
            if (!solve_sym(Mc, bc, NB)) continue;
            // apply the SAME policy to both blocks
            for (long p = 0; p < np; ++p) {
                const double s = S[off + p];
                const double ex = is_put ? c.K - s : s - c.K;
                if (ex <= 0.0) continue;
                basis(s, V[off + p], A[off + p], c.K, c.theta, f.data());
                double cont = 0.0;
                for (int i = 0; i < NB; ++i) cont += bc[i] * f[i];
                if (ex > cont) cash[p] = ex;
            }
        }
        const double d0 = std::exp(-c.r * dtau);        // discount tau_0 -> 0
        double sr = 0, sv = 0, sv2 = 0;
        for (long p = 0; p < nreg; ++p) sr += cash[p];
        for (long p = nreg; p < np; ++p) { sv += cash[p]; sv2 += cash[p] * cash[p]; }
        const double mr = d0 * sr / (double)nreg;
        const double mv = d0 * sv / (double)nval;
        const double se = d0 * std::sqrt(std::max(0.0, sv2 / (double)nval
                                        - (sv / (double)nval) * (sv / (double)nval))
                                        / (double)nval);
        // t = 0 exercise is always available
        const double ex0 = is_put ? std::max(c.K - c.S0, 0.0) : std::max(c.S0 - c.K, 0.0);
        if (is_put) { o.put = std::max(mv, ex0); o.put_se = se; o.put_insample = std::max(mr, ex0); }
        else        { o.call = std::max(mv, ex0); o.call_se = se; o.call_insample = std::max(mr, ex0); }
    }
    o.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return o;
}

// ------------------------------------------------------------------------ main
static std::vector<std::string> split(const std::string& s, char d) {
    std::vector<std::string> out;
    size_t a = 0;
    while (true) {
        size_t b = s.find(d, a);
        out.push_back(s.substr(a, b == std::string::npos ? std::string::npos : b - a));
        if (b == std::string::npos) break;
        a = b + 1;
    }
    return out;
}

int main() {
    std::printf("id,eu_put,eu_put_se,eu_call,eu_call_se,am_put,am_put_se,"
                "am_call,am_call_se,am_put_insample,am_call_insample,"
                "mean_ST,mean_VT,neg_hits,ex_dates,t_eu,t_am\n");
    std::fflush(stdout);
    char line[4096];
    bool header = true;
    while (std::fgets(line, sizeof line, stdin)) {
        std::string L(line);
        while (!L.empty() && (L.back() == '\n' || L.back() == '\r')) L.pop_back();
        if (L.empty()) continue;
        if (header) { header = false; if (L.find("id") == 0) continue; }
        auto f = split(L, ',');
        if (f.size() < 16) { std::fprintf(stderr, "bad line: %s\n", L.c_str()); continue; }
        Cfg c;
        int i = 0;
        c.id = f[i++];
        c.H = std::stod(f[i++]); c.V0 = std::stod(f[i++]); c.theta = std::stod(f[i++]);
        c.kappa = std::stod(f[i++]); c.eta = std::stod(f[i++]); c.rho = std::stod(f[i++]);
        c.r = std::stod(f[i++]); c.T = std::stod(f[i++]); c.S0 = std::stod(f[i++]);
        c.K = std::stod(f[i++]); c.steps = std::stoi(f[i++]);
        c.paths_eu = std::stol(f[i++]); c.paths_reg = std::stol(f[i++]);
        c.paths_val = std::stol(f[i++]); c.ex_stride = std::stoi(f[i++]);
        c.seed = (i < (int)f.size()) ? std::stoull(f[i++]) : 12345ULL;

        EuOut e = run_european(c);
        AmOut a = (c.paths_reg > 0 && c.paths_val > 0) ? run_american(c) : AmOut{};
        std::printf("%s,%.10f,%.3e,%.10f,%.3e,%.10f,%.3e,%.10f,%.3e,"
                    "%.10f,%.10f,%.6f,%.8f,%lld,%d,%.4f,%.4f\n",
                    c.id.c_str(), e.put, e.put_se, e.call, e.call_se,
                    a.put, a.put_se, a.call, a.call_se,
                    a.put_insample, a.call_insample,
                    e.mean_ST, e.mean_VT, (long long)e.neg, a.dates, e.secs, a.secs);
        std::fflush(stdout);
        std::fprintf(stderr, "%s  eu_put=%.4f  am_put=%.4f  t=%.2fs+%.2fs\n",
                     c.id.c_str(), e.put, a.put, e.secs, a.secs);
    }
    return 0;
}
