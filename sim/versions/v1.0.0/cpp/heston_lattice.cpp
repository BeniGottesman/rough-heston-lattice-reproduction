// heston_lattice.cpp — the PAPER'S lattice at h = 0, in C++, multithreaded.
//
// This is a straight port of `sim/heston_lattice.py`, not a new scheme.  Every
// formula, every clip and every index range is the same, so that the two can be
// required to agree to ~1e-12 before this file is used for anything.  The port
// exists for one reason: the Python version costs O(n^4) -- kmax grows like
// mref ~ sqrt(n) and the grid like kmax * n -- so a single price at n = 200 takes
// about seven minutes, and a table of 45 parameter sets at three step counts is
// out of reach.  Here the backward induction is threaded over the driver states
// and the inner loop over price nodes is a contiguous axpy.
//
// The construction, in one paragraph.  Classical Heston is inside the paper's
// model class because at h = 0 the operator K^0 is the identity, so
// v = v_0 + y with y autonomous.  A recombining walk needs a unit diffusion
// coefficient, so the variance goes through the Lamperti transform
// U = 2 sqrt(v) / nu, whose drift is (2 lambda theta / nu^2 - 1/2)/U - lambda U/2
// and whose diffusion is 1.  U therefore lives on u_0 + sqrt(delta) Z and its
// state after k steps is the number of up-moves; the drift is carried by the
// up-probability p = (1 + mu_U sqrt(delta))/2.  The price is coupled to that
// driver by the randomised rounding of Route A', on the finer grid a_X, which is
// what removes the admissibility constraint of the four-point kernel.
//
// All four vanilla payoffs are produced in ONE backward pass: European put,
// European call, American put, American call.  They share the lattice geometry
// and the transition kernel, which is the expensive part, so four prices cost far
// less than four runs.  The American call is a control rather than a product:
// with r >= 0 and no dividend its true value is the European call, so any gap is
// a defect of the induction, visible for free on every line of the table.
//
// Build:  clang++ -O3 -march=native -std=c++17 -pthread
// Usage:  heston_lattice <config.csv >out.csv    (see run_heston_ads_tables.py)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>

// Gauss--Hermite nodes and weights for the standard normal, 25 points, taken
// verbatim from numpy.polynomial.hermite_e.hermegauss(25) with the weights
// normalised to sum to one -- the same constants `route_aprime.interp_kernel`
// uses, so the quadrature cannot differ between the two implementations.
static const int GH = 25;
static const double GX[25] = {-8.71759767839958677e+00, -7.65603795539307441e+00, -6.76746496380971685e+00, -5.96601469060670198e+00, -5.21884809364427937e+00, -4.50892992296728501e+00, -3.82590056997249128e+00, -3.16277567938819271e+00, -2.51447330395220581e+00, -1.87705836994783870e+00, -1.24731197561678941e+00, -6.22462279186075995e-01, 0.00000000000000000e+00, 6.22462279186075995e-01, 1.24731197561678941e+00, 1.87705836994783870e+00, 2.51447330395220581e+00, 3.16277567938819271e+00, 3.82590056997249128e+00, 4.50892992296728501e+00, 5.21884809364427937e+00, 5.96601469060670198e+00, 6.76746496380971685e+00, 7.65603795539307441e+00, 8.71759767839958677e+00};
static const double GW[25] = {1.53003899799871729e-17, 7.10210303700400479e-14, 3.79115000047715087e-11, 5.73802386889937796e-09, 3.53015256024549414e-07, 1.06721949052025481e-05, 1.77766906926526275e-04, 1.75785040526379131e-03, 1.08567559914623038e-02, 4.33799701676450267e-02, 1.14880924303951568e-01, 2.04851025650340440e-01, 2.48169351176485448e-01, 2.04851025650340440e-01, 1.14880924303951568e-01, 4.33799701676450267e-02, 1.08567559914623038e-02, 1.75785040526379131e-03, 1.77766906926526275e-04, 1.06721949052025481e-05, 3.53015256024549414e-07, 5.73802386889937796e-09, 3.79115000047715087e-11, 7.10210303700400479e-14, 1.53003899799871729e-17};

struct Cfg {
    std::string id;
    int n;
    double V0, theta, lam, nu, rho, T, S0, K, r;
    int mref;
    double barrier_sd;      // <= 0 means "no barrier"
    int drift_floor;        // 0 = clip the probabilities, 1 = floor the band
    int walk;               // 0 = binomial +-sqrt(delta), 1 = trinomial Hull--White
};

// ---------------------------------------------------------------------------
// The two driver walks.
//
// BINOMIAL (walk = 0), the construction of the paper as written.  dU = +-sqrt(d)
// with p = (1 + mu sqrt(d))/2.  The mean is exact, but
//
//     Var[dU] = d - mu^2 d^2,
//
// so the walk loses the RELATIVE fraction mu^2 d of the driver's variance.  That
// is harmless when the Lamperti drift is small and ruinous when it is not: on the
// Beliaeva--Nawalkha grid with sqrt(V0) = 0.4 and theta = 0.04 one has mu = -9.06,
// so at n = 200 and T = 1/2 the walk throws away 20.5 per cent of the variance and
// the put comes out 0.26 too low.  Nothing signals it, because |mu| sqrt(d) = 0.45
// keeps the probability inside [0,1].
//
// TRINOMIAL (walk = 1), Hull--White with branch switching.  Nodes are spaced
// du = sqrt(3d) and, from a node u, the three successors are centred on the node
// NEAREST to the conditional mean u + mu d rather than on u itself -- that is the
// branch switching, and it is what lets an arbitrarily large drift be represented.
// Writing e for the residual (mean minus the central node) and A = e/du, B = d/du^2
// = 1/3, the probabilities
//
//     p_u = (B + A^2 + A)/2,   p_d = (B + A^2 - A)/2,   p_m = 1 - B - A^2
//
// match the mean AND the variance exactly.  Since branch switching gives
// |A| <= 1/2, one has p_m >= 2/3 - 1/4 and p_u, p_d >= (1/3 - 1/4)/2, so all three
// are strictly positive for ANY drift: the F048 failure mode cannot occur.
// ---------------------------------------------------------------------------

static inline double lamperti(double v, double nu) {
    return 2.0 * std::sqrt(std::max(v, 0.0)) / nu;
}
static inline double inv_lamperti(double u, double nu) {
    return 0.25 * nu * nu * u * u;
}
static inline double drift_U(double u, double lam, double theta, double nu) {
    const double c = 2.0 * lam * theta / (nu * nu) - 0.5;
    return c / std::max(u, 1e-12) - 0.5 * lam * u;
}

static void cir_moments(double V0, double theta, double lam, double nu, double t,
                        double& mean, double& sd) {
    const double e = std::exp(-lam * t);
    mean = theta + (V0 - theta) * e;
    const double var = V0 * (nu * nu / lam) * (e - e * e)
                     + theta * (nu * nu / (2.0 * lam)) * (1.0 - e) * (1.0 - e);
    sd = std::sqrt(std::max(var, 0.0));
}

// The band is taken over the WHOLE horizon, not at the final date alone: with
// V0 far from theta the mean reverts away from V0, so a band built around E[v_T]
// can fail to contain v_0 itself.  Scanning s in [0,T] and including s = 0 makes
// v_0 an endpoint by construction.  See the docstring of the Python twin.
static const int BARRIER_SCAN = 64;

static void barrier_U(const Cfg& c, double delta, bool use_delta,
                      double& u_lo, double& u_hi) {
    double v_lo = 1e300, v_hi = -1e300;
    for (int i = 0; i <= BARRIER_SCAN; ++i) {
        double mean, sd;
        cir_moments(c.V0, c.theta, c.lam, c.nu,
                    c.T * (double)i / (double)BARRIER_SCAN, mean, sd);
        v_lo = std::min(v_lo, mean - c.barrier_sd * sd);
        v_hi = std::max(v_hi, mean + c.barrier_sd * sd);
    }
    v_lo = std::max(0.0, v_lo);
    u_lo = lamperti(v_lo, c.nu);
    u_hi = lamperti(v_hi, c.nu);
    if (use_delta) {
        const double cc = 2.0 * c.lam * c.theta / (c.nu * c.nu) - 0.5;
        u_lo = std::max(u_lo, 2.0 * std::fabs(cc) * std::sqrt(delta));
    }
    u_lo = std::max(u_lo, 1e-6);
    u_hi = std::min(u_hi, 1e6);
}

static int n_threads() {
    unsigned h = std::thread::hardware_concurrency();
    if (h == 0) h = 4;
    if (h > 2) h -= 1;
    return (int)h;
}

struct Out {
    double eu_put = 0, eu_call = 0, am_put = 0, am_call = 0;
    long grid = 0, offsets = 0, driver_states = 0;
    long long violations = 0;
    double secs = 0;
    double v_at_upper = 0, feller = 0;
    // measured at every node of every step: how far the walk's own first two
    // moments are from the ones it is supposed to reproduce.  This is the direct
    // check on the mechanism, not an inference from prices.
    double max_var_err = 0;      // max |Var[dU]/d - 1|
    double max_mean_err = 0;     // max |E[dU] - mu d| / (|mu| d + d)
    double max_abs_drift = 0;    // max |mu| over the band, for the record
};

// The four payoffs carried through one backward pass.  Index into VAL:
//   0 European put, 1 European call, 2 American put, 3 American call
static const int NV = 4;
static const bool IS_AMERICAN[NV] = {false, false, true, true};
static const bool IS_PUT[NV] = {true, false, true, false};

static Out price(const Cfg& c) {
    const auto t0 = std::chrono::steady_clock::now();
    const double d = c.T / c.n;
    const double sqd = std::sqrt(d);
    const double u0 = lamperti(c.V0, c.nu);

    double u_lo, u_hi;
    if (c.barrier_sd <= 0.0) {
        u_lo = 1e-3;
        u_hi = u0 + sqd * c.n;
    } else {
        barrier_U(c, d, c.drift_floor != 0, u_lo, u_hi);
    }
    const double v_max = inv_lamperti(std::min(u_hi, u0 + sqd * c.n), c.nu);
    const double s_ref = std::sqrt(c.V0);
    const double a_X = s_ref * sqd / c.mref;
    const double sig_max = std::sqrt(std::max(v_max, c.V0));

    // trinomial geometry: a FIXED node grid u0 + j du over the absorbing band
    const bool tri = (c.walk == 1);
    const double du = tri ? std::sqrt(3.0 * d) : 0.0;
    int j_lo = 0, j_hi = 0, J = 0;
    if (tri) {
        j_lo = (int)std::ceil((u_lo - u0) / du);
        j_hi = (int)std::floor((u_hi - u0) / du);
        if (j_lo > 0) j_lo = 0;            // u0 is a node and must be inside
        if (j_hi < 0) j_hi = 0;
        if (j_hi - j_lo < 2) { j_lo = std::min(j_lo, -1); j_hi = std::max(j_hi, 1); }
        J = j_hi - j_lo + 1;
    }
    // the widest centred driver increment the price has to follow: sqrt(d)|zeta|
    // for the binomial, and 1.5 du = 1.5 sqrt(3) sqrt(d) for the trinomial, since
    // branch switching leaves |e| <= du/2 on top of the +-du branch
    const double zc_max = tri ? 1.5 * std::sqrt(3.0) : 1.0;
    const int kmax = (int)std::ceil(4.5 * sig_max * std::sqrt(1.0 - c.rho * c.rho)
                                    * sqd / a_X
                                    + std::fabs(c.rho) * sig_max * zc_max * sqd / a_X)
                     + 1;
    const long reach = (long)kmax * c.n;
    const long nx = 2 * reach + 1;
    const long ix0 = reach;

    std::vector<double> pay[2];              // 0 = put, 1 = call
    pay[0].resize(nx); pay[1].resize(nx);
    for (long i = 0; i < nx; ++i) {
        const double s = std::exp(std::log(c.S0) + (double)(i - ix0) * a_X);
        pay[0][i] = std::max(c.K - s, 0.0);
        pay[1][i] = std::max(s - c.K, 0.0);
    }

    // val[v][s * nx + ix]; the driver state s runs 0..kk at step kk for the
    // binomial (number of up-moves) and 0..J-1 for the trinomial (node index)
    const int max_states = tri ? J : (c.n + 1);
    std::vector<std::vector<double>> val(NV);
    for (int v = 0; v < NV; ++v) {
        val[v].resize((size_t)nx * max_states);
        for (int m = 0; m < max_states; ++m)
            std::copy(pay[IS_PUT[v] ? 0 : 1].begin(), pay[IS_PUT[v] ? 0 : 1].end(),
                      val[v].begin() + (size_t)m * nx);
    }
    std::vector<std::vector<double>> nxt(NV);
    for (int v = 0; v < NV; ++v) nxt[v].resize((size_t)nx * max_states);

    const double disc = std::exp(-c.r * d);
    const int nt = n_threads();
    const int noff = 2 * kmax + 1;
    long long violations = 0;
    double o_max_var_err = 0.0, o_max_mean_err = 0.0, o_max_drift = 0.0;

    const int NBR = tri ? 3 : 2;             // branches per driver state
    for (int kk = c.n - 1; kk >= 0; --kk) {
        const long lo = std::max(0L, ix0 - (long)kmax * kk);
        const long hi = std::min(nx, ix0 + (long)kmax * kk + 1);
        const long W = hi - lo;
        const int S = tri ? J : (kk + 1);    // driver states at this step

        // per-driver-state quantities and the branch list
        std::vector<double> sig(S), vh(S), mu(S);
        std::vector<double> bw((size_t)S * NBR), bdb((size_t)S * NBR);
        std::vector<int> bns((size_t)S * NBR);
        for (int m = 0; m < S; ++m) {
            double u = tri ? (u0 + (double)(j_lo + m) * du)
                           : (u0 + sqd * (2.0 * m - kk));
            u = std::min(std::max(u, u_lo), u_hi);
            vh[m] = inv_lamperti(u, c.nu);
            sig[m] = std::sqrt(vh[m]);
            mu[m] = drift_U(u, c.lam, c.theta, c.nu);
            if (std::fabs(mu[m]) > o_max_drift) o_max_drift = std::fabs(mu[m]);
            if (!tri) {
                double p = 0.5 * (1.0 + mu[m] * sqd);
                if (p < 0.0 || p > 1.0) ++violations;
                p = std::min(std::max(p, 0.0), 1.0);
                // order kept as (down, up), as in the version already validated
                bw[(size_t)m * 2 + 0] = 1.0 - p;
                bns[(size_t)m * 2 + 0] = m;
                bdb[(size_t)m * 2 + 0] = -sqd - mu[m] * d;
                bw[(size_t)m * 2 + 1] = p;
                bns[(size_t)m * 2 + 1] = m + 1;
                bdb[(size_t)m * 2 + 1] = sqd - mu[m] * d;
            } else {
                // branch switching: centre on the node nearest the conditional
                // mean, clamped so that both neighbours stay inside the band
                const double mean = u + mu[m] * d;
                int js = (int)std::lround((mean - u0) / du);
                if (js < j_lo + 1) js = j_lo + 1;
                if (js > j_hi - 1) js = j_hi - 1;
                const double e = mean - (u0 + (double)js * du);
                const double A = e / du;
                double p[3];
                p[0] = (1.0 / 3.0 + A * A - A) / 2.0;      // i = -1
                p[1] = 2.0 / 3.0 - A * A;                  // i =  0
                p[2] = (1.0 / 3.0 + A * A + A) / 2.0;      // i = +1
                bool bad = false;
                for (int i = 0; i < 3; ++i) if (p[i] < 0.0 || p[i] > 1.0) bad = true;
                if (bad) {
                    ++violations;
                    double s = 0.0;
                    for (int i = 0; i < 3; ++i) {
                        p[i] = std::min(std::max(p[i], 0.0), 1.0);
                        s += p[i];
                    }
                    for (int i = 0; i < 3; ++i) p[i] /= s;
                }
                for (int i = 0; i < 3; ++i) {
                    bw[(size_t)m * 3 + i] = p[i];
                    bns[(size_t)m * 3 + i] = (js + (i - 1)) - j_lo;
                    bdb[(size_t)m * 3 + i] = (double)(i - 1) * du - e;
                }
            }
            // the direct check on the walk's own first two moments
            double m1 = 0.0, m2 = 0.0;
            for (int b = 0; b < NBR; ++b) {
                const double w = bw[(size_t)m * NBR + b];
                const double x = bdb[(size_t)m * NBR + b];
                m1 += w * x;
                m2 += w * x * x;
            }
            const double var = m2 - m1 * m1;
            const double ve = std::fabs(var / d - 1.0);
            const double me = std::fabs(m1) / (std::fabs(mu[m]) * d + d);
            if (ve > o_max_var_err) o_max_var_err = ve;
            if (me > o_max_mean_err) o_max_mean_err = me;
        }

        // transition kernel P[branch][offset][state], by Gauss--Hermite quadrature
        std::vector<double> P((size_t)NBR * noff * S, 0.0);
        for (int b = 0; b < NBR; ++b) {
            for (int m = 0; m < S; ++m) {
                const double mu_M = sig[m] * c.rho * bdb[(size_t)m * NBR + b]
                                  - 0.5 * vh[m] * d + c.r * d;
                const double sd_M = sig[m] * std::sqrt(1.0 - c.rho * c.rho) * sqd;
                for (int gq = 0; gq < GH; ++gq) {
                    const double z = (mu_M + sd_M * GX[gq]) / a_X;
                    // Lambda(z - o) is non-zero only for the two o straddling z
                    const int olo = (int)std::floor(z) - 1, ohi = olo + 3;
                    for (int o = olo; o <= ohi; ++o) {
                        if (o < -kmax || o > kmax) continue;
                        const double t = 1.0 - std::fabs(z - (double)o);
                        if (t <= 0.0) continue;
                        P[((size_t)b * noff + (o + kmax)) * S + m] += GW[gq] * t;
                    }
                }
            }
        }

        // backward step, threaded over driver states
        auto work = [&](int tid) {
            for (int m = tid; m < S; m += nt) {
                double* dst[NV];
                for (int v = 0; v < NV; ++v)
                    dst[v] = nxt[v].data() + (size_t)m * nx + lo;
                for (int v = 0; v < NV; ++v) std::fill(dst[v], dst[v] + W, 0.0);
                for (int b = 0; b < NBR; ++b) {
                    const double w = bw[(size_t)m * NBR + b];
                    if (w == 0.0) continue;
                    const int mn = bns[(size_t)m * NBR + b];
                    for (int oi = 0; oi < noff; ++oi) {
                        const double pk = P[((size_t)b * noff + oi) * S + m];
                        if (pk == 0.0) continue;
                        const double wp = w * pk;
                        const long o = oi - kmax;
                        for (int v = 0; v < NV; ++v) {
                            const double* src =
                                val[v].data() + (size_t)mn * nx;
                            double* dd = dst[v];
                            // src index is clipped to the grid, as in the Python
                            for (long i = 0; i < W; ++i) {
                                long s = lo + i + o;
                                if (s < 0) s = 0;
                                else if (s >= nx) s = nx - 1;
                                dd[i] += wp * src[s];
                            }
                        }
                    }
                }
                for (int v = 0; v < NV; ++v) {
                    double* dd = dst[v];
                    const double* pp = pay[IS_PUT[v] ? 0 : 1].data() + lo;
                    for (long i = 0; i < W; ++i) {
                        double x = disc * dd[i];
                        if (IS_AMERICAN[v] && pp[i] > x) x = pp[i];
                        dd[i] = x;
                    }
                    // outside the reachable band the value is carried over
                    double* col = nxt[v].data() + (size_t)m * nx;
                    const double* old = val[v].data() + (size_t)m * nx;
                    for (long i = 0; i < lo; ++i) col[i] = old[i];
                    for (long i = hi; i < nx; ++i) col[i] = old[i];
                }
            }
        };
        std::vector<std::thread> th;
        for (int t = 0; t < nt; ++t) th.emplace_back(work, t);
        for (auto& t : th) t.join();
        for (int v = 0; v < NV; ++v) val[v].swap(nxt[v]);
    }

    Out o;
    const size_t root = (size_t)(tri ? -j_lo : 0) * nx + ix0;
    o.eu_put = val[0][root];
    o.eu_call = val[1][root];
    o.am_put = val[2][root];
    o.am_call = val[3][root];
    o.grid = nx;
    o.driver_states = tri ? J : (c.n + 1);
    o.max_var_err = o_max_var_err;
    o.max_mean_err = o_max_mean_err;
    o.max_abs_drift = o_max_drift;
    o.offsets = noff;
    o.violations = violations;
    o.v_at_upper = inv_lamperti(u_hi, c.nu);
    o.feller = 2.0 * c.lam * c.theta / (c.nu * c.nu);
    o.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
                 .count();
    return o;
}

static std::vector<std::string> split(const std::string& s, char dl) {
    std::vector<std::string> out;
    size_t a = 0;
    while (true) {
        size_t b = s.find(dl, a);
        out.push_back(s.substr(a, b == std::string::npos ? std::string::npos : b - a));
        if (b == std::string::npos) break;
        a = b + 1;
    }
    return out;
}

int main() {
    std::printf("id,walk,n,eu_put,eu_call,am_put,am_call,grid,offsets,"
                "driver_states,violations,max_var_err,max_mean_err,"
                "max_abs_drift,v_at_upper,feller,secs\n");
    std::fflush(stdout);
    char line[4096];
    bool header = true;
    while (std::fgets(line, sizeof line, stdin)) {
        std::string L(line);
        while (!L.empty() && (L.back() == '\n' || L.back() == '\r')) L.pop_back();
        if (L.empty()) continue;
        if (header) { header = false; if (L.find("id") == 0) continue; }
        auto f = split(L, ',');
        if (f.size() < 15) { std::fprintf(stderr, "bad line: %s\n", L.c_str()); continue; }
        Cfg c;
        int i = 0;
        c.id = f[i++];
        c.n = std::stoi(f[i++]);
        c.V0 = std::stod(f[i++]); c.theta = std::stod(f[i++]);
        c.lam = std::stod(f[i++]); c.nu = std::stod(f[i++]);
        c.rho = std::stod(f[i++]); c.T = std::stod(f[i++]);
        c.S0 = std::stod(f[i++]); c.K = std::stod(f[i++]);
        c.r = std::stod(f[i++]); c.mref = std::stoi(f[i++]);
        c.barrier_sd = std::stod(f[i++]); c.drift_floor = std::stoi(f[i++]);
        c.walk = std::stoi(f[i++]);
        Out o = price(c);
        std::printf("%s,%d,%d,%.12f,%.12f,%.12f,%.12f,%ld,%ld,%ld,%lld,"
                    "%.3e,%.3e,%.4f,%.8f,%.4f,%.4f\n",
                    c.id.c_str(), c.walk, c.n, o.eu_put, o.eu_call, o.am_put,
                    o.am_call, o.grid, o.offsets, o.driver_states, o.violations,
                    o.max_var_err, o.max_mean_err, o.max_abs_drift,
                    o.v_at_upper, o.feller, o.secs);
        std::fflush(stdout);
        std::fprintf(stderr, "%s walk=%d n=%d eu_put=%.6f am_put=%.6f "
                     "var_err=%.2e %.2fs\n", c.id.c_str(), c.walk, c.n,
                     o.eu_put, o.am_put, o.max_var_err, o.secs);
    }
    return 0;
}
