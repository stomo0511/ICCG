//
// Estimate the spectral condition number of an SPD CRS matrix from CG
// coefficients via the associated Lanczos tridiagonal matrix.
//
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include "crs_io.hpp"

#ifdef USEMKL
#include <mkl.h>
#endif

// Inner product
static inline double dot(const std::vector<double>& a, const std::vector<double>& b) {
    long double s = 0.0L;
    const size_t n = a.size();
    #pragma omp parallel for reduction(+:s) schedule(static)
    for (size_t i = 0; i < n; ++i)
        s += (long double)a[i] * b[i];
    return (double)s;
}

// 2-norm
static inline double nrm2(const std::vector<double>& a) {
    long double s = 0.0L;
    const size_t n = a.size();
    #pragma omp parallel for reduction(+:s) schedule(static)
    for (size_t i = 0; i < n; ++i)
        s += (long double)a[i] * a[i];
    return (double)std::sqrt((double)s);
}

// Sparse matrix-vector product: y = A*x
static void spmv(const CRS& A, const std::vector<double>& x, std::vector<double>& y) {
    int n = A.n;
    if ((int)y.size() != n) y.resize(n);

    #ifdef USEMKL
    const double alpha = 1.0, beta = 0.0;
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A.handle, A.descr, x.data(), beta, y.data());
    #else
    const int*    rowptr = A.rowptr.data();
    const int*    colind = A.colind.data();
    const double* aval   = A.val.data();
    const double* xp     = x.data();
    double*       yp     = y.data();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        const int rs = rowptr[i];
        const int re = rowptr[i+1];
        for (int k = rs; k < re; ++k)
            sum += aval[k] * xp[colind[k]];
        yp[i] = sum;
    }
    #endif
}

struct CGResult {
    int iters = 0;
    double rel_resid = NAN;
    bool converged = false;
    bool breakdown = false;
    std::string message;
};

struct CondEstimate {
    double lambda_min = NAN;
    double lambda_max = NAN;
    double cond = NAN;
};

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " matrix.mtx [tol=1e-08] [max_iter=10000] [interval=100]\n";
}

static bool parse_double_arg(const char* s, double& v) {
    char* end = nullptr;
    errno = 0;
    v = std::strtod(s, &end);
    return errno == 0 && end != s && *end == '\0' && std::isfinite(v);
}

static bool parse_int_arg(const char* s, int& v) {
    char* end = nullptr;
    errno = 0;
    long x = std::strtol(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0' || x > std::numeric_limits<int>::max())
        return false;
    v = (int)x;
    return true;
}

static void build_tridiagonal(
    const std::vector<double>& alpha_hist,
    const std::vector<double>& beta_hist,
    int m,
    std::vector<double>& d,
    std::vector<double>& e
) {
    if (m <= 0 || m > (int)alpha_hist.size())
        throw std::runtime_error("invalid tridiagonal size");
    if ((int)beta_hist.size() < m - 1)
        throw std::runtime_error("not enough beta coefficients for tridiagonal construction");

    d.assign(m, 0.0);
    e.assign(std::max(0, m - 1), 0.0);

    for (int k = 0; k < m; ++k) {
        const double alpha = alpha_hist[k];
        if (!(std::isfinite(alpha) && alpha > 0.0))
            throw std::runtime_error("invalid CG coefficient: alpha must be positive and finite");
    }
    for (int k = 0; k < m - 1; ++k) {
        const double beta = beta_hist[k];
        if (!(std::isfinite(beta) && beta >= 0.0))
            throw std::runtime_error("invalid CG coefficient: beta must be nonnegative and finite");
    }

    d[0] = 1.0 / alpha_hist[0];
    for (int k = 1; k < m; ++k)
        d[k] = 1.0 / alpha_hist[k] + beta_hist[k - 1] / alpha_hist[k - 1];
    for (int k = 0; k < m - 1; ++k)
        e[k] = std::sqrt(beta_hist[k]) / alpha_hist[k];

    for (double v : d) {
        if (!std::isfinite(v))
            throw std::runtime_error("invalid tridiagonal diagonal value");
    }
    for (double v : e) {
        if (!std::isfinite(v))
            throw std::runtime_error("invalid tridiagonal off-diagonal value");
    }
}

static int sturm_count_less_equal(
    const std::vector<double>& d,
    const std::vector<double>& e,
    double x
) {
    const int n = (int)d.size();
    int count = 0;
    const double tiny = std::numeric_limits<double>::min();

    double q = d[0] - x;
    if (q <= 0.0) ++count;
    if (std::abs(q) < tiny) q = (q < 0.0) ? -tiny : tiny;

    for (int i = 1; i < n; ++i) {
        q = d[i] - x - (e[i - 1] * e[i - 1]) / q;
        if (q <= 0.0) ++count;
        if (std::abs(q) < tiny) q = (q < 0.0) ? -tiny : tiny;
    }
    return count;
}

static void gershgorin_bounds(
    const std::vector<double>& d,
    const std::vector<double>& e,
    double& lower,
    double& upper
) {
    lower = std::numeric_limits<double>::infinity();
    upper = -std::numeric_limits<double>::infinity();
    const int n = (int)d.size();
    for (int i = 0; i < n; ++i) {
        double radius = 0.0;
        if (i > 0) radius += std::abs(e[i - 1]);
        if (i + 1 < n) radius += std::abs(e[i]);
        lower = std::min(lower, d[i] - radius);
        upper = std::max(upper, d[i] + radius);
    }
    if (!(std::isfinite(lower) && std::isfinite(upper) && lower < upper))
        throw std::runtime_error("failed to bracket tridiagonal eigenvalues");
}

static double tridiagonal_eigenvalue_by_index(
    const std::vector<double>& d,
    const std::vector<double>& e,
    int index
) {
    double lo = 0.0, hi = 0.0;
    gershgorin_bounds(d, e, lo, hi);

    const int n = (int)d.size();
    if (index < 0 || index >= n)
        throw std::runtime_error("invalid eigenvalue index");

    const int target = index + 1;
    const double span = hi - lo;
    lo -= 8.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(lo));
    hi += 8.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(hi));

    for (int iter = 0; iter < 100; ++iter) {
        double mid = 0.5 * (lo + hi);
        if (mid == lo || mid == hi)
            break;
        int count = sturm_count_less_equal(d, e, mid);
        if (count >= target)
            hi = mid;
        else
            lo = mid;
        if (hi - lo <= 16.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, span))
            break;
    }
    return 0.5 * (lo + hi);
}

static CondEstimate estimate_condition_from_prefix(
    const std::vector<double>& alpha_hist,
    const std::vector<double>& beta_hist,
    int m
) {
    std::vector<double> d, e;
    build_tridiagonal(alpha_hist, beta_hist, m, d, e);

    CondEstimate est;
    est.lambda_min = tridiagonal_eigenvalue_by_index(d, e, 0);
    est.lambda_max = tridiagonal_eigenvalue_by_index(d, e, m - 1);

    if (!(std::isfinite(est.lambda_min) && std::isfinite(est.lambda_max)))
        throw std::runtime_error("tridiagonal eigenvalue computation produced a non-finite value");
    if (est.lambda_min <= 0.0)
        throw std::runtime_error("estimated lambda_min is nonpositive; condition number was not computed");
    if (est.lambda_max < est.lambda_min)
        throw std::runtime_error("estimated lambda_max is smaller than lambda_min");

    // cond_est is an estimate based on extremal Ritz values obtained from the CG coefficients.
    est.cond = est.lambda_max / est.lambda_min;
    if (!std::isfinite(est.cond))
        throw std::runtime_error("condition estimate is not finite");
    return est;
}

static void print_checkpoint(
    int iteration,
    const std::vector<double>& alpha_hist,
    const std::vector<double>& beta_hist
) {
    CondEstimate est = estimate_condition_from_prefix(alpha_hist, beta_hist, iteration);
    std::cout << iteration << ','
              << est.lambda_min << ','
              << est.lambda_max << ','
              << est.cond << '\n';
}

static CGResult conjugate_gradient_no_precond(
    const CRS& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    int max_iter,
    double tol,
    int interval,
    std::vector<double>& alpha_hist,
    std::vector<double>& beta_hist
) {
    const int n = A.n;
    std::vector<double> r(n), p(n), Ap(n), Ax(n);

    spmv(A, x, Ax);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
        r[i] = b[i] - Ax[i];

    double normb = nrm2(b);
    if (normb == 0.0)
        normb = 1.0;

    p = r;
    double rr_old = dot(r, r);

    CGResult res;
    res.rel_resid = nrm2(r) / normb;
    if (res.rel_resid < tol) {
        res.converged = true;
        return res;
    }

    for (int k = 0; k < max_iter; ++k) {
        spmv(A, p, Ap);
        double pAp = dot(p, Ap);
        double eps = std::numeric_limits<double>::epsilon();
        double thr = eps * nrm2(p) * nrm2(Ap);
        if (!(std::isfinite(pAp)) || pAp <= thr) {
            res.iters = k;
            res.breakdown = true;
            res.message = "CG breakdown: nonpositive or invalid pAp";
            return res;
        }

        double alpha = rr_old / pAp;
        if (!(std::isfinite(alpha) && alpha > 0.0)) {
            res.iters = k;
            res.breakdown = true;
            res.message = "CG breakdown: invalid alpha";
            return res;
        }
        alpha_hist.push_back(alpha);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rel = nrm2(r) / normb;
        res.iters = k + 1;
        res.rel_resid = rel;

        if (res.iters % interval == 0)
            print_checkpoint(res.iters, alpha_hist, beta_hist);

        if (rel < tol) {
            res.converged = true;
            return res;
        }
        if (res.iters == max_iter)
            return res;

        double rr_new = dot(r, r);
        double beta = rr_new / rr_old;
        if (!(std::isfinite(beta) && beta >= 0.0)) {
            res.breakdown = true;
            res.message = "CG breakdown: invalid beta";
            return res;
        }
        beta_hist.push_back(beta);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i)
            p[i] = r[i] + beta * p[i];

        rr_old = rr_new;
    }

    return res;
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 5) {
        usage(argv[0]);
        return 1;
    }

    std::string path = argv[1];
    double tol = 1e-8;
    int max_iter = 10000;
    int interval = 100;

    if (argc >= 3 && !parse_double_arg(argv[2], tol)) {
        usage(argv[0]);
        return 1;
    }
    if (argc >= 4 && !parse_int_arg(argv[3], max_iter)) {
        usage(argv[0]);
        return 1;
    }
    if (argc >= 5 && !parse_int_arg(argv[4], interval)) {
        usage(argv[0]);
        return 1;
    }
    if (!(tol > 0.0) || max_iter <= 0 || interval <= 0) {
        usage(argv[0]);
        return 1;
    }

    try {
        CRS A = read_matrix_market_crs(path);
        #ifdef USEMKL
        A.build_mkl_handle();
        #endif

        if (A.n <= 0)
            throw std::runtime_error("invalid matrix size");

        std::vector<double> x(A.n, 0.0);
        // RHS is the all-one vector for consistency with the existing CG experiments.
        // Extremal Ritz estimates can be poor if the initial residual is nearly orthogonal
        // to an eigenvector associated with an extremal eigenvalue.
        std::vector<double> b(A.n, 1.0);

        std::vector<double> alpha_hist;
        std::vector<double> beta_hist;
        alpha_hist.reserve(max_iter);
        beta_hist.reserve(std::max(0, max_iter - 1));

        std::cout << "iteration,lambda_min_est,lambda_max_est,cond_est\n";
        CGResult out = conjugate_gradient_no_precond(
            A, b, x, max_iter, tol, interval, alpha_hist, beta_hist);

        if (out.breakdown)
            throw std::runtime_error(out.message);
        if (!out.converged)
            std::cerr << "Warning: CG did not converge within max_iter; estimating from available coefficients.\n";
        if (out.iters > 0 && out.iters < interval)
            std::cerr << "Warning: CG converged before the first checkpoint; condition estimate may be under-resolved.\n";
        if (out.iters > 0 && out.iters < 3)
            std::cerr << "Warning: very few CG iterations were available for condition estimation.\n";
        if (out.iters <= 0)
            throw std::runtime_error("no CG coefficients were generated; condition estimate is unavailable");

        if (out.iters % interval != 0)
            print_checkpoint(out.iters, alpha_hist, beta_hist);

        CondEstimate final_est = estimate_condition_from_prefix(alpha_hist, beta_hist, out.iters);
        std::cout << "FINAL,"
                  << out.iters << ','
                  << out.rel_resid << ','
                  << final_est.lambda_min << ','
                  << final_est.lambda_max << ','
                  << final_est.cond << '\n';
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
