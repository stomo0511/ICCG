#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "crs_io.hpp"
#include "color.hpp"
#include "block.hpp"
#include "ic0.hpp"

struct FrobeniusStats {
    double norm_A = 0.0;
    double norm_LLt = 0.0;
    double inner_A_LLt = 0.0;
    double norm_error = 0.0;
    double epsilon = 0.0;
};

static std::string matrix_name_from_path(const std::string& path) {
    size_t slash = path.find_last_of("/\\");
    std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
    size_t dot = name.find_last_of('.');
    if (dot != std::string::npos) name.erase(dot);
    return name;
}

static long double frobenius_lower_squared(const CRS& A_lower) {
    long double s = 0.0L;
    for (int i = 0; i < A_lower.n; ++i) {
        for (int k = A_lower.rowptr[i]; k < A_lower.rowptr[i + 1]; ++k) {
            const int j = A_lower.colind[k];
            const long double v = A_lower.val[k];
            s += (j == i) ? v * v : 2.0L * v * v;
        }
    }
    return s;
}

static long long symmetric_nnz_from_lower(const CRS& A_lower) {
    long long nnz = 0;
    for (int i = 0; i < A_lower.n; ++i) {
        for (int k = A_lower.rowptr[i]; k < A_lower.rowptr[i + 1]; ++k) {
            nnz += (A_lower.colind[k] == i) ? 1 : 2;
        }
    }
    return nnz;
}

static double row_dot_lower_L(const std::vector<int>& rowptr,
                              const std::vector<int>& colind,
                              const std::vector<double>& val,
                              int i, int j) {
    int pi = rowptr[i];
    int pj = rowptr[j];
    const int ei = rowptr[i + 1];
    const int ej = rowptr[j + 1];
    long double s = 0.0L;

    while (pi < ei && pj < ej) {
        const int ci = colind[pi];
        const int cj = colind[pj];
        if (ci == cj) {
            s += (long double)val[pi] * val[pj];
            ++pi;
            ++pj;
        } else if (ci < cj) {
            ++pi;
        } else {
            ++pj;
        }
    }
    return (double)s;
}

static long double inner_A_LLt_lower(const CRS& A_lower, const IC0& M) {
    const auto& rowptr = M.l_rowptr();
    const auto& colind = M.l_colind();
    const auto& val = M.l_val();

    long double s = 0.0L;
    for (int i = 0; i < A_lower.n; ++i) {
        for (int k = A_lower.rowptr[i]; k < A_lower.rowptr[i + 1]; ++k) {
            const int j = A_lower.colind[k];
            const long double lij = row_dot_lower_L(rowptr, colind, val, i, j);
            const long double aij = A_lower.val[k];
            s += (j == i) ? aij * lij : 2.0L * aij * lij;
        }
    }
    return s;
}

static void build_columns_of_L(const IC0& M,
                               std::vector<int>& colptr,
                               std::vector<int>& rowind,
                               std::vector<double>& cval) {
    const int n = M.size();
    const auto& rowptr = M.l_rowptr();
    const auto& colind = M.l_colind();
    const auto& val = M.l_val();
    const int nnz = (int)colind.size();

    colptr.assign(n + 1, 0);
    rowind.resize(nnz);
    cval.resize(nnz);

    for (int k = 0; k < nnz; ++k) ++colptr[colind[k] + 1];
    for (int c = 0; c < n; ++c) colptr[c + 1] += colptr[c];

    std::vector<int> cur = colptr;
    for (int i = 0; i < n; ++i) {
        for (int k = rowptr[i]; k < rowptr[i + 1]; ++k) {
            const int c = colind[k];
            const int p = cur[c]++;
            rowind[p] = i;
            cval[p] = val[k];
        }
    }
}

static long double frobenius_LLt_squared_exact(const IC0& M) {
    const int n = M.size();
    const auto& L_rowptr = M.l_rowptr();
    const auto& L_colind = M.l_colind();
    const auto& L_val = M.l_val();

    std::vector<int> colptr, rowind;
    std::vector<double> cval;
    build_columns_of_L(M, colptr, rowind, cval);

    std::vector<int> marker(n, -1);
    std::vector<double> acc(n, 0.0);
    std::vector<int> touched;
    long double norm2 = 0.0L;

    // ||LL^T||_F^2 = ||L^T L||_F^2.  For each column j of L, form the
    // nonzeros of column j of G=L^T L in a temporary sparse accumulator.
    // This is exact and includes fill outside Ap, without storing LL^T.
    // Work is O(sum_r nnz(L(r,:))^2); extra memory is O(n + max temporary nnz).
    for (int j = 0; j < n; ++j) {
        touched.clear();
        for (int p = colptr[j]; p < colptr[j + 1]; ++p) {
            const int r = rowind[p];
            const double lrj = cval[p];
            for (int q = L_rowptr[r]; q < L_rowptr[r + 1]; ++q) {
                const int c = L_colind[q];
                if (c > j) break;
                if (marker[c] != j) {
                    marker[c] = j;
                    acc[c] = 0.0;
                    touched.push_back(c);
                }
                acc[c] += lrj * L_val[q];
            }
        }

        for (int c : touched) {
            const long double g = acc[c];
            norm2 += (c == j) ? g * g : 2.0L * g * g;
        }
    }
    return norm2;
}

static FrobeniusStats evaluate_frobenius_error(const CRS& Ap_lower, const IC0& M) {
    FrobeniusStats stats;
    const long double norm_A2 = frobenius_lower_squared(Ap_lower);
    const long double norm_LLt2 = frobenius_LLt_squared_exact(M);
    const long double inner = inner_A_LLt_lower(Ap_lower, M);
    long double err2 = norm_A2 + norm_LLt2 - 2.0L * inner;
    if (err2 < 0.0L && err2 > -100.0L * std::numeric_limits<double>::epsilon() * norm_A2) {
        err2 = 0.0L;
    }
    err2 = std::max(err2, 0.0L);

    stats.norm_A = std::sqrt((double)norm_A2);
    stats.norm_LLt = std::sqrt((double)norm_LLt2);
    stats.inner_A_LLt = (double)inner;
    stats.norm_error = std::sqrt((double)err2);
    stats.epsilon = (norm_A2 > 0.0L) ? (stats.norm_error / stats.norm_A) : std::numeric_limits<double>::quiet_NaN();
    return stats;
}

static std::vector<std::vector<double>> dense_from_lower(const CRS& A_lower) {
    std::vector<std::vector<double>> A(A_lower.n, std::vector<double>(A_lower.n, 0.0));
    for (int i = 0; i < A_lower.n; ++i) {
        for (int k = A_lower.rowptr[i]; k < A_lower.rowptr[i + 1]; ++k) {
            const int j = A_lower.colind[k];
            A[i][j] += A_lower.val[k];
            if (j != i) A[j][i] += A_lower.val[k];
        }
    }
    return A;
}

static double dense_check_epsilon(const CRS& Ap_lower, const IC0& M) {
    const int n = Ap_lower.n;
    auto A = dense_from_lower(Ap_lower);
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

    const auto& rowptr = M.l_rowptr();
    const auto& colind = M.l_colind();
    const auto& val = M.l_val();
    for (int i = 0; i < n; ++i) {
        for (int k = rowptr[i]; k < rowptr[i + 1]; ++k) {
            L[i][colind[k]] = val[k];
        }
    }

    long double norm_A2 = 0.0L;
    long double err2 = 0.0L;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            long double llt = 0.0L;
            const int lim = std::min(i, j);
            for (int k = 0; k <= lim; ++k) llt += (long double)L[i][k] * L[j][k];
            const long double a = A[i][j];
            const long double d = a - llt;
            norm_A2 += a * a;
            err2 += d * d;
        }
    }
    return (norm_A2 > 0.0L) ? std::sqrt((double)(err2 / norm_A2)) : std::numeric_limits<double>::quiet_NaN();
}

static void usage(const char* prog) {
    std::cerr << "usage: " << prog << " [--debug] [--dense-check] matrix.mtx block.blk bcolor.bcol\n";
}

int main(int argc, char** argv) {
    try {
        bool debug = false;
        bool dense_check = false;
        std::vector<std::string> args;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--debug") {
                debug = true;
            } else if (a == "--dense-check") {
                dense_check = true;
            } else {
                args.push_back(a);
            }
        }
        if (args.size() != 3) {
            usage(argv[0]);
            return 1;
        }

        const std::string matrix_file = args[0];
        const std::string block_file = args[1];
        const std::string bcolor_file = args[2];

        CRS A = read_matrix_market_crs(matrix_file);
        int nb = 0;
        std::vector<int> block_of_old = ReadBlockFile_1Based(block_file, A.n, nb);
        int nc = 0;
        std::vector<int> block_color = ReadBlockColorFile_1Based(bcolor_file, nb, nc);

        std::vector<int> new_of_old = BuildPermutationByBlockColor(block_of_old, block_color);
        CRS Ap_lower = Permute_PtAP_LowerCRS_to_LowerCRS(A, new_of_old);

        const double shift = 0.0;
        IC0 M(Ap_lower, shift);
        FrobeniusStats stats = evaluate_frobenius_error(Ap_lower, M);

        if (dense_check) {
            if (Ap_lower.n > 512) {
                throw std::runtime_error("--dense-check is limited to n <= 512");
            }
            const double eps_dense = dense_check_epsilon(Ap_lower, M);
            const double denom = std::max({1.0, std::abs(eps_dense), std::abs(stats.epsilon)});
            const double rel = std::abs(eps_dense - stats.epsilon) / denom;
            if (rel > 1e-10) {
                throw std::runtime_error("dense check failed: sparse epsilon=" +
                                         std::to_string(stats.epsilon) +
                                         " dense epsilon=" + std::to_string(eps_dense));
            }
            if (debug) {
                std::cerr << "dense_epsilon," << std::scientific << std::setprecision(12)
                          << eps_dense << "\n";
            }
        }

        if (debug) {
            std::cerr << std::scientific << std::setprecision(12)
                      << "norm_A_F," << stats.norm_A << "\n"
                      << "norm_LLt_F," << stats.norm_LLt << "\n"
                      << "inner_A_LLt_F," << stats.inner_A_LLt << "\n"
                      << "norm_error_F," << stats.norm_error << "\n"
                      << "epsilon_F," << stats.epsilon << "\n";
        }

        std::cout << std::scientific << std::setprecision(12)
                  << matrix_name_from_path(matrix_file) << ", "
                  << Ap_lower.n << ", "
                  << symmetric_nnz_from_lower(Ap_lower) << ", "
                  << nb << ", "
                  << nc << ", "
                  << stats.norm_A << ", "
                  << stats.norm_error << ", "
                  << stats.epsilon << "\n";
    } catch (const std::exception& e) {
        std::cerr << "llt_eval: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
