#include "llt_err.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

long double FrobeniusLowerSquared(const CRS& A_lower) {
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

double RowDotL(const std::vector<int>& rowptr,
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

long double InnerALLtLower(const CRS& A_lower, const IC0& M) {
    const auto& rowptr = M.l_rowptr();
    const auto& colind = M.l_colind();
    const auto& val = M.l_val();

    long double s = 0.0L;
    for (int i = 0; i < A_lower.n; ++i) {
        for (int k = A_lower.rowptr[i]; k < A_lower.rowptr[i + 1]; ++k) {
            const int j = A_lower.colind[k];
            const long double llt_ij = RowDotL(rowptr, colind, val, i, j);
            const long double aij = A_lower.val[k];
            s += (j == i) ? aij * llt_ij : 2.0L * aij * llt_ij;
        }
    }
    return s;
}

void BuildColumnsOfL(const IC0& M,
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

long double FrobeniusLLtSquared(const IC0& M) {
    const int n = M.size();
    const auto& L_rowptr = M.l_rowptr();
    const auto& L_colind = M.l_colind();
    const auto& L_val = M.l_val();

    std::vector<int> colptr, rowind;
    std::vector<double> cval;
    BuildColumnsOfL(M, colptr, rowind, cval);

    std::vector<int> marker(n, -1);
    std::vector<double> acc(n, 0.0);
    std::vector<int> touched;
    long double norm2 = 0.0L;

    // ||LL^T||_F^2 = ||L^T L||_F^2.  Each column of G=L^T L is formed in
    // a temporary sparse accumulator, so fill outside A is counted without
    // storing the whole product.
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

std::vector<std::vector<double>> DenseFromLower(const CRS& A_lower) {
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

} // namespace

LLTErrorResult EvaluateLLTError(const CRS& A_lower, const IC0& M) {
    if (A_lower.n != M.size())
        throw std::runtime_error("EvaluateLLTError: matrix size mismatch");

    LLTErrorResult result;
    const long double norm_A2 = FrobeniusLowerSquared(A_lower);
    const long double norm_LLt2 = FrobeniusLLtSquared(M);
    const long double inner = InnerALLtLower(A_lower, M);
    long double err2 = norm_A2 + norm_LLt2 - 2.0L * inner;

    if (err2 < 0.0L && err2 > -100.0L * std::numeric_limits<double>::epsilon() * norm_A2)
        err2 = 0.0L;
    err2 = std::max(err2, 0.0L);

    result.norm_A_F = std::sqrt((double)norm_A2);
    result.norm_LLt_F = std::sqrt((double)norm_LLt2);
    result.inner_A_LLt = (double)inner;
    result.norm_error_F = std::sqrt((double)err2);
    result.epsilon_F = (norm_A2 > 0.0L)
        ? result.norm_error_F / result.norm_A_F
        : std::numeric_limits<double>::quiet_NaN();
    return result;
}

double DenseCheckLLTEpsilon(const CRS& A_lower, const IC0& M) {
    if (A_lower.n != M.size())
        throw std::runtime_error("DenseCheckLLTEpsilon: matrix size mismatch");

    const int n = A_lower.n;
    auto A = DenseFromLower(A_lower);
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

    const auto& rowptr = M.l_rowptr();
    const auto& colind = M.l_colind();
    const auto& val = M.l_val();
    for (int i = 0; i < n; ++i) {
        for (int k = rowptr[i]; k < rowptr[i + 1]; ++k)
            L[i][colind[k]] = val[k];
    }

    long double norm_A2 = 0.0L;
    long double err2 = 0.0L;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            long double llt = 0.0L;
            const int lim = std::min(i, j);
            for (int k = 0; k <= lim; ++k)
                llt += (long double)L[i][k] * L[j][k];
            const long double a = A[i][j];
            const long double d = a - llt;
            norm_A2 += a * a;
            err2 += d * d;
        }
    }
    return (norm_A2 > 0.0L)
        ? std::sqrt((double)(err2 / norm_A2))
        : std::numeric_limits<double>::quiet_NaN();
}

CRS LowerCRSWithDiag(const CRS& A) {
    CRS B;
    B.n = A.n;
    B.rowptr.assign(A.n + 1, 0);
    B.colind.reserve(A.colind.size());
    B.val.reserve(A.val.size());

    int nnz = 0;
    for (int i = 0; i < A.n; ++i) {
        B.rowptr[i] = nnz;
        for (int k = A.rowptr[i]; k < A.rowptr[i + 1]; ++k) {
            const int j = A.colind[k];
            if (j <= i) {
                B.colind.push_back(j);
                B.val.push_back(A.val[k]);
                ++nnz;
            }
        }
    }
    B.rowptr[A.n] = nnz;
    return B;
}

long long SymmetricNNZFromLower(const CRS& A_lower) {
    long long nnz = 0;
    for (int i = 0; i < A_lower.n; ++i) {
        for (int k = A_lower.rowptr[i]; k < A_lower.rowptr[i + 1]; ++k)
            nnz += (A_lower.colind[k] == i) ? 1 : 2;
    }
    return nnz;
}

std::string MatrixNameFromPath(const std::string& path) {
    size_t slash = path.find_last_of("/\\");
    std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
    size_t dot = name.find_last_of('.');
    if (dot != std::string::npos) name.erase(dot);
    return name;
}
