#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "crs_io.hpp"

// ---- 内部ユーティリティ ----
static std::string to_lower(std::string s) {
    for (auto& ch : s) ch = std::tolower(static_cast<unsigned char>(ch));
    return s;
}

static CRS coo_to_crs(int nrows, int ncols,
                      std::vector<int>& I, std::vector<int>& J, std::vector<double>& V)
{
    (void)ncols; // 未使用だが将来チェック用途に残す
    const size_t nnz = I.size();
    std::vector<size_t> ord(nnz);
    std::iota(ord.begin(), ord.end(), 0);

    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
        if (I[a] != I[b]) return I[a] < I[b];
        return J[a] < J[b];
    });

    std::vector<int> rows_u; rows_u.reserve(nnz);
    std::vector<int> cols_u; cols_u.reserve(nnz);
    std::vector<double> vals_u; vals_u.reserve(nnz);

    for (size_t idx = 0; idx < nnz; ) {
        int r = I[ord[idx]];
        int c = J[ord[idx]];
        double s = 0.0;

        do {
            s += V[ord[idx]];
            ++idx;
        } while (idx < nnz && I[ord[idx]] == r && J[ord[idx]] == c);

        rows_u.push_back(r);
        cols_u.push_back(c);
        vals_u.push_back(s);
    }

    CRS A;
    A.n = nrows;
    A.rowptr.assign(nrows+1, 0);
    A.colind.resize(cols_u.size());
    A.val.resize(vals_u.size());

    for (size_t k = 0; k < rows_u.size(); ++k)
        ++A.rowptr[rows_u[k]+1];
    for (int i = 0; i < nrows; ++i)
        A.rowptr[i+1] += A.rowptr[i];

    std::vector<int> ctr = A.rowptr;
    for (size_t k = 0; k < rows_u.size(); ++k) {
        int r = rows_u[k];
        int pos = ctr[r]++;
        A.colind[pos] = cols_u[k];
        A.val[pos] = vals_u[k];
    }
    return A;
}

// ---- 公開関数 ----
CRS read_matrix_market_crs(const std::string& filepath) {
    std::ifstream fin(filepath);
    if (!fin)
        throw std::runtime_error("Failed to open: " + filepath);

    std::string banner, mtx, storage, field, symmetry;
    {
        std::string line;
        if (!std::getline(fin, line))
            throw std::runtime_error("Empty file");
        std::istringstream iss(line);
        if (!(iss >> banner >> mtx >> storage >> field >> symmetry))
            throw std::runtime_error("Invalid MatrixMarket header");
        if (banner != "%%MatrixMarket")
            throw std::runtime_error("Not a MatrixMarket file");
        mtx = to_lower(mtx);
        storage = to_lower(storage);
        field = to_lower(field);
        symmetry = to_lower(symmetry);
        if (mtx != "matrix" || storage != "coordinate")
            throw std::runtime_error("Only 'matrix coordinate' is supported");
        if (!(field == "real" || field == "integer" || field == "pattern"))
            throw std::runtime_error("Unsupported field: " + field);
        if (!(symmetry == "general" || symmetry == "symmetric"))
            throw std::runtime_error("Unsupported symmetry: " + symmetry);
    }

    auto skip_comments = [&](std::istream& in) {
        std::streampos pos;
        std::string line;
        while (true) {
            pos = in.tellg();
            if (!std::getline(in, line)) break;
            if (!line.empty() && line[0] == '%') continue;
            in.seekg(pos);
            break;
        }
    };
    skip_comments(fin);

    int nrows=0, ncols=0, nnz=0;
    {
        std::string line;
        if (!std::getline(fin, line))
            throw std::runtime_error("Missing size line");
        std::istringstream iss(line);
        if (!(iss >> nrows >> ncols >> nnz))
            throw std::runtime_error("Invalid size line");
    }

    std::vector<int> I; I.reserve(nnz * (symmetry=="symmetric" ? 2 : 1));
    std::vector<int> J; J.reserve(I.capacity());
    std::vector<double> V; V.reserve(I.capacity());

    auto push_entry = [&](int i, int j, double v) {
        if (i < 0 || i >= nrows || j < 0 || j >= ncols)
            throw std::runtime_error("Index out of range");
        I.push_back(i); J.push_back(j); V.push_back(v);
    };

    for (int k = 0; k < nnz; ++k) {
        std::string line;
        if (!std::getline(fin, line))
            throw std::runtime_error("Unexpected EOF");
        if (line.empty() || line[0] == '%') {
            --k;
            continue;
        }
        std::istringstream iss(line);
        long long ii, jj;
        if (!(iss >> ii >> jj))
            throw std::runtime_error("Invalid entry line (indices)");
        double vv = 1.0;
        if (field == "real" || field == "integer") {
            if (!(iss >> vv))
                throw std::runtime_error("Invalid entry line (value)");
        } else {
            vv = 1.0;
        }
        int i = static_cast<int>(ii - 1);
        int j = static_cast<int>(jj - 1);
        push_entry(i, j, vv);
        if (symmetry == "symmetric" && i != j) push_entry(j, i, vv);
    }

    return coo_to_crs(nrows, ncols, I, J, V);
}

// 下三角CRS (j<=i) → 上下両三角CRSへ展開（対称）
// 対角は1回、非対角は(i,j)と(j,i)の2つを出力
CRS expand_lower_to_full(const CRS& L) {
    CRS F; F.n = L.n;
    // まず出力nnzを数える
    int n = L.n;
    std::vector<int> nnz_row(n, 0);
    int nnzF = 0;

    for (int i = 0; i < n; ++i) {
        for (int k = L.rowptr[i]; k < L.rowptr[i+1]; ++k) {
            int j = L.colind[k];
            if (j == i) {                 // 対角は一つ
                nnz_row[i] += 1; nnzF += 1;
            } else {                       // 非対角は2倍（(i,j),(j,i)）
                nnz_row[i] += 1;
                nnz_row[j] += 1;
                nnzF += 2;
            }
        }
    }

    F.rowptr.assign(n+1, 0);
    for (int i = 0; i < n; ++i) F.rowptr[i+1] = F.rowptr[i] + nnz_row[i];
    F.colind.resize(nnzF);
    F.val.resize(nnzF);

    std::vector<int> cursor = F.rowptr;
    for (int i = 0; i < n; ++i) {
        for (int k = L.rowptr[i]; k < L.rowptr[i+1]; ++k) {
            int j = L.colind[k];
            double v = L.val[k];
            // (i,j)
            int pos = cursor[i]++;
            F.colind[pos] = j;
            F.val[pos]    = v;
            if (j != i) {
                // (j,i)
                int pos2 = cursor[j]++;
                F.colind[pos2] = i;
                F.val[pos2]    = v;
            }
        }
    }

    // 行内の列を昇順に整える（任意だが、以後の処理が楽になる）
    for (int i = 0; i < n; ++i) {
        int rs = F.rowptr[i], re = F.rowptr[i+1];
        auto first = F.colind.begin() + rs;
        auto last  = F.colind.begin() + re;
        // valも同時に並び替えたいのでインデックスソート
        std::vector<int> idx(re - rs);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b){
            return F.colind[rs + a] < F.colind[rs + b];
        });
        std::vector<int>   newC(idx.size());
        std::vector<double> newV(idx.size());
        for (size_t t=0; t<idx.size(); ++t) {
            newC[t] = F.colind[rs + idx[t]];
            newV[t] = F.val[rs + idx[t]];
        }
        std::copy(newC.begin(), newC.end(), F.colind.begin()+rs);
        std::copy(newV.begin(), newV.end(), F.val.begin()+rs);
    }
    return F;
}