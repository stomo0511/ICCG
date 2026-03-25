#include "ic0.hpp"
#include "color.hpp"

// ic0.cpp：本体
IC0::IC0(const CRS& A, double shift) {
    // 下＋対角を抽出（値も持ってくる）
    CRS Alo = extract_lower_with_diag(A);
    n = Alo.n;
    build_ic0(Alo, shift);
    build_csc_from_crs();
}

CRS IC0::extract_lower_with_diag(const CRS& A) {
    CRS B; B.n = A.n;
    B.rowptr.assign(A.n+1, 0);
    std::vector<int> cols;
    std::vector<double> vals;
    cols.reserve(A.colind.size());
    vals.reserve(A.val.size());
    int nnz = 0;

    for (int i = 0; i < A.n; ++i) {
        B.rowptr[i] = nnz;
        for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
            int j = A.colind[k];
            if (j <= i) { // 下三角＋対角のみ
                cols.push_back(j);
                vals.push_back(A.val[k]);
                ++nnz;
            }
        }
    }
    B.rowptr[A.n] = nnz;
    B.colind.swap(cols);
    B.val.swap(vals);
    return B;
}

void IC0::build_ic0(const CRS& Alo, double shift) {
    // L のパターン＝Alo のパターン（同一）
    L_rowptr = Alo.rowptr;
    L_colind = Alo.colind;
    L_val    = Alo.val;   // いったん A の値をコピー。上書きで L にしていく。
    n = Alo.n;

    diag_pos.assign(n, -1);
    for (int i = 0; i < n; ++i) {
        for (int k = L_rowptr[i]; k < L_rowptr[i+1]; ++k) {
            if (L_colind[k] == i) { diag_pos[i] = k; break; }
        }
        if (diag_pos[i] < 0)
            throw std::runtime_error("IC(0): diagonal missing at row " + std::to_string(i));
    }

    // 行ごとに L を求める
    for (int i = 0; i < n; ++i) {
        // 1) 行iの各 j<i について L(i,j) を更新
        for (int kij = L_rowptr[i]; kij < L_rowptr[i+1]; ++kij) {
            int j = L_colind[kij];
            if (j >= i)
                break; // 行は昇順なのでここで止められる

            // L(i,0..j-1) と L(j,0..j-1) のドットを「二本ポインタ」で計算
            double s = 0.0;
            int pi = L_rowptr[i], pj = L_rowptr[j];
            const int ei = L_rowptr[i+1], ej = L_rowptr[j+1];

            while (pi < ei && pj < ej) {
                int ci = L_colind[pi];
                int cj = L_colind[pj];
                if (ci >= j || cj >= j)
                    break;      // どちらかが j に到達したら終了（<j のみ対象）
                if (ci == cj) {
                    s += L_val[pi] * L_val[pj]; ++pi; ++pj;
                }
                else if (ci < cj)
                    ++pi;
                else
                    ++pj;
            }

            const int dj = diag_pos[j];
            const double Ljj = L_val[dj];

            if (std::abs(Ljj) < 1e-30)
                throw std::runtime_error("IC(0): zero/near-zero pivot at row " + std::to_string(j));
            L_val[kij] = (L_val[kij] - s) / Ljj; // overwrite with L(i,j)
        }

        // 2) 対角 L(i,i)
        double s = 0.0;
        for (int k = L_rowptr[i]; k < L_rowptr[i+1]; ++k) {
            int j = L_colind[k];
            if (j >= i)
                break;
            double lij = L_val[k];
            s += lij * lij;
        }
        int di = diag_pos[i];
        double diag = L_val[di] - s;
        diag += shift;                 // 安全シフト
        if (diag <= 0.0)
            throw std::runtime_error("IC(0): non-positive pivot at row " + std::to_string(i));
        L_val[di] = std::sqrt(diag);
    }
}

void IC0::build_csc_from_crs() {
    // L(CRS) → L^T の CSC を構築（同じ値を列構造に並べる）
    LT_colptr.assign(n+1, 0);
    const int nnz = (int)L_colind.size();
    LT_rowind.resize(nnz);
    LT_val.resize(nnz);
    #ifdef PIC
    LT_diag_pos.assign(n, -1);
    #endif

    // 列ごとの数をカウント
    for (int k = 0; k < nnz; ++k) {
        int c = L_colind[k];
        ++LT_colptr[c+1];
    }
    // 累積和
    for (int c = 0; c < n; ++c)
        LT_colptr[c+1] += LT_colptr[c];

    // 一時カウンタ
    std::vector<int> ctr = LT_colptr;
    for (int i = 0; i < n; ++i) {
        for (int k = L_rowptr[i]; k < L_rowptr[i+1]; ++k) {
            int c = L_colind[k];
            int pos = ctr[c]++;
            LT_rowind[pos] = i;        // 列c に行i が来る（= L^T の (c,i)）
            LT_val[pos]    = L_val[k];
            #ifdef PIC
            if (i == c)
                LT_diag_pos[i] = pos;  // 対角位置
            #endif
        }
    }
}

void IC0::forward_solve(const std::vector<double>& b, std::vector<double>& y) const {
    // y.assign(n, 0.0);
    if ((int)y.size() != n) y.resize(n);    // y.assign(n, 0.0); をやめる
    for (int i = 0; i < n; ++i) {
        double s = b[i];
        double di = 1.0;
        for (int k = L_rowptr[i]; k < L_rowptr[i+1]; ++k) {
            int j = L_colind[k];
            if (j < i) s -= L_val[k] * y[j];
            else if (j == i) di = L_val[k];
            else break;
        }
        y[i] = s / di;
    }
}

void IC0::backward_solve(const std::vector<double>& y, std::vector<double>& x) const {
    // x.assign(n, 0.0);
    if ((int)x.size() != n) x.resize(n);    // x.assign(n, 0.0); をやめる
    for (int i = n-1; i >= 0; --i) {
        double s = y[i];
        double di = 1.0;
        // L^T の行 i は、L の列 i。CSC なら列 i がそのままアクセス可能。
        for (int k = LT_colptr[i]; k < LT_colptr[i+1]; ++k) {
            int j = LT_rowind[k];      // これは L^T の列 index に当たる
            if (j > i) s -= LT_val[k] * x[j];  // 上三角項
            else if (j == i) di = LT_val[k];   // 対角
            // j < i は関与しない
        }
        x[i] = s / di;
    }
}

void IC0::apply(const std::vector<double>& r, std::vector<double>& z) const {
    std::vector<double> y;
    forward_solve(r, y);          // L y = r
    backward_solve(y, z);         // L^T z = y
}

#ifdef PIC
// --- multi-color forward substitution: L y = r ---
void IC0::forward_solve_mc(int nc, const int* color_ptr, const int* rows_of_color,
                           const double* r, double* y) const
{
    #pragma omp parallel
    {
        for (int c = 0; c < nc; ++c) {
            #pragma omp for schedule(static)
            for (int k = color_ptr[c]; k < color_ptr[c+1]; ++k) {
                const int i   = rows_of_color[k];
                const int rs  = L_rowptr[i];
                const int re  = L_rowptr[i+1];
                const double di = L_val[diag_pos[i]];   // ★ O(1) で対角取得

                double acc = r[i];
                for (int jj = rs; jj < re; ++jj) {
                    const int j = L_colind[jj];
                    if      (j < i) acc -= L_val[jj] * y[j];
                    else if (j == i) /*skip*/;
                    else             break;           // 列は昇順
                }
                y[i] = acc / di;
            }
            // #pragma omp barrier;  // forの末尾に暗黙バリアがあるため不要
        }
    }
}

// --- multi-color backward substitution: L^T x = y ---
void IC0::backward_solve_mc(int nc, const int* color_ptr, const int* rows_of_color,
                            const double* y, double* x) const
{
    #pragma omp parallel
    {
        for (int c = nc - 1; c >= 0; --c) {
            #pragma omp for schedule(static)
            for (int k = color_ptr[c]; k < color_ptr[c+1]; ++k) {
                const int i   = rows_of_color[k];
                const int cs  = LT_colptr[i];
                const int ce  = LT_colptr[i+1];
                const double di = LT_val[LT_diag_pos[i]];  // ★ O(1) 対角

                double acc = y[i];
                for (int jj = cs; jj < ce; ++jj) {
                    const int j = LT_rowind[jj];     // L^T の列 index に相当
                    if (j > i)
                        acc -= LT_val[jj] * x[j];  // 上三角
                    else if (j == i) /*skip*/
                        ;
                    // j < i は依存しない
                }
                x[i] = acc / di;
            }
            // #pragma omp barrier;  // forの末尾に暗黙バリアがあるため不要
        }
    }
}

// void IC0::apply_mc(int nc, const int* color_ptr, const int* rows_of_color,
//                    const std::vector<double>& r, std::vector<double>& z) const
// {
//     if ((int)r.size() != n) throw std::runtime_error("IC0::apply_mc: r.size()!=n");
//     if ((int)z.size() != n) z.resize(n);

//     std::vector<double> y(n);
//     forward_solve_mc(nc, color_ptr, rows_of_color, r.data(), y.data());  // L y = r
//     backward_solve_mc(nc, color_ptr, rows_of_color, y.data(), z.data()); // L^T z = y
// }

// --- ABMC forward: 色×ブロックで並列、ブロック内は逐次 ---
void IC0::forward_solve_abmc(int nc,
                             const int* color_ptr_blk,
                             const int* blocks_of_color,
                             const int* block_ptr,
                             const int* rows_of_block,
                             const double* r, double* y) const
{
    #pragma omp parallel
    {
        for (int c = 0; c < nc; ++c) {
            #pragma omp for schedule(static)
            for (int kb = color_ptr_blk[c]; kb < color_ptr_blk[c+1]; ++kb) {
                const int b   = blocks_of_color[kb];
                const int rsb = block_ptr[b];
                const int reb = block_ptr[b+1];
                // ブロック内は new-index 昇順に逐次
                for (int t = rsb; t < reb; ++t) {
                    const int i   = rows_of_block[t];
                    const int rs  = L_rowptr[i];
                    const int re  = L_rowptr[i+1];
                    const double di = L_val[diag_pos[i]];

                    double acc = r[i];
                    for (int jj = rs; jj < re; ++jj) {
                        const int j = L_colind[jj];
                        if      (j < i)
                            acc -= L_val[jj] * y[j];
                        else if (j == i) /*skip*/
                            ;
                        else
                            break;
                    }
                    y[i] = acc / di;
                }
            }
        }
    }
}

// --- ABMC backward: 色降順×ブロックで並列、ブロック内は降順逐次 ---
void IC0::backward_solve_abmc(int nc,
                              const int* color_ptr_blk,
                              const int* blocks_of_color,
                              const int* block_ptr,
                              const int* rows_of_block,
                              const double* y, double* x) const
{
    #pragma omp parallel
    {
        for (int c = nc-1; c >= 0; --c) {
            #pragma omp for schedule(static)
            for (int kb = color_ptr_blk[c]; kb < color_ptr_blk[c+1]; ++kb) {
                const int b   = blocks_of_color[kb];
                const int rsb = block_ptr[b];
                const int reb = block_ptr[b+1];
                // ブロック内は new-index 降順に逐次
                for (int t = reb-1; t >= rsb; --t) {
                    const int i    = rows_of_block[t];
                    const int cs   = LT_colptr[i];
                    const int ce   = LT_colptr[i+1];
                    const double di = LT_val[LT_diag_pos[i]];

                    double acc = y[i];
                    for (int jj = cs; jj < ce; ++jj) {
                        const int j = LT_rowind[jj];
                        if      (j > i)
                            acc -= LT_val[jj] * x[j];
                        else if (j == i) /*skip*/
                            ;
                        // j<i は関与しない
                    }
                    x[i] = acc / di;
                }
            }
        }
    }
}
#endif // end of PIC

// 統一 apply: sched があれば MC, なければ通常
// void IC0::apply(const std::vector<double>& r, std::vector<double>& z, const ColorSchedule* sched) const
// {
//     // スケジュール未指定（nullptr）や不完全なら逐次版へフォールバック
//     if (!sched || sched->nc <= 0 || !sched->color_ptr || !sched->rows_of_color) {
//         apply(r, z);
//         return;
//     }
// #ifdef PIC
//     if ((int)z.size() != n) z.resize(n);
//     std::vector<double> y(n);
//     forward_solve_mc(sched->nc, sched->color_ptr, sched->rows_of_color, r.data(), y.data());
//     backward_solve_mc(sched->nc, sched->color_ptr, sched->rows_of_color, y.data(), z.data());

// #ifdef ABMC

// #endif // end of ABMC

// #else 
//     // PIC 非ビルドでも安全に動くよう逐次にフォールバック
//     apply(r, z);
// #endif
// }

void IC0::apply(const std::vector<double>& r, std::vector<double>& z, const ColorSchedule* sched) const
{
    if (!sched || sched->nc <= 0) { apply(r, z); return; }

#ifdef PIC
    if ((int)z.size() != n) z.resize(n);
    std::vector<double> y(n);

    if (sched->mode == ColorSchedule::ABMC_BLOCK &&
        sched->color_ptr_blk && sched->blocks_of_color &&
        sched->block_ptr && sched->rows_of_block)
    {
        forward_solve_abmc(sched->nc, sched->color_ptr_blk, sched->blocks_of_color,
                           sched->block_ptr, sched->rows_of_block, r.data(), y.data());
        backward_solve_abmc(sched->nc, sched->color_ptr_blk, sched->blocks_of_color,
                            sched->block_ptr, sched->rows_of_block, y.data(), z.data());
        return;
    }

    // 既存: 行粒度MC（PIC）
    if (sched->color_ptr && sched->rows_of_color) {
        forward_solve_mc(sched->nc, sched->color_ptr, sched->rows_of_color, r.data(), y.data());
        backward_solve_mc(sched->nc, sched->color_ptr, sched->rows_of_color, y.data(), z.data());
        return;
    }
#endif
    apply(r, z); // フォールバック
}
