//
// CRS形式の疎行列に対する並列前処理付き共役勾配(CG)法
//
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <limits>
#include "crs_io.hpp"

#ifndef IC
#include "precond.hpp"
#else
#include "ic0.hpp"
#endif

#ifdef PIC
#include "color.hpp"
#endif

#ifdef ABMC
#include "block.hpp"
#endif

#ifdef USEMKL
#include <mkl.h>
#endif

using namespace std::chrono;

// 内積
static inline double dot(const std::vector<double>& a, const std::vector<double>& b) {
    long double s = 0.0L;
    const size_t n = a.size();
    #pragma omp parallel for reduction(+:s) schedule(static)
    for (size_t i = 0; i < n; ++i)
        s += (long double)a[i] * b[i];
    return (double)s;
}

// 2ノルム
static inline double nrm2(const std::vector<double>& a) {
    long double s = 0.0L;
    const size_t n = a.size();
    #pragma omp parallel for reduction(+:s) schedule(static)
    for (size_t i = 0; i < n; ++i)
        s += (long double)a[i] * a[i];
    return (double)std::sqrt((double)s);
}

// 疎行列ベクトル積
// y = A*x
static void spmv(const CRS& A, const std::vector<double>& x, std::vector<double>& y) {
    int n = A.n;
    // y.assign(n, 0.0);
    if ((int)y.size() != n) y.resize(n);  // 再割当て/ゼロ埋めを避ける

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
        for (int k = rs; k < re; ++k) {
            sum += aval[k] * xp[colind[k]];
        }
        yp[i] = sum;
    }
    #endif
}

// CG法の結果
struct CGResult {
    int iters = 0;
    double rel_resid = NAN;
    bool converged = false;
};

// CG法
template<class Precond>
CGResult conjugate_gradient(
    const CRS& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    const Precond& M,           // 前処理行列
    int max_iter,
    double tol,
    const ColorSchedule* sched
) {
    const int n = A.n;
    std::vector<double> r(n), p(n), z(n), Ap(n), Ax(n);

    spmv(A, x, Ax);
    // r0 = b - A x0
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
        r[i] = b[i] - Ax[i];

    double normb = nrm2(b);
    if (normb == 0.0)
        normb = 1.0;

    M.apply(r, z, sched);   // z0 = M^{-1} r0

    p = z;
    double rz_old = dot(r, z);

    CGResult res; 
    res.rel_resid = nrm2(r)/normb;
    if (res.rel_resid < tol) { 
        res.converged = true; 
        res.iters = 0; 
        return res; 
    }

    for (int k = 0; k < max_iter; ++k) {
        spmv(A, p, Ap);
        double pAp = dot(p, Ap);
        double eps = std::numeric_limits<double>::epsilon();
        double thr = eps * nrm2(p) * nrm2(Ap);
        if (pAp <= thr) {
            res.iters = k;
            return res; // 方向喪失（丸め）を検知
        }

        double alpha = rz_old / pAp;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rel = nrm2(r) / normb;
        res.iters = k+1;
        res.rel_resid = rel;
        if (rel < tol) {
            res.converged = true;
            return res;
        }

        M.apply(r, z, sched);   // z = M^{-1} r
        double rz_new = dot(r, z);
        double beta   = rz_new / rz_old;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i)
            p[i] = z[i] + beta * p[i];

        rz_old = rz_new;
    }
    return res;
} // End of conjugate_gradient

////////////////////////////////////////////////////////////
// デバッグ用関数（デバッグ後は削除可）
// ===== 2.1: 同色ブロック間エッジ検査 =====
// A: もとの行列（full/upper/lower どれでもOK。対称なら i<j のみ数える）
// block_of_old: 旧ノード -> ブロック (0-based)
// block_color : ブロック -> 色 (0-based)
// nb: ブロック数, nc: 色数
void DebugCheck_InterBlockSameColorEdges(
    const CRS& A,
    const std::vector<int>& block_of_old,
    const std::vector<int>& block_color,
    int nb, int nc)
{
    const int n = A.n;
    auto in_range = [&](int x, int lo, int hi){ return (x>=lo && x<hi); };
    if ((int)block_of_old.size()!=n) throw std::runtime_error("block_of_old size mismatch");
    if ((int)block_color.size()!=nb) throw std::runtime_error("block_color size mismatch");

    std::vector<long long> conflict_by_color(nc, 0);
    long long total_conflicts = 0;

    long long undirected_edges_seen = 0; // optional: sanity
    for (int i=0;i<n;++i){
        for (int k=A.rowptr[i]; k<A.rowptr[i+1]; ++k){
            int j = A.colind[k];
            if (j<=i) continue; // i<j の片側だけ数える（対称想定）

            int bi = block_of_old[i];
            int bj = block_of_old[j];
            if (!in_range(bi,0,nb) || !in_range(bj,0,nb))
                throw std::runtime_error("block id out of range while scanning edges");

            if (bi==bj) continue; // 同一ブロック内は対象外（ブロック間だけ見る）

            int ci = block_color[bi];
            int cj = block_color[bj];
            if (!in_range(ci,0,nc) || !in_range(cj,0,nc))
                throw std::runtime_error("color id out of range while scanning edges");

            ++undirected_edges_seen;

            if (ci==cj){
                ++conflict_by_color[ci];
                ++total_conflicts;
            }
        }
    }

    std::cerr << "[ABMC check 2.1] same-color inter-block edges per color:\n";
    for (int c=0;c<nc;++c){
        std::cerr << "  color " << c << " : conflicts=" << conflict_by_color[c] << "\n";
    }
    std::cerr << "  TOTAL conflicts=" << total_conflicts
              << " (undirected edges seen=" << undirected_edges_seen << ")\n";

    if (total_conflicts==0){
        std::cerr << "  ==> OK: No edges between blocks that share the same color.\n";
    }else{
        std::cerr << "  ==> NG: There exist edges between blocks of the same color.\n";
    }
}

// 下＋対角を抽出（IC0::extract_lower_with_diag と同等の処理）
static CRS extract_lower_with_diag_like(const CRS& A_full_or_lower) {
    CRS B; B.n = A_full_or_lower.n;
    B.rowptr.assign(B.n+1, 0);
    std::vector<int> cols;
    std::vector<double> vals;
    cols.reserve(A_full_or_lower.colind.size());
    vals.reserve(A_full_or_lower.val.size());
    int nnz = 0;

    for (int i = 0; i < B.n; ++i) {
        B.rowptr[i] = nnz;
        for (int k = A_full_or_lower.rowptr[i]; k < A_full_or_lower.rowptr[i+1]; ++k) {
            int j = A_full_or_lower.colind[k];
            if (j <= i) {
                cols.push_back(j);
                vals.push_back(A_full_or_lower.val[k]);
                ++nnz;
            }
        }
    }
    B.rowptr[B.n] = nnz;
    B.colind.swap(cols);
    B.val.swap(vals);
    return B;
}

// ===== 3: スケジュール危険依存（同色 j<i 依存）の検査 =====
// Ap_full : 置換後の行列（上下対称に展開済みを推奨）
// new_of_old : 旧->新 置換
// color_old  : 旧ノード -> 色 (0-based)
// （ABMCの診断を強める場合）block_of_old, block_color を渡すと、
//   危険依存が「同ブロック内」か「異ブロック間」かを分類する
void DebugCheck_ScheduleHazards_OnLower(
    const CRS& Ap_full,
    const std::vector<int>& new_of_old,
    const std::vector<int>& color_old,
    const std::vector<int>* block_of_old_opt,   // nullptr 可
    const std::vector<int>* block_color_opt)    // 未使用でも可
{
    const int n = Ap_full.n;
    if ((int)new_of_old.size()!=n || (int)color_old.size()!=n)
        throw std::runtime_error("size mismatch in DebugCheck_ScheduleHazards_OnLower");

    // new->old
    std::vector<int> old_of_new(n, -1);
    for (int old_id=0; old_id<n; ++old_id){
        int nid = new_of_old[old_id];
        if (nid<0 || nid>=n) throw std::runtime_error("new_of_old out of range");
        old_of_new[nid] = old_id;
    }

    // new 空間での色/ブロック（0-based）
    std::vector<int> color_new(n, -1);
    std::vector<int> block_new(n, -1);
    const bool has_block = (block_of_old_opt && (int)block_of_old_opt->size()==n);
    for (int i=0;i<n;++i){
        int old_i = old_of_new[i];
        color_new[i] = color_old[old_i];
        if (has_block) block_new[i] = (*block_of_old_opt)[old_i];
    }

    // L のパターン（Ap の下＋対角）を作る
    CRS Lpat = extract_lower_with_diag_like(Ap_full); // IC(0)の L と同パターンでOK

    // 色数を推定（最大+1）
    int nc = 0;
    for (int c : color_new) nc = std::max(nc, c+1);

    std::vector<long long> hazard_by_color(nc, 0);
    std::vector<long long> hazard_same_block_by_color(nc, 0);
    std::vector<long long> hazard_cross_block_by_color(nc, 0);

    long long total_hazard = 0;

    for (int i=0;i<n;++i){
        int ci = color_new[i];
        for (int kk=Lpat.rowptr[i]; kk<Lpat.rowptr[i+1]; ++kk){
            int j = Lpat.colind[kk];
            if (j>=i) break;      // j<i のみ（前進代入の依存）
            int cj = color_new[j];
            if (cj==ci){
                ++hazard_by_color[ci];
                ++total_hazard;
                if (has_block){
                    if (block_new[i]==block_new[j]) ++hazard_same_block_by_color[ci];
                    else                             ++hazard_cross_block_by_color[ci]; // これは本来あり得ない（ABMC前提）
                }
            }
        }
    }

    std::cerr << "[ABMC check 3] same-color lower deps (j<i) by color:\n";
    for (int c=0;c<nc;++c){
        std::cerr << "  color " << c
                  << " : hazards=" << hazard_by_color[c];
        if (has_block){
            std::cerr << " (same-block=" << hazard_same_block_by_color[c]
                      << ", cross-block=" << hazard_cross_block_by_color[c] << ")";
        }
        std::cerr << "\n";
    }
    std::cerr << "  TOTAL hazards=" << total_hazard << "\n";

    if (total_hazard==0){
        std::cerr << "  ==> OK: MC行粒度スケジュールでも依存ハザードなし（PIC互換の色）\n";
    }else{
        std::cerr << "  ==> 注意: 行粒度MCでは同色内に j<i 依存あり。ABMCでは"
                     "『色×ブロック並列・ブロック内逐次』が必要。\n";
    }
}
////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
    // cg_crs, dcg_crs
    #ifndef PIC
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " matrix.mtx [tol=1e-14] [max_iter=10000]\n";
        return 1;
    }
    std::string path = argv[1];  // matrix data file name
    double tol = (argc >= 3) ? std::atof(argv[2]) : 1e-14;    // 収束判定
    int max_iter = (argc >= 4) ? std::atoi(argv[3]) : 10000;  // 最大反復回数
    #endif  // end of PIC

    // iccg_crs, piccg_crs
    #ifdef PIC
    #ifndef ABMC
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " matrix.mtx color.col [tol=1e-14] [max_iter=10000]\n";
        return 1;
    }
    std::string path = argv[1];  // matrix data file name
    std::string cpath = argv[2]; // color define file name
    double tol = (argc >= 4) ? std::atof(argv[3]) : 1e-14;    // 収束判定
    int max_iter = (argc >= 5) ? std::atoi(argv[4]) : 10000;  // 最大反復回数

    // abmc_crs
    #else
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " matrix.mtx block.blk bcolor.bcol [tol=1e-14] [max_iter=10000]\n";
        return 1;
    }
    std::string path = argv[1];  // matrix data file name
    std::string bpath = argv[2]; // block define file name
    std::string cpath = argv[3]; // block-color define file name
    double tol = (argc >= 5) ? std::atof(argv[4]) : 1e-14;    // 収束判定
    int max_iter = (argc >= 6) ? std::atoi(argv[5]) : 10000;  // 最大反復回数
    #endif // end of ABMC
    #endif // end of PIC

    // 行列データの読み込み (MatrixMarket形式ファイル -> CRS形式)
    CRS A = read_matrix_market_crs(path);
    #ifdef USEMKL
    A.build_mkl_handle();
    #endif

    if (A.n <= 0) {
        std::cerr << "Invalid matrix size.\n"; return 1;
    }

    // 色ファイルの読み込み
    #ifdef PIC
    int nc = 0;  // # of colors

    // iccg_crs, piccg_crs
    #ifndef ABMC
    // ノード彩色情報の読み込み
    std::vector<int> color = ReadColorFile_1Based(cpath, A.n, nc);

    // abmc_crs
    #else
    int nbs = 0;   // # of blocks
    int nbc = 0;  // # of block colors

    // ブロック情報の読み込み
    std::vector<int> block_of_old = ReadBlockFile_1Based(bpath, A.n, nbs);
    // ブロック彩色情報の読み込み
    std::vector<int> block_color  = ReadBlockColorFile_1Based(cpath, nbs, nbc);

    // ノード彩色情報の構築
    std::vector<int> color(A.n);
    for (int i = 0; i < A.n; ++i)
        color[i] = block_color[ block_of_old[i] ];

    nc = nbc;
    #endif // end of ABMC

    // 置換と PtAP を適用
    #ifndef ABMC
    std::vector<int> new_of_old = BuildPermutationByColor(color);
    #else
    std::vector<int> new_of_old = BuildPermutationByBlockColor(block_of_old, block_color);
    #endif

    CRS Ap = Permute_PtAP_LowerCRS_to_LowerCRS(A, new_of_old);
    Ap = expand_lower_to_full(Ap);  // 上三角部分に要素を充填

    // ABMC スケジュールを new 空間で構築
    #ifdef ABMC
    std::vector<int> color_ptr_blk, blocks_of_color, block_ptr, rows_of_block;
    BuildABMCSchedule(new_of_old, block_of_old, block_color, nbs, nbc,
                      color_ptr_blk, blocks_of_color, block_ptr, rows_of_block);
    #endif // end of ABMC

    #ifdef USEMKL
    Ap.build_mkl_handle();
    #endif

    std::vector<int> color_ptr;
    std::vector<int> rows_of_color;
    BuildColorPtrRows_FromPermutationAndColor(color, nc, new_of_old, color_ptr, rows_of_color);
    #endif // end of PIC

    std::vector<double> x(A.n, 0.0);  // 解ベクトル（初期値はすべて 0）
    std::vector<double> b(A.n, 1.0);  // RHSベクトル（すべて 1）

    const CRS& A_run =
    #ifdef PIC
    Ap;
    #else
    A;
    #endif // end of PIC

    #if defined(NOPRE)
    Identity M;           // 前処理なし
    #elif defined(JAC)
    Jacobi M(A_run);      // 対角スケーリング
    #elif defined(IC)
    double shift = 0.0;
    IC0 M(A_run, shift);  // IC(0) 前処理
    #endif

    // std::cout << M.name() << " applied.\n";
    // #ifdef IC
    // std::cout << "IC shift=" << shift << "\n";
    // #endif

    const ColorSchedule* sched_ptr = nullptr;
    ColorSchedule sched;
    
    #ifdef PIC
    #ifndef ABMC
    sched.nc = nc;
    sched.color_ptr = color_ptr.data();
    sched.rows_of_color = rows_of_color.data();
    sched_ptr = &sched;
    // std::cout << "Number of colors (MC): " << nc << "\n";
    #else
    sched.mode = ColorSchedule::ABMC_BLOCK;
    sched.nc = nbc;
    sched.nb = nbs;
    sched.color_ptr_blk  = color_ptr_blk.data();
    sched.blocks_of_color= blocks_of_color.data();
    sched.block_ptr      = block_ptr.data();
    sched.rows_of_block  = rows_of_block.data();
    sched_ptr = &sched;
    // std::cout << "Number of colors (ABMC): " << nbc << "\n";
    #endif
    #endif // end of PIC

    // #ifdef ABMC
    /////////////////////////////////////////////
    // デバッグ用チェック
    // 2.1: 同色ブロック間エッジ検査（元の A でOK）
    // DebugCheck_InterBlockSameColorEdges(
    //     A,                // 読み込み直後の元行列
    //     block_of_old,     // 旧ノード -> ブロック
    //     block_color,      // ブロック -> 色
    //     nbs, nbc);

    // 3: スケジュール危険依存検査（Ap の下＋対角で j<i 依存を見る）
    // DebugCheck_ScheduleHazards_OnLower(
    //     Ap,               // expand_lower_to_full 済みの Ap（上下あり）
    //     new_of_old,
    //     color,            // 旧ノード -> 色（ABMCでは block_colorをノードに展開したもの）
    //     &block_of_old,
    //     &block_color);
    // #endif
    /////////////////////////////////////////////

    // CG反復部分
    auto start = steady_clock::now();  // 計測開始
    auto out = conjugate_gradient(A_run, b, x, M, max_iter, tol, sched_ptr);
    auto end = steady_clock::now();    // 計測終了
    auto elapsed = duration_cast<milliseconds>(end - start).count();

    std::cout << "cg_time [ms]=" << elapsed
              << ", converged=" << out.converged
              << ", iters=" << out.iters
              << ", rel_resid=" << out.rel_resid;

    // 残差 ||Ax - b||/||b|| のチェック
    std::vector<double> Ax;
    spmv(A_run, x, Ax);

    double nr=0.0, nb=0.0;
    #pragma omp parallel for reduction(+:nr,nb) schedule(static)
    for (int i = 0; i < (int)Ax.size(); ++i) {
        double d = Ax[i]-b[i];
        nr += d*d;
        nb += b[i]*b[i];
    }
    std::cout << ", ||Ax-b||/||b|| = " << std::sqrt(nr)/std::sqrt(nb) << "\n";

    return 0;
}
