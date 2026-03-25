#pragma once
#include <vector>
#include <string>

#ifdef USEMKL
#include <mkl.h>
#endif

// CRS形式の疎行列を表す構造体
struct CRS {
    int n;                        // matrix size (n x n)
    std::vector<int> rowptr;      // size n+1
    std::vector<int> colind;      // size nnz
    std::vector<double> val;      // size nnz

    #ifdef USEMKL
    // oneMKL handle
    sparse_matrix_t handle = nullptr;
    matrix_descr   descr;

    // 0-based CRS から MKL ハンドル作成
    void build_mkl_handle() {
        if (handle) mkl_sparse_destroy(handle);
        // 0-based で OK（row_start=ia[0..n-1], row_end=ia[1..n]）
        auto status = mkl_sparse_d_create_csr(
            &handle,
            SPARSE_INDEX_BASE_ZERO,
            (MKL_INT)n,
            (MKL_INT)n,
            (MKL_INT*)rowptr.data(),            // row_start[i] = ia[i]
            (MKL_INT*)(rowptr.data() + 1),      // row_end[i]   = ia[i+1]
            (MKL_INT*)colind.data(),
            (double*)val.data()
        );
        if (status != SPARSE_STATUS_SUCCESS) {
            throw std::runtime_error("mkl_sparse_d_create_crs failed");
        }

        // y = A x 用の記述子（GENERAL で十分）
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        descr.mode = SPARSE_FILL_MODE_FULL;
        descr.diag = SPARSE_DIAG_NON_UNIT;

        // ヒント & 最適化
        mkl_sparse_set_mv_hint(handle, SPARSE_OPERATION_NON_TRANSPOSE, descr, /*expected_calls=*/1000);
        mkl_sparse_optimize(handle);
    }
    #endif
};

// Matrix Market (*.mtx) を読み込んで CRS を返す（general/symmetric, real/integer/pattern）
CRS read_matrix_market_crs(const std::string& filepath);

// 下三角CRS (j<=i) → 上下両三角CRSへ展開（対称）
// 対角は1回、非対角は(i,j)と(j,i)の2つを出力
CRS expand_lower_to_full(const CRS& L);
