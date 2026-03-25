#pragma once

#include <iostream>
#include <vector>
#include "color.hpp"

// 前処理なし
struct Identity {
    Identity() {}
    explicit Identity(const CRS&) {}

    void apply(const std::vector<double>& r, std::vector<double>& z, const ColorSchedule* /*sched*/) const {
        z = r;
    }

    std::string name() const { return "No preconditioner"; }
};

// Jacobi 前処理（対角スケーリング）: z = M^{-1} r with M = diag(A)
struct Jacobi {
    std::vector<double> inv_diag;  // 対角要素の逆数

    explicit Jacobi(const CRS& A) {
        inv_diag.assign(A.n, 1.0);
        for (int i = 0; i < A.n; ++i) {
            double d = 0.0;
            for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
                if (A.colind[k] == i) {
                    d = A.val[k];
                    break;
                }
            }
            // 念のための保護（ゼロ対角はCGの前提を満たさない）
            inv_diag[i] = (fabs(d) > 0.0) ? 1.0 / d : 1.0;
        }
    }

    void apply(const std::vector<double>& r, std::vector<double>& z, const ColorSchedule* /*sched*/) const {
        z.resize(r.size());
        for (size_t i = 0; i < r.size(); ++i)
            z[i] = inv_diag[i] * r[i];
    }

    std::string name() const { return "Jacobi preconditioner"; }
};