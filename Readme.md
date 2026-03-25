# ICCG
- CRS形式の疎行列に対する並列前処理付き共役勾配(CG)法
- 前処理なしのCG法から、複数の並列前処理を実装

## cg_crs
シンプルな共役勾配法

## dcg_crs
ヤコビ前処理（対角スケーリング）付きCG

## iccg_crs
ICCG

## picg_crs
並列前処理（MC法）付きICCG

## abmc_crs
ブロック並列前処理付きICCG

## 並列前進・後退代入
``` C++
Input:
  L      : lower triangular matrix (IC(0))
  r      : right-hand side vector
  C[1..m]: color sets (C_k ⊂ {1,…,n})

Output:
  y      : solution of Ly = r

for k = 1 to m do                    // color loop (sequential)
    parallel for each i ∈ C[k] do    // nodes in same color
        sum = 0
        for each j < i such that L[i,j] ≠ 0 do
            sum += L[i,j] * y[j]     // j ∈ C[1] ∪ … ∪ C[k−1]
        end for
        y[i] = (r[i] − sum) / L[i,i]
    end parallel for
end for
```