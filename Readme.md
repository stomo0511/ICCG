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

## 並列前進代入
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

## 並列ブロック前進代入
``` C++
#pragma omp paralel {
#pragma omp single {
  for p = 0 to nc-1 do
    for each block a in color p do
      // rhs^p_a = r^p_a を初期化する task
      #pragma omp task depend(out: rhs[p][a])
      {
        rhs[p][a] = r[p][a]
      }

      // 過去色 q < p からの寄与を rhs^p_a から引く
      for q = 0 to p-1 do
        for each block b in color q do
          if L^{p,q}_{a,b} is nonzero then

            #pragma omp task depend(in: y[q][b]) depend(inout: rhs[p][a])
            {
              temp = SpMV(L^{p,q}_{a,b}, y[q][b])
              rhs[p][a] = rhs[p][a] - temp
            }
          end if
        end for
      end for

      // 対角ブロックの前進代入
      #pragma omp task depend(inout: rhs[p][a]) depend(out: y[p][a])
      {
        FS(L^{p,p}_{a,a}, rhs[p][a], y[p][a])
      }
    end for
  end for
}  
}
```