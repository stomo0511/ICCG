#!/usr/bin/sh

export OMP_NUM_THREADS=32

echo "time, converged, iters, rel_resid, resid"

# 実行回数
N=10

CMD="./piccg_crs"
TOL="1.0E-08"
MIT="10000"
DAT="./0_florida/525825_parabolic_fem.mtx"
COL="./20250915/525825_parabolic_fem_greedy.col"

for i in $(seq 1 $N)
do
    $CMD $DAT $COL $TOL $MIT
done