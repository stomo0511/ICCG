#!/bin/bash

export OMP_NUM_THREADS=32

echo "time, converged, iters, rel_resid, resid"

# 実行回数
N=10

CMD="./abmc_crs"
TOL="1.0E-08"
MIT="10000"
DAT="./0_florida/0525825_parabolic_fem.mtx"

Bs=(128)

for B in "${Bs[@]}"
do
    BLK="./0525825_parabolic_fem_abmc_B${B}_p1.blk"
    BCOL="./0525825_parabolic_fem_abmc_B${B}_p1.bcol"

    for i in $(seq 1 $N)
    do
        $CMD $DAT $BLK $BCOL $TOL $MIT
    done
done
