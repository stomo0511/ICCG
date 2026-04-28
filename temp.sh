#!/usr/bin/bash

export OMP_NUM_THREADS=32

echo "mat, cnv, time, iters, rel_resid, resid"

# 実行回数
N=10

CMD="./abmc_crs"
TOL="1.0E-08"
MIT="10000"

for i in $(seq 1 $N); do
    echo -n "Emilia_923, "
    $CMD ../0_florida/Emilia_923.mtx ../ABMC/Emilia_923_rblock_s14000.blk ../ABMC/Emilia_923_rblock_s14000.bcol $TOL $MIT
done

for i in $(seq 1 $N); do
    echo -n "Hook_1498, "
    $CMD ../0_florida/Hook_1498.mtx ../ABMC/Hook_1498_rblock_s3000.blk ../ABMC/Hook_1498_rblock_s3000.bcol $TOL $MIT
done

for i in $(seq 1 $N); do
    echo -n "thermal2, "
    $CMD ../0_florida/thermal2.mtx ../ABMC/thermal2_rblock_s1200.blk ../ABMC/thermal2_rblock_s1200.bcol $TOL $MIT
done