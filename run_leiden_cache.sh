#!/bin/bash

export OMP_NUM_THREADS=32

#echo "mat, nb, nc, cnv, time, iters, rel_resid, resid"

# 実行回数
N=5

EVENTS="l2_cache_accesses_from_dc_misses,\
l2_cache_hits_from_dc_misses,\
l2_cache_misses_from_dc_misses"

ABMC="../ABMC/leiden"

CMD="./abmc_crs"
TOL="1.0E-08"
MIT="10000"

TARGET_DIR="../0_florida"

for file in $TARGET_DIR/*.mtx; do
    mat=$(basename "$file" .mtx)
    if [ -f "$file" ]; then
        abmc_out=$($ABMC $file)
        nb_val=$(echo "$abmc_out" | awk '{print $3}')
        nc_val=$(echo "$abmc_out" | awk '{print $6}')
            
        BLK="$(basename "$file" .mtx)_leiden.blk"
        BCOL="$(basename "$file" .mtx)_leiden.bcol"

        for i in $(seq 1 $N); do
            echo -n "$mat, $nb_val, $nc_val, "
            perf stat -e $EVENTS $CMD $file $BLK $BCOL $TOL $MIT
        done

        rm -f $BLK $BCOL
    fi
done
