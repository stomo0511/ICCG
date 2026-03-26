#!/usr/bin/sh

export OMP_NUM_THREADS=32

TARGET_DIR="../0_florida"

GMC="../ABMC/gmc"

CMD="./piccg_crs"
TOL="1.0E-08"
MIT="10000"

# 実行回数
N=10

echo "mat, nc, cnv, time, iters, rel_resid, resid"

for file in $TARGET_DIR/*.mtx; do
    mat=$(basename "$file" .mtx)
    if [ -f "$file" ]; then
        gmc_out=$($GMC $file)
    fi
    nc_val=$(echo "$gmc_out" | awk '{print $3}')
    COL="$(basename "$file" .mtx)_gmc.col"

    for i in $(seq 1 $N); do
        echo -n "$mat, $nc_val, "
        $CMD $file $COL $TOL $MIT
    done

    rm -f $COL
done
