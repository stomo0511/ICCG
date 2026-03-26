#!/usr/bin/sh

export OMP_NUM_THREADS=32

echo "mat, cnv, time, iters, rel_resid, resid"

# 実行回数
N=10

CMD="./iccg_crs"
TOL="1.0E-08"
MIT="10000"

TARGET_DIR="../0_florida"

for file in $TARGET_DIR/*.mtx; do
    mat=$(basename "$file" .mtx)
    if [ -f "$file" ]; then
        for i in $(seq 1 $N); do
            echo -n "$mat, "
            $CMD $file $COL $TOL $MIT
        done
    fi
done