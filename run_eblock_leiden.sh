#!/bin/bash

CMD="./llt_bmcic"

MATRIX_DIR="../0_florida"
BLOCK_DIR="../0_blocks/Leiden"

#MATRIX="af_shell3 ecology2 Emilia_923 G3_circuit Geo_1438 Hook_1498 thermal2 tmt_sym"
MATRIX="tmt_sym"
GAMMA=(05 075 1 125 15 2)
COLS=(1 2 3)

for matname in ${MATRIX[@]}; do
    file="$MATRIX_DIR/$matname.mtx"
    mat=$(basename "$file" .mtx)
    if [ -f "$file" ]; then
        for g in "${GAMMA[@]}"; do
            for c in "${COLS[@]}"; do
                BLK_FILE="$(basename "$file" .mtx)_leiden_cpm_r0p00${g}_c${c}.blk"
                COL_FILE="$(basename "$file" .mtx)_leiden_cpm_r0p00${g}_c${c}.bcol"

                $CMD $file $BLOCK_DIR/$BLK_FILE $BLOCK_DIR/$COL_FILE
            done
        done
    fi
done
