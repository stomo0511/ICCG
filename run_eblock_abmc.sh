#!/bin/bash

CMD="./llt_bmcic"

#MATRIX="af_shell3 ecology2 Emilia_923 G3_circuit Geo_1438 Hook_1498 thermal2 tmt_sym"
MATRIX="Hook_1498 thermal2 tmt_sym"
MATRIX_DIR="../0_florida"
BLOCK_DIR="../0_blocks/ABMC"

BS=(32 64 128 256 512 1024)
POLYC=(1 2 3)
COLS=(1 2 3)

for matname in ${MATRIX[@]}; do
    file="$MATRIX_DIR/$matname.mtx"
    mat=$(basename "$file" .mtx)
    if [ -f "$file" ]; then
        for bs in "${BS[@]}"; do
			for p in "${POLYC[@]}"; do
                for c in "${COLS[@]}"; do
                    BLK_FILE="$(basename "$file" .mtx)_abmc_B${bs}_p${p}_c${c}.blk"
                    COL_FILE="$(basename "$file" .mtx)_abmc_B${bs}_p${p}_c${c}.bcol"

                    $CMD $file $BLOCK_DIR/$BLK_FILE $BLOCK_DIR/$COL_FILE
                done
            done
        done
    fi
done
