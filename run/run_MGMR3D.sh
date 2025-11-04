#!/bin/bash

FIRST_RUN=$1
LAST_RUN=$2
OUTPUT_FILE_PATH=$3

cd /user/vitaldehenau/Python_venvs/jupyter_lab_9/bin
source activate

cd /user/vitaldehenau/Github/pyMGMR3D/run/

if [ -z "$FIRST_RUN" ] || [ -z "$LAST_RUN" ]; then
    echo "Usage: $0 FIRST_RUN LAST_RUN"
    exit 1
fi

for ((RUN=$FIRST_RUN; RUN<=$LAST_RUN; RUN++)); do
    FILENAME=$OUTPUT_FILE_PATH/$(printf "SIM%06d.in" $RUN)
    if [ -f "$FILENAME" ]; then
        echo "Running MGMR3D.sh for $FILENAME..."
        ./pyMGMR3D.sh "$FILENAME"
        mv $(printf "*SIM%06d.*" $RUN) "$OUTPUT_FILE_PATH"
    else
        echo "Warning: $FILENAME not found, skipping."
    fi
done
