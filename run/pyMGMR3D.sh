#!/bin/bash

# --- Validate input ---
# Checks if at least one argument (the input file) is provided.
# If not, it prints a usage message and exits with an error.
if [ $# -lt 1 ]; then
    echo "Usage: $0 <inputfile.in> [--no-hdf5] [--plot <plot_output_directory>] [--result-dir <final_results_directory>]"
    exit 1
fi

# Assign the first argument as the input file.
# Shift arguments so that $1 now refers to the next argument (flags).
INPUT_FILE=$1
shift

# Extract the base name of the input file (without extension).
# Get the current working directory, which is considered the run folder.
# Determine the program directory by going up one level from the run folder.
BASENAME=$(basename "$INPUT_FILE" .in)
RUN_FOLDER=$(pwd)
PROG_DIR=$(dirname "$RUN_FOLDER")/program

# Initialize flags for plotting and HDF5 conversion as enabled by default.
# Initialize the directory where plots will be directly saved. This is the output for --plot.
# Initialize the target directory for moving the temporary 'plot' folder (results).
# Initialize HDF5 filename. By default, it will be the same as the input file with .hdf5 extension.
PLOT_ENABLED=false
HDF5_ENABLED=true
PLOT_OUTPUT_DIR=""
RESULT_DIR=""
HDF5_FILENAME="${RUN_FOLDER}/${BASENAME}.hdf5"

# --- Parse flags ---
# Loop through all provided command-line arguments to parse flags.
while [[ $# -gt 0 ]]; do
    case $1 in
        --plot)
            PLOT_ENABLED=true
            if [ -z "$2" ] || [[ "$2" == --* ]]; then
                echo "Error: --plot requires a directory argument."
                exit 1
            fi
            PLOT_OUTPUT_DIR=$2
            shift 2
            ;;
        --no-hdf5)
            HDF5_ENABLED=false
            shift
            ;;
        --result-dir)
            if [ -z "$2" ] || [[ "$2" == --* ]]; then
                echo "Error: --result-dir requires a directory argument."
                exit 1
            fi
            RESULT_DIR=$2
            shift 2
            ;;
        *)
            echo "Warning: Unrecognized argument '$1'. Ignoring."
            shift
            ;;
    esac
done

# Create the plot output directory if plotting is enabled and it doesn't exist.
if $PLOT_ENABLED && [ -n "$PLOT_OUTPUT_DIR" ]; then
    mkdir -p "$PLOT_OUTPUT_DIR"
fi

# Print important paths for user information.
echo "Program folder = ${PROG_DIR}"
echo "Run folder     = ${RUN_FOLDER}"
echo "Input file     = ${INPUT_FILE}"
if $PLOT_ENABLED; then
    echo "Plot output to = ${PLOT_OUTPUT_DIR}"
fi
echo "HDF5 output to = ${HDF5_FILENAME}"
if [ -n "$RESULT_DIR" ]; then
    echo "Temporary 'plot' folder will be moved to = ${RESULT_DIR}/${BASENAME}_results"
fi

# --- Compile MGMR ---
# Change directory to the program directory.
# If changing directory fails, exit the script.
# Compile the MGMR program using the specified makefile.
# If compilation fails, exit the script.
# Change back to the original run folder.
# If changing directory fails, exit the script.
cd "${PROG_DIR}" || { echo "Error: Could not change to program directory ${PROG_DIR}. Exiting."; exit 1; }
make -f "MGMR3D_fit-makefile-v5.make" || { echo "Error: MGMR compilation failed. Exiting."; exit 1; }
cd "${RUN_FOLDER}" || { echo "Error: Could not change back to run directory ${RUN_FOLDER}. Exiting."; exit 1; }

# --- Run simulation ---
# Create a temporary 'plot' directory within the run folder for intermediate output.
# Execute the MGMR program, feeding it the input file's content.
mkdir -p "${RUN_FOLDER}/plot"
"${PROG_DIR}/MGMR3D_fit-v5" < "$INPUT_FILE"

# --- Plotting with GLE and Python ---
# Change directory to the temporary 'plot' folder to handle output files.
cd "${RUN_FOLDER}/plot" || { echo "Error: Could not change to temporary plot directory. Exiting."; exit 1; }

# Check if plotting is enabled.
if $PLOT_ENABLED; then
    # Check if GLE (Graphics Layout Engine) command is available.
    if command -v gle &> /dev/null; then
        gle -d pdf -o "${PLOT_OUTPUT_DIR}/${BASENAME}_FitStokes.pdf" "${PROG_DIR}/FitStokes.GLE" "${RUN_FOLDER}/plot/FitResult"
        gle -d jpg -r 200 -o "${PLOT_OUTPUT_DIR}/${BASENAME}_FitStokes-map.jpg" "${PROG_DIR}/FitStokes-map.GLE" "${RUN_FOLDER}/plot/"
        gle -d pdf -o "${PLOT_OUTPUT_DIR}/${BASENAME}_sh-current.pdf" "${PROG_DIR}/sh-current.GLE" "${RUN_FOLDER}/plot/"
    else
        echo "GLE not found. Skipping GLE plots."
    fi
    # Check if python3 command is available and the required Python libraries.
    if command -v python3 &> /dev/null; then
        python3 - <<EOF
try:
    import pandas
    import matplotlib
    import numpy
except ImportError as e:
    print(f"Missing library: {e.name}")
    exit(1)
EOF
        if [ $? -eq 0 ]; then
            python3 "${PROG_DIR}/FitStokes.py" "${RUN_FOLDER}/plot/FitResult.dat" "${PLOT_OUTPUT_DIR}"
            python3 "${PROG_DIR}/FitStokes-map.py" "${RUN_FOLDER}/plot/" "${PLOT_OUTPUT_DIR}"
            python3 "${PROG_DIR}/sh-current.py" "${RUN_FOLDER}/plot/" "${PLOT_OUTPUT_DIR}"
        else
            echo "Skipping Python plotting due to missing libraries."
        fi
    else
        echo "Python3 not installed. Skipping plot generation."
    fi
fi

# --- Convert to HDF5 ---
# Check if HDF5 conversion is enabled and the required Python libraries.
if $HDF5_ENABLED; then
    if command -v python3 &> /dev/null; then
        python3 - <<EOF
try:
    import pandas as pd
    import numpy as np
    import h5py
    from pathlib import Path
    import re
except ImportError as e:
    print(f"Missing library: {e.name}")
    exit(1)
EOF
        if [ $? -eq 0 ]; then
            python3 "${RUN_FOLDER}/MGMR_HDF5.py" "${RUN_FOLDER}/${BASENAME}.hdf5" "${RUN_FOLDER}/${INPUT_FILE}" "${RUN_FOLDER}/plot/FitResult.dat" "${RUN_FOLDER}/plot/sh_Current.dat" "${RUN_FOLDER}/plot/"
        else
            echo "Skipping HDF5 conversion due to missing libraries."
        fi
    else
        echo "Python3 not installed. Skipping HDF5 generation."
    fi
fi

# --- Handle Result Directory and Cleanup ---
# Change back to the original run folder to ensure correct path for mv/rm.
cd "${RUN_FOLDER}" || { echo "Error: Could not change back to run directory ${RUN_FOLDER}. Exiting."; exit 1; }

# If --result-dir was provided, move the 'plot' directory.
if [ -n "$RESULT_DIR" ]; then
    echo "Moving temporary results folder to: ${RESULT_DIR}/${BASENAME}_results"
    mkdir -p "$RESULT_DIR"
    mv "${RUN_FOLDER}/plot" "${RESULT_DIR}/${BASENAME}_results" || \
        { echo "Error: Failed to move temporary 'plot' directory to ${RESULT_DIR}/${BASENAME}_results. Please check permissions or path."; exit 1; }
else
    echo "Cleaning up temporary 'plot' directory."
    rm -r "${RUN_FOLDER}/plot"
fi

echo "Run complete."
exit 0
