#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np
import h5py
import sys
import re

def initialize_mgmr_hdf5(output_path):
    """
    Initializes the HDF5 structure for storing MGMR output.

    Parameters:
        output_path (str or Path): Path to the output .hdf5 file to create.
        input_dict (dict, optional): Dictionary of input parameters to store as attributes under '/inputs'.

    Returns:
        h5py.File object (open in 'a' mode).
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists():
        response = input(f"Warning: {output_path} already exists. Overwrite? [y/n]: ").strip().lower()
        if response != 'y':
            print("Aborting.")
            sys.exit(1)
        else:
            output_path.unlink()  # Delete the old file

    hdf = h5py.File(output_path, 'w')
    
    # Create groups
    hdf.create_group('inputs')            # For run configuration parameters
    hdf.create_group('observers')                     # From FitResult.dat
    hdf.create_group('atmosphere')                  # From sh_Current.dat

    print(f"HDF5 file created at: {output_path.resolve()}")
    return hdf

def parse_mgmr_infile(file_path):
    """
    Extracts input parameters from a MGMR `.in` file.
    
    Parameters:
        file_path (str): Path to the `.in` file.

    Returns:
        dict: Dictionary of input parameters as key-value pairs.
    """
    inputs = {}
    in_block = False

    with open(file_path, 'r') as f:
        for line in f:
            # Remove comments and trim whitespace
            line = re.split(r'!', line)[0].strip()
            if not line:
                continue

            # Detect start and end of &ShPars block
            if '&ShPars' in line:
                in_block = True
                line = line.replace('&ShPars', '').strip()
            elif '&end' in line:
                in_block = False
                line = line.replace('&end', '').strip()

            # Skip if no content after block directives
            if not line and (in_block or not in_block):
                continue

            # Extract all key=value pairs in the line
            parts = [p.strip() for p in line.split(',') if '=' in p]
            for part in parts:
                try:
                    key, val = part.split('=')
                    key = key.strip()
                    val = val.strip()

                    # Remove surrounding quotes for strings
                    if val.startswith('"') and val.endswith('"'):
                        val = val.strip('"')
                        inputs[key] = val
                    elif val.lower() in ['.true.', '.false.']:
                        inputs[key] = val.lower() == '.true.'
                    elif '.' in val or 'e' in val.lower():
                        inputs[key] = float(val)
                    else:
                        inputs[key] = int(val)
                except ValueError:
                    pass  # Silently skip malformed

    return inputs

def write_inputs_to_hdf5(hdf5_path, inputs_dict):
    """
    Writes the MGMR input parameters as attributes under the `/inputs` group in an HDF5 file.

    Parameters:
        hdf5_path (str): Path to the HDF5 file to write or modify.
        inputs_dict (dict): Dictionary of MGMR simulation input parameters.
    """
    with h5py.File(hdf5_path, "a") as f:
        # Create the /inputs group if it doesn't exist
        if "inputs" not in f:
            inputs_group = f.create_group("inputs")
        else:
            inputs_group = f["inputs"]

        # Store each input parameter as an attribute
        for key, value in inputs_dict.items():
            inputs_group.attrs[key] = value

    print(f"Stored {len(inputs_dict)} input parameters under '/inputs' in {hdf5_path}")

def add_observers_to_hdf5(hdf5_path, fitresult_path):
    """
    Reads FitResult.dat and stores Stokes parameters per antenna under /observers
    Each antenna is stored as a group containing one dataset with all Stokes-related values.

    Parameters:
        hdf5_path (str): Path to the HDF5 file to update.
        fitresult_path (str): Path to the FitResult.dat file.
    """
    # Define all columns for clarity
    column_names = [
        "R", "phi", 
        "I", "I_model", "I_sigma", 
        "Q/I", "Q/I_model", "Q/I_sigma", 
        "U/I", "U/I_model", "U/I_sigma", 
        "V/I", "V/I_model", "V/I_sigma", 
        "x", "y", 
        "Psi"
    ]

    df = pd.read_csv(fitresult_path, sep=r'\s+', comment='!', names=column_names, engine='python')

    with h5py.File(hdf5_path, "a") as f:
        fit_grp = f.require_group("observers")

        for idx, row in df.iterrows():
            ant_grp = fit_grp.create_group(f"antenna_{idx:03d}")

            # Store antenna metadata as attributes
            for attr in ["R", "phi", "x", "y", "Psi"]:
                ant_grp.attrs[attr] = float(row[attr])

            # Collect all relevant data into a single array
            stokes_data = np.array([
                row["I"], row["I_model"], row["I_sigma"],
                row["Q/I"], row["Q/I_model"], row["Q/I_sigma"],
                row["U/I"], row["U/I_model"], row["U/I_sigma"],
                row["V/I"], row["V/I_model"], row["V/I_sigma"]
            ], dtype='f8')

            dset = ant_grp.create_dataset("stokes", data=stokes_data)
            dset.attrs["columns"] = np.array([
                "I", "I_model", "I_sigma",
                "Q/I", "Q/I_model", "Q/I_sigma",
                "U/I", "U/I_model", "U/I_sigma",
                "V/I", "V/I_model", "V/I_sigma"
            ], dtype='S')  # Store as byte strings

    print(f"Added {len(df)} antenna entries to '/observers' in '{hdf5_path}'")

def add_atmosphere_to_hdf5(hdf5_path, sh_current_path):
    """
    Reads sh_Current.dat and stores the shower profile in the HDF5 file under /atmosphere,
    grouped into Geometry, Currents, and Force.

    Parameters:
        hdf5_path (str): Path to the HDF5 file.
        sh_current_path (str): Path to the sh_Current.dat file.
    """

    # Define all columns from the .dat file
    columns = [
        "z_km", "X_gcm2", "refractivity", "Ix", "Iy", "charge_excess",
        "dxi", "alpha_tr", "Fx", "Fy", "F_mag_keVpm", "phi"
    ]

    df = pd.read_csv(sh_current_path, sep=r'\s+', comment='!', names=columns, engine='python')

    with h5py.File(hdf5_path, "a") as f:
        atmosphere_group = f.require_group("atmosphere")

        # ---- Geometry group ----
        geometry_cols = ["z_km", "X_gcm2", "refractivity", "alpha_tr"]
        geom_data = df[geometry_cols].to_numpy(dtype='f8')
        dset_geom = atmosphere_group.create_dataset("Geometry", data=geom_data)
        dset_geom.attrs["columns"] = np.array(geometry_cols, dtype='S')

        # ---- Currents group ----
        current_cols = ["Ix", "Iy", "charge_excess", "dxi"]
        current_data = df[current_cols].to_numpy(dtype='f8')
        dset_curr = atmosphere_group.create_dataset("Currents", data=current_data)
        dset_curr.attrs["columns"] = np.array(current_cols, dtype='S')

        # ---- Force group ----
        force_cols = ["Fx", "Fy", "F_mag_keVpm", "phi"]
        force_data = df[force_cols].to_numpy(dtype='f8')
        dset_force = atmosphere_group.create_dataset("Force", data=force_data)
        dset_force.attrs["columns"] = np.array(force_cols, dtype='S')

    print(f"Stored atmosphere data in grouped form under '/atmosphere'. Rows: {len(df)}.")

def add_timetraces_to_observers(hdf5_path, trace_dir):

    trace_dir = Path(trace_dir)
    pattern = re.compile(r"gridttrace-(\d+)-(\d+)\.csv")

    with h5py.File(hdf5_path, "a") as f:
        observer_grp = f.require_group("observers")
        for trace_file in trace_dir.glob("gridttrace-*.csv"):

            match = pattern.match(trace_file.name)
            if not match:
                print(f"Skipping unexpected file format: {trace_file.name}")
                continue

            d, theta = match.groups()
            antenna_name = f"pos_{d}_{theta}"

            df = pd.read_csv(trace_file, sep=r'\s*,\s*', engine='python', comment='!', header=None)
            df.columns = ["t_us", "Re_Ex", "Re_Ey", "Im_Ex", "Im_Ey"]
            df[["Re_Ex", "Re_Ey", "Im_Ex", "Im_Ey"]]*=1e-5 #convert 100kV/m to V/m 
            data = df.to_numpy(dtype='f8')
            dset = observer_grp.create_dataset(f"{antenna_name}", data=data)
            x, y = np.cos(np.deg2rad(float(theta))) * int(d), np.sin(np.deg2rad(float(theta))) * int(d)
            dset.attrs['position'] = np.array([x, y], dtype=float)
            dset.attrs["columns"] = np.array(df.columns, dtype='S')
            dset.attrs["units"] = np.array(["nano s", "V/m", "V/m", "V/m", "V/m"], dtype='S')

    print(f"Stored time traces for all antennas in: '{hdf5_path}'")

if __name__ == "__main__":
    
    hdf5_path = Path(sys.argv[1])
    file_path_inputs = Path(sys.argv[2])
    fitresult_path = Path(sys.argv[3])
    sh_current_path = Path(sys.argv[4])
    trace_dir = Path(sys.argv[5])
    
    if not file_path_inputs.exists():
        print(f"Error: File '{file_path_inputs}' not found.")
        sys.exit(1)
        
    if not fitresult_path.exists():
        print(f"Error: File '{fitresult_path}' not found.")
        sys.exit(1)

    if not sh_current_path.exists():
        print(f"Error: File '{fitresult_path}' not found.")
        sys.exit(1)

    if not trace_dir.exists():
        print(f"Error: File '{trace_dir}' not found.")
        sys.exit(1)
        
    initialize_mgmr_hdf5(hdf5_path)
    
    inputs_dict = parse_mgmr_infile(file_path_inputs)
    write_inputs_to_hdf5(hdf5_path, inputs_dict)
            
    #add_observers_to_hdf5(hdf5_path, fitresult_path)
    add_timetraces_to_observers(hdf5_path, trace_dir)
    add_atmosphere_to_hdf5(hdf5_path, sh_current_path)

    