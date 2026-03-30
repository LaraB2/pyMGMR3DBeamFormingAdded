from radiotools.atmosphere import models as atm
import numpy as np
import sys
import re

def export_for_mgmr(model_id, atm_models, rh0):
    model = atm_models[model_id]
    
    # Extract the arrays
    a_vals = model['a']
    b_vals = model['b']
    c_vals = model['c']
    h_tops = model['h']  # These are the boundaries (km to m already done in your dict)

    filename = 'current_atm.dat'
    with open(filename, 'w') as f:
        # 1. Write rh0
        f.write(f"{rh0}\n")
        
        # 2. Handle the Sea Level Buffer
        # If the first h_top is > 0, we prepend a layer starting at 0.0
        starts_above_zero = h_tops[0] > 0
        num_layers = len(a_vals)
        f.write(f"{num_layers}\n")
        
        # 3. Write layers
        # Note: Your 'h' array usually has 4 values for 5 layers (the last is implicit infinity)
        for i in range(num_layers):
            a = a_vals[i]*1e-4
            b = b_vals[i]*1e-4
            c = c_vals[i]
            # Use the boundary from 'h', or a very high number for the last layer
            h_boundary = h_tops[i] if i < len(h_tops) else 1.e9
            
            # If it's the very first layer and it starts above 0, 
            # we should technically define the bottom as 0 in our Fortran reader.
            # But for now, just exporting the boundaries is enough.
            f.write(f"{a} {b} {c} {h_boundary}\n")
            
    print(f"Atmosphere model {model_id} exported to {filename}")

def get_val(text, key):
    # This regex looks for: key = value, but ignores everything after a '!'
    # It handles cases like 'atm_model_id = 15, ! comment'
    pattern = rf"^\s*{key}\s*=\s*([^,! \n]+)"
    match = re.search(pattern, text, re.IGNORECASE | re.MULTILINE)
    if match:
        return match.group(1).strip().replace("'", "").replace('"', "")
    return None

if __name__ == "__main__":
    
   with open(sys.argv[1], 'r') as f:
        content = f.read()

        # Find 'atm_model_id = X' in the text
        model_id_str = get_val(content, "atm_model_id")
        model_id = int(model_id_str) if model_id_str else 1
        
        # Find 'rh0 = X'
        rh0_str = get_val(content, "rh0")
        rh0 = float(rh0_str) if rh0_str else 0.000292
    
        if model_id in atm.atm_models:
            print(f"Selecting existing model: {model_id}")
            export_for_mgmr(model_id, atm.atm_models, rh0)
        else:
            print(f"Error: Model ID {model_id} not found in dictionary!")
            sys.exit(1)