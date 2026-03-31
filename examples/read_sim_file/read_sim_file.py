"""
Example runner: read a MEGAlib .sim/.sim.gz file and produce per-panel outputs.

Purpose:
- Load the MEGAlib geometry and simulated events
- Build a DataFrame with energy deposits for each ACS/BGO panel and timestamps
- Save the main CSV and a lightweight .lc counts file per panel

Requirements:
- MEGAlib available in ROOT (MEGALIB environment set)
- numpy
- The function `read_sim_file_with_megalib` must be importable in PYTHONPATH

Usage:
    python examples/read_sim_file/read_sim_file.py <geometry_name> <filename> <output_file> <shared> <noise>

Args:
- geometry_name: path to MEGAlib geometry setup file
- filename: path to input .sim or .sim.gz
- output_file: path to output .evt CSV
- shared: "true" to apply ASIC-sharing for Z0_4/Z1_4, otherwise "false"
- noise: "0" to disable energy noising, any other value keeps default noising
"""

import sys
import numpy as np

from cosiburstpy.megalib.read_sim_file_with_megalib import read_sim_file_with_megalib

if __name__ == "__main__":

    geometry_name = sys.argv[1]
    filename = sys.argv[2]
    output_file = sys.argv[3]
    shared = sys.argv[4]
    noise = sys.argv[5]

    df = read_sim_file_with_megalib(filename, geometry_name, shared, noise)
    
    if shared == "true":
        output_file = output_file.replace(".evt","_shared.evt")
    df.to_csv(output_file, index=True)

    # create file with counts for each of the 6 panels

    bgo_data_seconds = {}

    for col in df.columns:
        if col.startswith("bgo_"):
            filtered_data = df[df[col] > 0][[col]]
            times_in_seconds = filtered_data.index.map(lambda x: x.timestamp())
            bgo_data_seconds[col] = np.array(times_in_seconds)


    output_file_lc = open(output_file+".lc","w")

    for col, times_in_seconds in bgo_data_seconds.items():

        output_file_lc.write(col+" "+str(len(times_in_seconds))+"\n")

    output_file_lc.close()

    