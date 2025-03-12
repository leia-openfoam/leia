#!/usr/bin/env python3

import os
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rcParams

# Set default font and figure size for plots 
rcParams.update({'font.size': 11, 'figure.figsize': (6.4, 4.8)})  


def read_parameters_from_file(file_path):
    """
    Reads the parameters from the PyFoamPrepareCaseParameters file.
    """
    params = {}
    with open(file_path, "r") as file:
        for line in file:
            match = re.match(r"(\w+)\s+(\S+);", line)
            if match:
                params[match.group(1)] = match.group(2)
    return params

def agglomerate_data_files(folder_pattern, csv_file_name):
    """
    Agglomerates data from CSV files in folders matching the pattern, 
    extending them with parameters from PyFoamPrepareCaseParameters.
    """
    # Find folders matching the pattern and sort them
    folder_names = [f for f in os.listdir('.') if os.path.isdir(f) and folder_pattern in f]
    folder_names.sort()

    global_dataframe = pd.DataFrame()

    for folder in folder_names:
        parameters_file = os.path.join(folder, "PyFoamPrepareCaseParameters")
        if not os.path.exists(parameters_file):
            print(f"Warning: Parameters file not found in folder {folder}")
            continue

        # Read parameters from PyFoamPrepareCaseParameters
        parameters = read_parameters_from_file(parameters_file)

        # Read the CSV file in the folder
        data_file_path = os.path.join(folder, csv_file_name)
        if not os.path.exists(data_file_path):
            print(f"Warning: CSV file '{csv_file_name}' not found in folder {folder}")
            continue

        df = pd.read_csv(data_file_path)

        # Extend the DataFrame with additional columns
        #for key in ["DIV_SCHEME", "N_CELLS", "solver", "N_NON_ORTH", "NU"]:
        #    df[key] = parameters.get(key, None)

        # Extend the DataFrame with all additional columns 
        # - Use all parameter variation parameters, ignore OpenFOAM metadata
        # - OpenFOAM metadata is ignored because the CSV files must still be
        #   archived alongside the whole input setup for OpenFOAM cases, all
        #   solution parameters must be archived along the CSV otherwise the
        #   results are not reproducible.
       
        exclude_keys = {"casePath", "caseName", "foamFork", 
                        "numberOfProcessors", "foamVersion"}

        # Assign only the allowed keys
        for key in parameters.keys():
                if key not in exclude_keys:
                            df[key] = parameters.get(key, None)

        # Append to the global DataFrame
        global_dataframe = pd.concat([global_dataframe, df], ignore_index=True)

    return global_dataframe 

def plot_convergence_rate(global_dataframe):
    """
    Plots UerrLinf, UerrL1, UerrL2 over h for each unique DIV_SCHEME in the global dataframe.
    Annotates each plot with the convergence rate at the midpoint of the linear fit.

    Parameters:
        global_dataframe (pd.DataFrame): The dataframe containing the aggregated data.
    """
    # Get unique NU values
    nus = global_dataframe["NU"].unique()

    for nu in nus:

        # Get unique DIV_SCHEME values
        div_schemes = global_dataframe["DIV_SCHEME"].unique()

        for div_scheme in div_schemes:

            # Filter the dataframe for the current div_scheme
            is_divscheme = global_dataframe["DIV_SCHEME"] == div_scheme
            is_nu = global_dataframe["NU"] == nu 
            scheme_dataframe = global_dataframe[is_divscheme & is_nu]
            print (scheme_dataframe)

            # Sort by 'h' to ensure proper plotting
            scheme_dataframe = scheme_dataframe.sort_values(by="h")

            # Create a log-log plot for each error type
            plt.figure(figsize=(10, 6))

            for error_type in ["UerrLinf", "UerrL1", "UerrL2"]:
                # Plot data
                h = scheme_dataframe["h"]
                error = scheme_dataframe[error_type]
                print(h)
                print(error)
                plt.loglog(h, error, marker='o', label=f"{error_type}")

                # Fit a linear line on log-log scale
                log_h = np.log(h)
                log_error = np.log(error)
                slope, intercept = np.polyfit(log_h, log_error, 1)

                # Generate the linear fit line
                fit_line = np.exp(intercept) * h**slope
                plt.loglog(h, fit_line, linestyle='--', label=f"{error_type} Fit")
                # global_dataframe[f"{error_type}Convergence"] = slope

                # Annotate the plot with the slope at the midpoint of the fit line
                mid_x = np.sqrt(h.iloc[0] * h.iloc[-1])  # Geometric midpoint on the x-axis
                mid_y = np.exp(intercept) * mid_x**slope
                plt.annotate(f"convergence rate: {slope:.2f}", xy=(mid_x, mid_y),
                             xytext=(10, 10), textcoords="offset points", fontsize=10)

            # Configure plot
            plt.title(f"Convergence Plot for DIV_SCHEME: {div_scheme}, NU: {nu}")
            plt.xlabel("h")
            plt.ylabel("Error")
            plt.legend()
            plt.grid(True, which="both", linestyle="--", linewidth=0.5)
            plt.tight_layout()

            # Save or show the plot
            plt.savefig(f"convergence-rate-{div_scheme}-nu-{nu}.pdf")
            plt.show()

def plot_elapsed_cpu_time(global_dataframe, folder_pattern):
    """
    Creates a column plot for ELAPSED_CPU_TIME(h) grouped by DIV_SCHEME.

    Parameters:
        global_dataframe (pd.DataFrame): The dataframe containing the aggregated data.
    """

    plt.clf()

    # Get unique h values and DIV_SCHEME values
    unique_h = sorted(global_dataframe["h"].unique())
    div_schemes = global_dataframe["DIV_SCHEME"].unique()
    nus = global_dataframe["NU"].unique()

    bar_width = 0.2
    colors = plt.cm.tab10.colors  # Use a colormap for better differentiation

    for i, h in enumerate(unique_h):
        h_dataframe = global_dataframe[global_dataframe["h"] == h]

        for j, div_scheme in enumerate(div_schemes):
            # Filter by h and DIV_SCHEME
            scheme_dataframe = h_dataframe[h_dataframe["DIV_SCHEME"] == div_scheme]

            if not scheme_dataframe.empty:
                # Get ELAPSED_CPU_TIME value
                cpu_time = scheme_dataframe["ELAPSED_CPU_TIME"].values[0]

                # Calculate bar position
                bar_position = i + j * bar_width

                # Plot bar
                plt.bar(bar_position, cpu_time, bar_width, 
                        #label=div_scheme if h == unique_h[0] else "", 
                        color=colors[j % len(colors)])

    # Create a separate legend using div_schemes and colors
    legend_handles = [plt.Line2D([0], [0], color=colors[i % len(colors)], lw=4) for i in range(len(div_schemes))]
    plt.legend(legend_handles, div_schemes, title="DIV_SCHEME")

    # Configure x-axis with h labels
    bar_positions = np.arange(len(unique_h)) + (len(div_schemes) - 1) * bar_width / 2
    plt.xticks(bar_positions, unique_h, rotation=45, ha='right')

    # Configure the plot
    plt.title("Elapsed CPU time by h for each DIV_SCHEME")
    #plt.semilogy()
    plt.xlabel("h")
    plt.ylabel("Elapsed CPU Time")
    plt.grid(axis='y', linestyle="--", linewidth=0.5)
    plt.tight_layout()

    # Save and show the plot
    output_file = f"elapsed-cpu-time-plot-{folder_pattern}.pdf"
    plt.savefig(output_file)
    print(f"Plot saved as {output_file}")
    plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Agglomerate data files with extended parameters.")
    parser.add_argument("--folder-pattern", type=str, required=True, help="Regex pattern to match folder names.")
    parser.add_argument("--data-file", type=str, required=True, help="Name of the CSV file to read in each folder.")
    
    args = parser.parse_args()

    global_dataframe = agglomerate_data_files(args.folder_pattern, args.data_file)

    plot_convergence_rate(global_dataframe)

    plot_elapsed_cpu_time(global_dataframe, args.folder_pattern)

    # Save the global DataFrame as a CSV file
    output_file_name = args.folder_pattern.rstrip('0').rstrip("_")
    output_file = f"{output_file_name}.csv"
    global_dataframe.to_csv(output_file, index=False)
    print(f"Global data saved to {output_file}")
    
