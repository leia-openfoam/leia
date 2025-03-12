#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile


def load_fv_schemes_entry():
    """Load the convective scheme from the fvSchemes file."""
    fv_schemes_file = ParsedParameterFile("system/fvSchemes")
    div_phi_U = fv_schemes_file["divSchemes"]["div(phi,U)"]
    return div_phi_U.lstrip("$")


def load_latest_sample_data(post_processing_dir):
    """Load the latest time step sample data."""
    latest_time_step = max(
        [int(folder) for folder in os.listdir(post_processing_dir) if folder.isdigit()]
    )
    latest_folder = os.path.join(post_processing_dir, str(latest_time_step))
    sample_file = os.path.join(latest_folder, 'centreLine_p_U.xy')
    data = pd.read_csv(sample_file, sep='\s+', header=None, names=['r', 'p', 'Ux', 'Uy', 'Uz'])
    return data


def analytical_utheta(r, A, B):
    """Calculate analytical velocity."""
    return A * r + B / r


def analytical_p(r, A, B, C):
    """Calculate analytical pressure."""
    return 0.5 * A**2 * r**2 + 2 * A * B * np.log(r) - 0.5 * B**2 / r**2 - C


def load_case_parameters(file_path):
    """Load, clean, and prefix case parameters from the PyFoamPrepareCaseParameters file."""
    keys, values = [], []
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            for line in file:
                if not line.strip().startswith('/'):  # Ignore lines starting with '/'
                    parts = line.strip().split()
                    if len(parts) >= 2:  # Ensure at least a key and a value
                        key = "METADATA_" + parts[0]  # Add "METADATA_" prefix to the key
                        value = parts[1].strip('";')  # Remove wrapping quotes and trailing semicolon
                        keys.append(key)
                        values.append(value)
    return keys, values



def add_parameters_to_dataframe(data, keys, values):
    """Add case parameters to the dataframe."""
    for key, value in zip(keys, values):
        data[key] = value
    return data


def plot_comparison(data, x_col, y_cols, analytical_col, labels, title, output_file):
    """Plot comparison of numerical and analytical results."""
    plt.figure(figsize=(8, 6))
    plt.plot(data[x_col], data[y_cols[0]], label=labels[0], color="red", marker="o", linestyle="")
    plt.plot(data[x_col], data[analytical_col], label=labels[1], color="black", linestyle="-")
    plt.xlabel("Radius, r")
    plt.ylabel(labels[2])
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(output_file, dpi=300)
    plt.show()


def main():
    # Load fvSchemes entry
    div_phi_U = load_fv_schemes_entry()
    print("div(phi,U) entry:", div_phi_U)

    # Load latest sample data
    post_processing_dir = 'postProcessing/sample1'
    data = load_latest_sample_data(post_processing_dir)

    # Analytical solution parameters
    Omega1, Omega2, R1, R2 = 100.0, 0.0, 1.0, 2.0
    mu = Omega2 / Omega1
    A = Omega1 * (1.0 - R2**2 * mu / R1**2) / (1.0 - R2**2 / R1**2)
    B = R1**2 * Omega1 * (1.0 - mu) / (1.0 - R1**2 / R2**2)
    C = 0.5 * A**2 * R1**2 + 2 * A * B * np.log(R1) - 0.5 * B**2 / R1**2

    # Add analytical solutions to dataframe
    data['analytical_p'] = analytical_p(data['r'], A, B, C)
    data['analytical_utheta'] = analytical_utheta(data['r'], A, B)

    # Load case parameters and add them to the dataframe
    file_path = 'PyFoamPrepareCaseParameters'
    keys, values = load_case_parameters(file_path)
    data = add_parameters_to_dataframe(data, keys, values)

    # Save data
    data_folder = "data"
    os.makedirs(data_folder, exist_ok=True)
    data.to_csv(f"{data_folder}/{div_phi_U}_data.csv", index=False)
    print(data)

    # Plot comparisons
    plot_comparison(data, 'r', ['p'], 'analytical_p',
                    [f"{div_phi_U} Pressure", "Analytical Pressure", "Pressure, p"],
                    "Pressure Comparison", f"{data_folder}/p_{div_phi_U}.png")

    plot_comparison(data, 'r', ['Uy'], 'analytical_utheta',
                    [f"{div_phi_U} Utheta", "Analytical Utheta", "Uθ (rad/s)"],
                    "Velocity Comparison", f"{data_folder}/Utheta_{div_phi_U}.png")


if __name__ == "__main__":
    main()

