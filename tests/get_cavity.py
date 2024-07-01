import numpy as np
import MDAnalysis as mda
import yaml
import matplotlib.pyplot as plt
import pandas as pd
from ase.io import read, write
from multiprocessing import Pool
import os
from functools import partial
import ctypes
import time

T1 = time.time()
# Load the compiled shared library
lib = ctypes.CDLL('./cavity.so')

# Define C function parameter and return types
lib.compute_for_z.argtypes = [
    ctypes.POINTER(ctypes.c_double), ctypes.c_int, 
    ctypes.POINTER(ctypes.c_double), ctypes.c_int, 
    ctypes.POINTER(ctypes.c_double), ctypes.c_double, 
    ctypes.POINTER(ctypes.c_double)
]

def compute_for_z_python(i, lattices, universe, selected_atoms, x, y, rcutoff, k_B, T):
    """
    Compute the cavity formation energy for a specific z coordinate.

    Parameters:
    i (float): The z coordinate value.
    lattices (numpy.ndarray): The lattice dimensions of the simulation box.
    universe (MDAnalysis.Universe): The MDAnalysis universe containing the trajectory data.
    selected_atoms (MDAnalysis.core.groups.AtomGroup): The selected atoms for the computation.
    x (numpy.ndarray): The x coordinates grid.
    y (numpy.ndarray): The y coordinates grid.
    rcutoff (float): The cutoff distance for interactions.
    k_B (float): Boltzmann constant.
    T (float): Temperature in Kelvin.

    Returns:
    float: The change in chemical potential due to cavity formation.
    """
    # Create a mesh grid for x, y, and the given z slice
    xx, yy, zz = np.meshgrid(x, y, [i])
    centers = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    
    N_values = []
    
    # Convert centers and lattice arrays to C-compatible types
    centers_c = centers.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    box_c = lattices.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # Iterate over each timestep in the trajectory
    for ts in universe.trajectory:
        positions = selected_atoms.positions.astype(np.float64)
        positions_c = positions.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        N_values_temp = np.zeros(len(centers), dtype=np.float64)
        N_values_temp_c = N_values_temp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        # Call the C function
        lib.compute_for_z(centers_c, len(centers), positions_c, len(selected_atoms.positions), box_c, rcutoff, N_values_temp_c)
        
        N_values.extend(N_values_temp)
    
    # Calculate mean and variance of N values
    if N_values:
        average_N = np.mean(N_values)
        variance_N = np.var(N_values)
    else:
        average_N = 0
        variance_N = 0

    # Calculate probability P(v=0) and delta mu for the cavity
    if variance_N > 0:
        P_v_0 = 1 / np.sqrt(2 * np.pi * variance_N) * np.exp(-(average_N ** 2) / (2 * variance_N))
    else:
        P_v_0 = 0

    if P_v_0 > 0:
        delta_mu_v = -k_B * T * np.log(P_v_0)
    else:
        delta_mu_v = float('nan')

    return delta_mu_v

def compute_cavity_formation_energy(universe, volume_size, selection_str, cutoff, temperature, z_range, step_size=1.0):
    """
    Compute the cavity formation energy across a range of z coordinates.

    Parameters:
    universe (MDAnalysis.Universe): The MDAnalysis universe containing the trajectory data.
    volume_size (list): The dimensions of the simulation box [x, y, z].
    selection_str (str): The selection string for selecting atoms.
    cutoff (float): The cutoff distance for interactions.
    temperature (float): Temperature in Kelvin.
    z_range (numpy.ndarray): The range of z coordinates to compute over.
    step_size (float): The step size for the grid in x and y directions.

    Returns:
    tuple: The z coordinates and corresponding energies.
    """
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    T = temperature  # Temperature in Kelvin
    
    # Select atoms based on the selection string
    selected_atoms = universe.select_atoms(selection_str)
    
    # Define the grid in the x and y directions based on step_size
    num_steps_x = int(volume_size[0] // step_size)
    num_steps_y = int(volume_size[1] // step_size)
    x, y = np.linspace(0, volume_size[0], num_steps_x, endpoint=False), np.linspace(0, volume_size[1], num_steps_y, endpoint=False)
    
    rcutoff = cutoff

    # Get the number of CPU cores available
    num_cores = os.cpu_count()
    print(f"Number of CPU cores: {num_cores}")

    # Use multiprocessing pool for parallel computation
    with Pool(processes=num_cores) as pool:
        partial_compute = partial(compute_for_z_python, lattices=volume_size, universe=universe, selected_atoms=selected_atoms, x=x, y=y, rcutoff=rcutoff, k_B=k_B, T=T)
        energies = pool.map(partial_compute, z_range)

    return z_range, energies

def read_config(config_file):
    """
    Load configuration from a YAML file.

    Parameters:
    config_file (str): Path to the configuration file.

    Returns:
    dict: The loaded configuration.
    """
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

def save_results_to_txt(output_file, results):
    """
    Save results to a text file.

    Parameters:
    output_file (str): Path to the output text file.
    results (list): The list of results to save.
    """
    with open(output_file + '.txt', 'w') as f:
        for result in results:
            f.write(f'File: {result["file"]}\n')
            f.write('z,Energy\n')
            for z, energy in zip(result["z"], result["energies"]):
                f.write(f'{z},{energy}\n')
            f.write('\n')

def save_results_to_excel(output_file, results):
    """
    Save results to an Excel file.

    Parameters:
    output_file (str): Path to the output Excel file.
    results (list): The list of results to save.
    """
    writer = pd.ExcelWriter(output_file + '.xlsx', engine='xlsxwriter')
    for result in results:
        df = pd.DataFrame({'z': result["z"], 'Energy': result["energies"]})
        df.to_excel(writer, sheet_name=result["file"], index=False)
    writer.close()  # Use close() instead of save()

def plot_results(results, output_file):
    """
    Plot the results and save as PNG.

    Parameters:
    results (list): The list of results to plot.
    output_file (str): Path to the output image file.
    """
    plt.figure(figsize=(10, 6))
    for result in results:
        plt.plot(result["z"], result["energies"], label=result["file"])
    plt.xlabel('z')
    plt.ylabel('Cavity Formation Energy (J)')
    plt.title('Cavity Formation Energy vs. z')
    plt.legend()
    plt.savefig(output_file + '.png')
    plt.show()

def preprocess_xyz_file_with_ase(xyz_file, temp_xyz_file):
    """
    Preprocess an XYZ file using ASE and extract volume information.

    Parameters:
    xyz_file (str): Path to the input XYZ file.
    temp_xyz_file (str): Path to the temporary XYZ file for MDAnalysis.

    Returns:
    list: The dimensions of the simulation box [x, y, z].
    """
    # Read XYZ file using ASE
    atoms_list = read(xyz_file, ':')
    
    # Write to a temporary file to ensure MDAnalysis can read it
    write(temp_xyz_file, atoms_list)
    
    # Extract volume information
    volume_size = atoms_list[0].cell.cellpar()[:3]
    
    return volume_size

# Load configuration from YAML file
config = read_config('config.yaml')

# Process each XYZ file
results = []
for xyz_file in config['xyz_files']:
    try:
        # Preprocess XYZ file using ASE
        temp_xyz_file = 'temp.xyz'
        volume_size = preprocess_xyz_file_with_ase(xyz_file, temp_xyz_file)
        
        # Load the preprocessed XYZ file with MDAnalysis
        universe = mda.Universe(temp_xyz_file, permissive=True)
        
        # Parameters from the config file or use default values
        selection_str = config.get('selection_str', 'name O')  # Use the correct selection string
        cutoff = config.get('cutoff', 2.0)
        temperature = config.get('temperature', 298.15)
        z_range_config = config.get('z_range', {'start': 3.4, 'end': 15.4, 'step': 0.25})
        z_range = np.arange(z_range_config['start'], z_range_config['end'], z_range_config['step'])
        step_size = config.get('step_size', 1.0)

        # Compute the cavity formation energy for each z interval
        z_values, energies = compute_cavity_formation_energy(universe, volume_size, selection_str, cutoff, temperature, z_range, step_size)
        
        # Store results
        results.append({
            "file": xyz_file,
            "z": z_values,
            "energies": energies
        })
    except Exception as e:
        print(f"Error processing file {xyz_file}: {e}")

T2 = time.time()
print(T2-T1)

# Output file name
output_file = config.get('output_file', 'results')

# Save results to txt and Excel
save_results_to_txt(output_file, results)
save_results_to_excel(output_file, results)

# Plot results
plot_results(results, output_file)
