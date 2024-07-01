# Cavity Analysis Project

This project is designed to analyze cavity formation energies in molecular dynamics simulations. It includes a Python script to process data, a configuration file for settings, and a C program for computational efficiency.

## Requirements

- Python 3.x
- Poetry (for dependency management)
- GCC (for compiling the C code)
- Libraries: `numpy`, `MDAnalysis`, `matplotlib`, `pandas`, `ase`, `pyyaml`

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/WanluLigroupUCSD/Cavity.git
   cd Cavity
   ```

2. **Install dependencies using Poetry:**

    ```bash
    poetry install
    ```

3. **Compile the C code:**
    ```bash
    gcc -shared -o cavity.so -fPIC cavity.c -lm
    ```
## Configuration

The configuration settings for the project are defined in `config.yaml`, which you can find in the `tests` directory.

## Usage

1. **Compile the `cavity.c` with gcc or icc:**

    Make sure the `cavity.so` file in the directory that you will run the python script.

2. **Run the python script:**

    ```bash
    poetry run python get_cavity.py
    ```

    The script performs the following tasks:

    -   Loads the configuration from config.yaml
    -   Initializes MDAnalysis with the specified XYZ files
    -   Computes cavity formation energies using the compiled C library (cavity.so)
    -   Outputs results to the specified file

# Project Structure

- **get_cavity.py:** Main Python script for data processing and analysis
- **config.yaml:** Configuration file with parameters for the analysis
- **cavity.c:** C code for computational routines to speed up distance calculations

# License
This project is licensed under the MIT License. See the `LICENSE` file for details.

# Acknowledgements
This project is developed as part of the research work at the University of California, San Diego, under the supervision of Dr. Li.