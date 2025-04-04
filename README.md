﻿# GO/PEEK/iC Molecular Dynamics Simulation

This project contains scripts and input files for conducting molecular dynamics simulations of graphene oxide (GO) surface segregation in poly(ether ether ketone) (PEEK) and the grafting of iota-carrageenan (iC) onto the GO/PEEK surface.

## Files Included

1. `input.gograft` - LAMMPS input script for running the simulation
2. `generate_data.py` - Python script to generate the LAMMPS data file with the molecular structures
3. `analyze_simulation.py` - Python script to analyze the simulation results
4. `run_simulation.sh` - Shell script to run the complete workflow

## Requirements

- LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator)
- Python 3 with the following packages:
  - NumPy
  - Matplotlib
  - SciPy

## Usage

1. Ensure all prerequisites are installed
2. Make the shell script executable:
   ```
   chmod +x run_simulation.sh
   ```
3. Run the simulation:
   ```
   ./run_simulation.sh
   ```

The script will:
1. Generate the data file with molecular structures
2. Run the LAMMPS simulation
3. Analyze the results and generate plots
4. Create a summary report

## Simulation Details

The simulation follows these steps:

1. System setup with PEEK chains and randomly distributed GO sheets
2. Equilibration at 380°C (processing temperature)
3. Production run to observe GO segregation
4. Cooling to room temperature
5. Surface activation and iC grafting
6. Analysis of the final structure

## Customization

You can modify simulation parameters by editing:

- `generate_data.py`: Change molecular parameters, system size, concentrations
- `input.gograft`: Modify simulation conditions, timescales, temperature

## Analysis

The analysis script produces several plots:

- Density profiles
- Radial distribution functions
- Gyration radius over time
- GO sheet orientation distribution
- End-to-end distance distributions

All plots and raw data are saved in the `simulation_results` directory.
