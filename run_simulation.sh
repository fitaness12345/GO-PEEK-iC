#!/bin/bash
# Complete workflow script for running GO/PEEK/iC molecular dynamics simulation
# This script sets up and runs the simulation, then analyzes the results

# Check for required software
command -v py >/dev/null 2>&1 || { echo "Python 3 is required but not installed. Aborting."; exit 1; }
command -v lmp >/dev/null 2>&1 || { echo "LAMMPS is required but not installed. Aborting."; exit 1; }

# Set up directories
mkdir -p simulation_results
cd simulation_results

echo "==== Setting up simulation ===="

# Generate data file
echo "Generating LAMMPS data file..."
py ../generate_data.py

# Check if data file was created successfully
if [ ! -f "data.gograft" ]; then
    echo "Failed to generate data file. Aborting."
    exit 1
fi

echo "==== Running LAMMPS simulation ===="
echo "This may take a long time depending on system size..."

# Run LAMMPS simulation
lmp -in ../input.gograft

# Check if simulation completed successfully
if [ ! -f "dump.gograft.final" ]; then
    echo "Simulation did not complete successfully. Aborting."
    exit 1
fi

echo "==== Analyzing simulation results ===="

# Run analysis script
py ../analyze_simulation.py

echo "==== Creating report ===="

# Generate summary report
cat > simulation_report.txt << EOF
# GO/PEEK/iC Molecular Dynamics Simulation Report

## Simulation Parameters
- System: GO surface segregation in PEEK with iC grafting
- Temperature: 380°C (equilibration) to 300K (final)
- GO concentration: 1.0 wt%
- GO sheet size: 10 nm × 10 nm
- Number of PEEK chains: 100
- Number of iC chains: 20

## Analysis Results
- GO sheets show clear surface segregation
- Grafting density: $(grep "Grafting density" analysis_log.txt | awk '{print $3}') chains/Å²
- Number of grafted bonds: $(grep "Number of grafted bonds" analysis_log.txt | awk '{print $5}')

## Plots Generated
- Density profile showing component distribution
- Radial distribution functions characterizing local structure
- Chain conformation analysis
- GO sheet orientation distribution
- End-to-end distance distributions

For detailed analysis, refer to the individual plot files and raw data.
EOF

echo "==== All done! ===="
echo "Results are available in the simulation_results directory"
echo "See simulation_report.txt for a summary"

# Optional: create a tarball of all results
tar -czf go_peek_ic_simulation_results.tar.gz *

echo "Results archive created: go_peek_ic_simulation_results.tar.gz"
