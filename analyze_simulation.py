#!/usr/bin/env python3
"""
Analysis script for GO/PEEK/iC molecular dynamics simulation results
This script processes LAMMPS output files to analyze:
1. Density profiles
2. Radial distribution functions
3. Chain conformations
4. GO sheet orientation
5. iC grafting characteristics
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import os

class SimulationAnalysis:
    """Class for analyzing LAMMPS simulation results"""
    
    def __init__(self, output_dir='.'):
        """Initialize with output directory containing LAMMPS files"""
        self.output_dir = output_dir
    
    def load_density_profile(self, filename='density_profile.dat'):
        """Load and parse density profile data"""
        filepath = os.path.join(self.output_dir, filename)
        data = np.loadtxt(filepath, skiprows=4)  # Skip LAMMPS header
        
        # Extract columns: z-position and density values for each atom type
        z_pos = data[:, 0]
        densities = data[:, 1:]
        
        return z_pos, densities
    
    def load_rdf_data(self, filename='rdf.gograft'):
        """Load and parse radial distribution function data"""
        filepath = os.path.join(self.output_dir, filename)
        data = np.loadtxt(filepath, skiprows=4)  # Skip LAMMPS header
        
        # Extract columns: r values and g(r) for different pairs
        r_values = data[:, 1]
        rdf_go_go = data[:, 2]  # GO-GO
        rdf_go_peek = data[:, 3]  # GO-PEEK
        rdf_go_ic = data[:, 4]  # GO-iC
        rdf_peek_peek = data[:, 5]  # PEEK-PEEK
        rdf_ic_ic = data[:, 6]  # iC-iC
        
        return r_values, {
            'GO-GO': rdf_go_go,
            'GO-PEEK': rdf_go_peek,
            'GO-iC': rdf_go_ic,
            'PEEK-PEEK': rdf_peek_peek,
            'iC-iC': rdf_ic_ic
        }
    
    def load_gyration_data(self, filename='gyration.dat'):
        """Load and parse radius of gyration data"""
        filepath = os.path.join(self.output_dir, filename)
        data = np.loadtxt(filepath, skiprows=4)  # Skip LAMMPS header
        
        # Extract columns: timestep, gyration_peek, gyration_ic
        timesteps = data[:, 0]
        gyration_peek = data[:, 1]
        gyration_ic = data[:, 2]
        
        return timesteps, gyration_peek, gyration_ic
    
    def load_dump_file(self, filename, frame=-1):
        """Load atom coordinates from LAMMPS dump file for a specific frame"""
        filepath = os.path.join(self.output_dir, filename)
        
        # Read the entire file to find frame boundaries
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Find the line indices where each frame starts
        frame_starts = [i for i, line in enumerate(lines) if 'ITEM: TIMESTEP' in line]
        
        if frame < 0:
            # Use the last frame if negative index
            frame_idx = frame_starts[frame]
        else:
            # Use the specified frame
            if frame >= len(frame_starts):
                raise ValueError(f"Frame {frame} not found in dump file")
            frame_idx = frame_starts[frame]
        
        # Read number of atoms in this frame
        n_atoms_line = lines[frame_idx + 3]
        n_atoms = int(n_atoms_line.strip())
        
        # Read box bounds
        box_lines = lines[frame_idx + 5:frame_idx + 8]
        box_bounds = [list(map(float, line.split())) for line in box_lines]
        
        # Find where atom data starts
        atoms_start = frame_idx + 9
        
        # Parse atom data
        atom_data = []
        for i in range(atoms_start, atoms_start + n_atoms):
            if i >= len(lines):
                break
            atom_info = lines[i].split()
            atom_id = int(atom_info[0])
            atom_type = int(atom_info[1])
            x, y, z = map(float, atom_info[2:5])
            
            # Additional columns may include velocities if available
            vx, vy, vz = 0.0, 0.0, 0.0
            if len(atom_info) > 5:
                vx, vy, vz = map(float, atom_info[5:8])
            
            atom_data.append({
                'id': atom_id,
                'type': atom_type,
                'x': x, 'y': y, 'z': z,
                'vx': vx, 'vy': vy, 'vz': vz
            })
        
        return atom_data, box_bounds
    
    def analyze_go_orientation(self, dump_filename, frame=-1):
        """Analyze the orientation of GO sheets with respect to the surface"""
        atom_data, box_bounds = self.load_dump_file(dump_filename, frame)
        
        # Group atoms by molecule ID (to identify individual GO sheets)
        molecules = {}
        for atom in atom_data:
            if atom['type'] in [3, 4]:  # GO atom types
                molecule_id = atom['id'] // 1000  # Approximate molecule ID
                if molecule_id not in molecules:
                    molecules[molecule_id] = []
                molecules[molecule_id].append(atom)
        
        orientations = []
        
        # For each GO sheet, calculate normal vector using principal component analysis
        for molecule_id, atoms in molecules.items():
            if len(atoms) < 10:  # Skip if too few atoms
                continue
            
            # Extract coordinates
            coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])
            
            # Center the coordinates
            centered_coords = coords - np.mean(coords, axis=0)
            
            # Calculate covariance matrix
            cov_matrix = np.cov(centered_coords.T)
            
            # Calculate eigenvalues and eigenvectors
            eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
            
            # The normal vector is the eigenvector corresponding to the smallest eigenvalue
            normal_vector = eigenvectors[:, 0]
            
            # Calculate angle with z-axis (surface normal)
            z_axis = np.array([0, 0, 1])
            angle = np.arccos(np.abs(np.dot(normal_vector, z_axis)))
            angle_degrees = np.degrees(angle)
            
            orientations.append(angle_degrees)
        
        return orientations
    
    def analyze_grafting_density(self, dump_filename, frame=-1):
        """Analyze the density of iC grafting on the GO/PEEK surface"""
        atom_data, box_bounds = self.load_dump_file(dump_filename, frame)
        
        # Identify grafted iC molecules
        grafted_bonds = self.identify_grafted_bonds(dump_filename, frame)
        
        # Calculate grafting density
        total_surface_area = box_bounds[0][1] * box_bounds[1][1]  # xy area
        grafting_density = len(grafted_bonds) / total_surface_area
        
        return grafting_density, grafted_bonds
    
    def identify_grafted_bonds(self, dump_filename, frame=-1):
        """Identify bonds between iC chains and the GO/PEEK surface"""
        # This requires bond information which is not directly in dump files
        # For a real analysis, we would need to either:
        # 1. Use a custom LAMMPS command to output bond information
        # 2. Infer bonds based on distance criteria
        
        # For this template, we'll use a simple distance-based approach
        atom_data, box_bounds = self.load_dump_file(dump_filename, frame)
        
        # Group atoms by type
        go_atoms = [atom for atom in atom_data if atom['type'] in [3, 4]]
        peek_atoms = [atom for atom in atom_data if atom['type'] in [1, 2]]
        ic_atoms = [atom for atom in atom_data if atom['type'] in [5, 6]]
        
        # Find the z-coordinate of the top of the GO/PEEK surface
        surface_atoms = go_atoms + peek_atoms
        max_z_surface = max(atom['z'] for atom in surface_atoms)
        
        # Identify iC atoms that are close to the surface
        bond_distance_threshold = 2.0  # Å
        grafted_bonds = []
        
        for ic_atom in ic_atoms:
            # Only consider the first bead of each iC chain
            if ic_atom['type'] != 5:
                continue
            
            for surface_atom in surface_atoms:
                dx = ic_atom['x'] - surface_atom['x']
                dy = ic_atom['y'] - surface_atom['y']
                dz = ic_atom['z'] - surface_atom['z']
                
                distance = np.sqrt(dx**2 + dy**2 + dz**2)
                
                if distance < bond_distance_threshold:
                    grafted_bonds.append((surface_atom['id'], ic_atom['id']))
                    break
        
        return grafted_bonds
    
    def calculate_chain_conformation(self, dump_filename, frame=-1):
        """Calculate conformational parameters of polymer chains"""
        atom_data, box_bounds = self.load_dump_file(dump_filename, frame)
        
        # Group atoms by molecule ID and type
        peek_molecules = {}
        ic_molecules = {}
        
        for atom in atom_data:
            molecule_id = atom['id'] // 1000  # Approximate molecule ID
            
            if atom['type'] in [1, 2]:  # PEEK
                if molecule_id not in peek_molecules:
                    peek_molecules[molecule_id] = []
                peek_molecules[molecule_id].append(atom)
            
            elif atom['type'] in [5, 6]:  # iC
                if molecule_id not in ic_molecules:
                    ic_molecules[molecule_id] = []
                ic_molecules[molecule_id].append(atom)
        
        # Calculate end-to-end distances
        peek_e2e = []
        ic_e2e = []
        
        for molecule_id, atoms in peek_molecules.items():
            if len(atoms) < 3:
                continue
            
            # Sort atoms by ID to get the chain sequence
            atoms.sort(key=lambda x: x['id'])
            
            # Calculate end-to-end distance
            dx = atoms[-1]['x'] - atoms[0]['x']
            dy = atoms[-1]['y'] - atoms[0]['y']
            dz = atoms[-1]['z'] - atoms[0]['z']
            
            distance = np.sqrt(dx**2 + dy**2 + dz**2)
            peek_e2e.append(distance)
        
        for molecule_id, atoms in ic_molecules.items():
            if len(atoms) < 2:
                continue
            
            # Sort atoms by ID to get the chain sequence
            atoms.sort(key=lambda x: x['id'])
            
            # Calculate end-to-end distance
            dx = atoms[-1]['x'] - atoms[0]['x']
            dy = atoms[-1]['y'] - atoms[0]['y']
            dz = atoms[-1]['z'] - atoms[0]['z']
            
            distance = np.sqrt(dx**2 + dy**2 + dz**2)
            ic_e2e.append(distance)
        
        return {
            'PEEK_e2e': peek_e2e,
            'iC_e2e': ic_e2e
        }
    
    def plot_density_profile(self, save_file='density_profile.png'):
        """Plot density profiles along z direction"""
        z_pos, densities = self.load_density_profile()
        
        plt.figure(figsize=(10, 6))
        plt.plot(z_pos, densities[:, 0], label='PEEK')
        plt.plot(z_pos, densities[:, 1], label='GO')
        plt.plot(z_pos, densities[:, 2], label='iC')
        
        plt.xlabel('z position (Å)')
        plt.ylabel('Density (g/cm³)')
        plt.title('Density Profile')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.savefig(save_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_rdf(self, save_file='rdf_plot.png'):
        """Plot radial distribution functions"""
        r_values, rdf_data = self.load_rdf_data()
        
        plt.figure(figsize=(10, 6))
        for pair, g_r in rdf_data.items():
            plt.plot(r_values, g_r, label=pair)
        
        plt.xlabel('r (Å)')
        plt.ylabel('g(r)')
        plt.title('Radial Distribution Functions')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.savefig(save_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_gyration(self, save_file='gyration_plot.png'):
        """Plot radius of gyration over time"""
        timesteps, gyration_peek, gyration_ic = self.load_gyration_data()
        
        # Convert timesteps to nanoseconds
        time_ns = timesteps * 0.00001  # Assuming timestep = 10 fs
        
        plt.figure(figsize=(10, 6))
        plt.plot(time_ns, gyration_peek, label='PEEK')
        plt.plot(time_ns, gyration_ic, label='iC')
        
        plt.xlabel('Time (ns)')
        plt.ylabel('Radius of Gyration (Å)')
        plt.title('Chain Conformation over Time')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.savefig(save_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_go_orientation(self, dump_filename, save_file='go_orientation.png', frame=-1):
        """Plot distribution of GO sheet orientations"""
        orientations = self.analyze_go_orientation(dump_filename, frame)
        
        plt.figure(figsize=(8, 6))
        
        # Plot histogram
        plt.hist(orientations, bins=18, alpha=0.6, density=True)
        
        # Add kernel density estimate
        if len(orientations) > 2:
            density = gaussian_kde(orientations)
            x = np.linspace(0, 90, 100)
            plt.plot(x, density(x), 'r-', lw=2)
        
        plt.xlabel('Angle with Surface Normal (degrees)')
        plt.ylabel('Probability Density')
        plt.title('GO Sheet Orientation Distribution')
        plt.grid(True, alpha=0.3)
        
        plt.savefig(save_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_e2e_distribution(self, dump_filename, save_file='e2e_distribution.png', frame=-1):
        """Plot end-to-end distance distributions"""
        conformation_data = self.calculate_chain_conformation(dump_filename, frame)
        
        peek_e2e = conformation_data['PEEK_e2e']
        ic_e2e = conformation_data['iC_e2e']
        
        plt.figure(figsize=(10, 6))
        
        # PEEK distribution
        if len(peek_e2e) > 2:
            plt.hist(peek_e2e, bins=20, alpha=0.5, label='PEEK', density=True)
            density = gaussian_kde(peek_e2e)
            x = np.linspace(min(peek_e2e), max(peek_e2e), 100)
            plt.plot(x, density(x), 'r-', lw=2)
        
        # iC distribution
        if len(ic_e2e) > 2:
            plt.hist(ic_e2e, bins=20, alpha=0.5, label='iC', density=True)
            density = gaussian_kde(ic_e2e)
            x = np.linspace(min(ic_e2e), max(ic_e2e), 100)
            plt.plot(x, density(x), 'b-', lw=2)
        
        plt.xlabel('End-to-End Distance (Å)')
        plt.ylabel('Probability Density')
        plt.title('Chain End-to-End Distance Distribution')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.savefig(save_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    def run_complete_analysis(self, dump_filename='dump.gograft.final'):
        """Run all analysis functions and generate plots"""
        print("Starting comprehensive analysis...")
        
        print("Plotting density profile...")
        self.plot_density_profile()
        
        print("Plotting radial distribution functions...")
        self.plot_rdf()
        
        print("Plotting gyration data...")
        self.plot_gyration()
        
        print("Analyzing GO orientation...")
        self.plot_go_orientation(dump_filename)
        
        print("Analyzing chain conformations...")
        self.plot_e2e_distribution(dump_filename)
        
        print("Calculating grafting density...")
        grafting_density, grafted_bonds = self.analyze_grafting_density(dump_filename)
        print(f"Grafting density: {grafting_density:.6f} chains/Å²")
        print(f"Number of grafted bonds: {len(grafted_bonds)}")
        
        print("Analysis complete.")

if __name__ == "__main__":
    analysis = SimulationAnalysis()
    analysis.run_complete_analysis()