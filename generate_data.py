#!/usr/bin/env python3
"""
Data file generator for GO/PEEK/iC coarse-grained molecular dynamics simulation
This script creates a LAMMPS data file with proper molecular structures for
PEEK polymer chains, graphene oxide sheets, and iota-carrageenan chains.
"""

import numpy as np
import random

# System parameters
box_x = 200.0  # A
box_y = 200.0  # A
box_z = 300.0  # A

# PEEK parameters
n_peek_chains = 100
n_peek_units_per_chain = 50
peek_bond_length = 1.53  # A
peek_angle = 120.0  # degrees

# GO parameters
go_concentrations = [0.1, 0.5, 1.0, 2.0]  # wt%
go_conc = go_concentrations[2]  # Select 1.0 wt% for this simulation
go_sheet_sizes = [(50, 50), (100, 100)]  # A
go_sheet_size = go_sheet_sizes[1]  # Select 10 nm × 10 nm
go_bond_length = 1.40  # A

# iC parameters
n_ic_chains = 20
n_ic_units_per_chain = 10
ic_bond_length = 1.45  # A
ic_angle = 110.0  # degrees

class Atom:
    """Class representing an atom in the simulation"""
    def __init__(self, atom_id, molecule_id, atom_type, charge, x, y, z):
        self.atom_id = atom_id
        self.molecule_id = molecule_id
        self.atom_type = atom_type
        self.charge = charge
        self.x = x
        self.y = y
        self.z = z
    
    def __str__(self):
        return f"{self.atom_id} {self.molecule_id} {self.atom_type} {self.charge:.1f} {self.x:.3f} {self.y:.3f} {self.z:.3f}"

class Bond:
    """Class representing a bond in the simulation"""
    def __init__(self, bond_id, bond_type, atom_id1, atom_id2):
        self.bond_id = bond_id
        self.bond_type = bond_type
        self.atom_id1 = atom_id1
        self.atom_id2 = atom_id2
    
    def __str__(self):
        return f"{self.bond_id} {self.bond_type} {self.atom_id1} {self.atom_id2}"

class Angle:
    """Class representing an angle in the simulation"""
    def __init__(self, angle_id, angle_type, atom_id1, atom_id2, atom_id3):
        self.angle_id = angle_id
        self.angle_type = angle_type
        self.atom_id1 = atom_id1
        self.atom_id2 = atom_id2
        self.atom_id3 = atom_id3
    
    def __str__(self):
        return f"{self.angle_id} {self.angle_type} {self.atom_id1} {self.atom_id2} {self.atom_id3}"

def generate_peek_chain(chain_id, start_atom_id, start_x, start_y, start_z):
    """Generate a PEEK chain with specified parameters"""
    atoms = []
    bonds = []
    angles = []
    
    atom_id = start_atom_id
    prev_atom_id = None
    prev_prev_atom_id = None
    bond_id = start_atom_id
    angle_id = start_atom_id
    
    x, y, z = start_x, start_y, start_z
    
    for unit in range(n_peek_units_per_chain):
        # Each PEEK unit has 3 beads: two ether-linked (type 1) and one ketone-linked (type 2)
        for bead in range(3):
            atom_type = 2 if bead == 1 else 1  # Middle bead is ketone-linked
            atoms.append(Atom(atom_id, chain_id, atom_type, 0.0, x, y, z))
            
            # Create bond to previous atom
            if prev_atom_id is not None:
                bonds.append(Bond(bond_id, 1, prev_atom_id, atom_id))
                bond_id += 1
            
            # Create angle with two previous atoms
            if prev_prev_atom_id is not None:
                angles.append(Angle(angle_id, 1, prev_prev_atom_id, prev_atom_id, atom_id))
                angle_id += 1
            
            prev_prev_atom_id = prev_atom_id
            prev_atom_id = atom_id
            atom_id += 1
            
            # Move to next position (simple linear chain for now)
            theta = np.radians(120.0)  # 120-degree angle
            if bead == 0:
                x += peek_bond_length * np.cos(theta)
                y += peek_bond_length * np.sin(theta)
            elif bead == 1:
                x += peek_bond_length
            else:
                x += peek_bond_length * np.cos(theta)
                y -= peek_bond_length * np.sin(theta)
        
        # Add some random displacement for realistic amorphous structure
        x += random.uniform(-0.5, 0.5)
        y += random.uniform(-0.5, 0.5)
        z += random.uniform(-0.5, 0.5)
    
    return atoms, bonds, angles

def generate_go_sheet(sheet_id, start_atom_id, start_x, start_y, start_z, width, height):
    """Generate a GO sheet with specified parameters"""
    atoms = []
    bonds = []
    angles = []
    
    atom_id = start_atom_id
    bond_id = start_atom_id
    angle_id = start_atom_id
    
    # Calculate number of beads in each dimension
    nx = int(width / go_bond_length)
    ny = int(height / go_bond_length)
    
    # Create 2D array of atom IDs for bond and angle generation
    atom_id_grid = np.zeros((nx, ny), dtype=int)
    
    # Generate atoms in a 2D grid
    for i in range(nx):
        for j in range(ny):
            x = start_x + i * go_bond_length
            y = start_y + j * go_bond_length
            z = start_z
            
            # Randomly assign oxygen-containing groups (type 4) with 10% probability
            # Rest are carbon beads (type 3)
            # Edges have higher probability (40%) of oxygen groups
            is_edge = (i == 0 or i == nx-1 or j == 0 or j == ny-1)
            prob_oxygen = 0.4 if is_edge else 0.1
            atom_type = 4 if random.random() < prob_oxygen else 3
            
            atoms.append(Atom(atom_id, sheet_id, atom_type, 0.0, x, y, z))
            atom_id_grid[i, j] = atom_id
            atom_id += 1
    
    # Generate bonds between adjacent atoms
    for i in range(nx):
        for j in range(ny):
            current_id = atom_id_grid[i, j]
            
            # Bond to right neighbor
            if i < nx - 1:
                right_id = atom_id_grid[i+1, j]
                bonds.append(Bond(bond_id, 2, current_id, right_id))
                bond_id += 1
            
            # Bond to bottom neighbor
            if j < ny - 1:
                bottom_id = atom_id_grid[i, j+1]
                bonds.append(Bond(bond_id, 2, current_id, bottom_id))
                bond_id += 1
    
    # Generate angles (simplified - only considering some key angles)
    for i in range(1, nx-1):
        for j in range(ny):
            left_id = atom_id_grid[i-1, j]
            current_id = atom_id_grid[i, j]
            right_id = atom_id_grid[i+1, j]
            angles.append(Angle(angle_id, 2, left_id, current_id, right_id))
            angle_id += 1
    
    for i in range(nx):
        for j in range(1, ny-1):
            top_id = atom_id_grid[i, j-1]
            current_id = atom_id_grid[i, j]
            bottom_id = atom_id_grid[i, j+1]
            angles.append(Angle(angle_id, 2, top_id, current_id, bottom_id))
            angle_id += 1
    
    return atoms, bonds, angles

def generate_ic_chain(chain_id, start_atom_id, start_x, start_y, start_z):
    """Generate an iota-carrageenan chain with specified parameters"""
    atoms = []
    bonds = []
    angles = []
    
    atom_id = start_atom_id
    prev_atom_id = None
    prev_prev_atom_id = None
    bond_id = start_atom_id
    angle_id = start_atom_id
    
    x, y, z = start_x, start_y, start_z
    
    for unit in range(n_ic_units_per_chain):
        # Each iC unit has 2 beads: b-D-galactopyranose (type 5) and a-D-galactopyranose with sulfate (type 6)
        for bead in range(2):
            atom_type = 5 if bead == 0 else 6
            charge = -0.5 if bead == 0 else -1.0  # Sulfate groups have negative charge
            
            atoms.append(Atom(atom_id, chain_id, atom_type, charge, x, y, z))
            
            # Create bond to previous atom
            if prev_atom_id is not None:
                bonds.append(Bond(bond_id, 3, prev_atom_id, atom_id))
                bond_id += 1
            
            # Create angle with two previous atoms
            if prev_prev_atom_id is not None:
                angles.append(Angle(angle_id, 3, prev_prev_atom_id, prev_atom_id, atom_id))
                angle_id += 1
            
            prev_prev_atom_id = prev_atom_id
            prev_atom_id = atom_id
            atom_id += 1
            
            # Move to next position with helical structure
            theta = np.radians(ic_angle)
            phi = np.radians(60.0 * unit + 180.0 * bead)  # Create helical structure
            x += ic_bond_length * np.sin(theta) * np.cos(phi)
            y += ic_bond_length * np.sin(theta) * np.sin(phi)
            z += ic_bond_length * np.cos(theta)
    
    return atoms, bonds, angles

def calculate_go_sheets_count(peek_mass, go_wt_percent, sheet_width, sheet_height):
    """Calculate number of GO sheets needed for specified weight percentage"""
    # Approximate mass of a GO sheet
    nx = int(sheet_width / go_bond_length)
    ny = int(sheet_height / go_bond_length)
    atoms_per_sheet = nx * ny
    
    # Assuming average atom mass of 65 Da (mix of carbon and oxygen)
    sheet_mass = atoms_per_sheet * 65.0
    
    # Total GO mass based on weight percentage
    go_mass = peek_mass * go_wt_percent / (100 - go_wt_percent)
    
    # Number of sheets needed
    n_sheets = max(1, int(go_mass / sheet_mass))
    return n_sheets

def generate_lammps_data_file():
    """Generate complete LAMMPS data file"""
    all_atoms = []
    all_bonds = []
    all_angles = []
    
    # Initialize counters
    atom_id_counter = 1
    molecule_id_counter = 1
    
    # Generate PEEK chains
    for i in range(n_peek_chains):
        start_x = random.uniform(10, box_x - 10)
        start_y = random.uniform(10, box_y - 10)
        start_z = random.uniform(10, box_z - 10)
        
        atoms, bonds, angles = generate_peek_chain(molecule_id_counter, atom_id_counter, start_x, start_y, start_z)
        all_atoms.extend(atoms)
        all_bonds.extend(bonds)
        all_angles.extend(angles)
        
        atom_id_counter += n_peek_units_per_chain * 3
        molecule_id_counter += 1
    
    # Calculate approximate mass of PEEK
    peek_mass = n_peek_chains * n_peek_units_per_chain * 3 * 80.0  # 80 Da average per bead
    
    # Calculate number of GO sheets
    n_go_sheets = calculate_go_sheets_count(peek_mass, go_conc, go_sheet_size[0], go_sheet_size[1])
    
    # Generate GO sheets
    for i in range(n_go_sheets):
        start_x = random.uniform(10, box_x - go_sheet_size[0] - 10)
        start_y = random.uniform(10, box_y - go_sheet_size[1] - 10)
        start_z = random.uniform(10, box_z - 10)
        
        atoms, bonds, angles = generate_go_sheet(molecule_id_counter, atom_id_counter, 
                                                start_x, start_y, start_z,
                                                go_sheet_size[0], go_sheet_size[1])
        all_atoms.extend(atoms)
        all_bonds.extend(bonds)
        all_angles.extend(angles)
        
        # Update counters based on actual number of atoms generated
        atom_id_counter += len(atoms)
        molecule_id_counter += 1
    
    # Generate iC chains above the system (for later grafting)
    for i in range(n_ic_chains):
        start_x = random.uniform(10, box_x - 10)
        start_y = random.uniform(10, box_y - 10)
        start_z = box_z - 50.0  # Start near the top of the box
        
        atoms, bonds, angles = generate_ic_chain(molecule_id_counter, atom_id_counter, start_x, start_y, start_z)
        all_atoms.extend(atoms)
        all_bonds.extend(bonds)
        all_angles.extend(angles)

        atom_id_counter += n_ic_units_per_chain * 2
        molecule_id_counter += 1
    
    # Create reactive sites on GO sheets (to be activated during simulation)
    # This is just a placeholder - actual sites will be determined during simulation
    
    # Write the data file
    with open('data.gograft', 'w', encoding='utf-8') as f:
        f.write("# LAMMPS data file for GO/PEEK/iC system\n")
        f.write("# Generated by Python script\n\n")
        
        f.write(f"{len(all_atoms)} atoms\n")
        f.write(f"{len(all_bonds)} bonds\n")
        f.write(f"{len(all_angles)} angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n\n")
        
        f.write("7 atom types\n")
        f.write("4 bond types\n")
        f.write("3 angle types\n")
        f.write("0 dihedral types\n")
        f.write("0 improper types\n\n")
        
        f.write(f"0.0 {box_x} xlo xhi\n")
        f.write(f"0.0 {box_y} ylo yhi\n")
        f.write(f"0.0 {box_z} zlo zhi\n\n")
        
        f.write("Masses\n\n")
        f.write("1 76.09   # Ether-linked phenyl bead (PEEK)\n")
        f.write("2 90.08   # Ketone-linked phenyl bead (PEEK)\n")
        f.write("3 60.05   # GO carbon bead\n")
        f.write("4 76.05   # GO oxygen-containing group bead\n")
        f.write("5 162.14  # 3-linked b-D-galactopyranose bead (iC)\n")
        f.write("6 242.20  # 4-linked a-D-galactopyranose with sulfate bead (iC)\n")
        f.write("7 76.05   # Reactive site bead\n\n")
        
        f.write("Atoms  # full atom style: atom-ID molecule-ID atom-type q x y z\n\n")
        for atom in all_atoms:
            f.write(f"{atom}\n")
        
        f.write("\nBonds\n\n")
        for bond in all_bonds:
            f.write(f"{bond}\n")
        
        f.write("\nAngles\n\n")
        for angle in all_angles:
            f.write(f"{angle}\n")
    
    print(f"Generated LAMMPS data file with {len(all_atoms)} atoms, {len(all_bonds)} bonds, and {len(all_angles)} angles")

if __name__ == "__main__":
    generate_lammps_data_file()