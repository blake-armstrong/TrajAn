import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import distances

# Load the PDB file
u = mda.Universe("test2.pdb")

# Get all atoms
atoms = u.atoms

# Get coordinates of all atoms
coords = atoms.positions

# Calculate pairwise distances between all atoms
dist_array = distances.distance_array(coords, coords)

# Count pairs within 6 Angstroms (excluding self-pairs)
# Only count upper triangle to avoid counting each pair twice
mask = np.logical_and(dist_array <= 6.0, np.triu(np.ones_like(dist_array), k=1) > 0)
pair_count = np.sum(mask)

print(f"Number of atom pairs within 6 Angstroms: {pair_count}")
