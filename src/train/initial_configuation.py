import autode as ade
import numpy as np


class Initial:
	def __init__(self, xyz_file, method, configs=10):
		self.xyz = xyz_file
		self.method = method
		self.num_configs = configs
		self.configs = []

	def generate(self):
		if self.method == 'mmff':
			# Load an XYZ file with the associated system to simulate 
			mol = ade.Molecule(self.xyz)
			num_atoms = len(mol.atoms)
			# Check if the molecule was loaded successfully
			if mol is None:
					print(f"Failed to load molecule from {self.xyz_file}")
			else:
				# Do something with the loaded molecule, e.g. print its SMILES string
					print('System load succesful')
			# embed multiple conformations of the polymer
			mol.populate_conformers(self.num_configs)
			self.configs=mol.conformers
	def get_xyz(self):
		for config in self.configs:
			config.print_xyz_file()
		
	

