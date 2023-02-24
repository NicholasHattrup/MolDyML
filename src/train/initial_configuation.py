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
			num_atoms = mol.GetNumAtoms()
			positions = np.empty(shape=(num_atoms,3)) # Generate numpy array to contain atomic positions 
			# Check if the molecule was loaded successfully
			if mol is None:
					print(f"Failed to load molecule from {self.xyz_file}")
			else:
				# Do something with the loaded molecule, e.g. print its SMILES string
					print('System load succesful')
			# embed multiple conformations of the polymer
			conformer_generator = ade.conformers.ConformerGenerator()
			conformer_generator.generate_conformers(mol, self.num_configs,mmff=True)
			self.configs=mol.conformers
		
				
				
		
	

