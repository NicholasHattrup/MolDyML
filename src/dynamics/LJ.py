from ase import Atoms
from ase.io import read, write
from ase.calculators.lj import LennardJones
from nequip.ase import NequIPCalculator
from ase.constraints import FixAtoms
import numpy as np
from ase.geometry import distance
from scipy.spatial.distance import pdist, squareform
from ase.calculators.lj import LennardJones
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from xtb.ase.calculator import XTB
import matplotlib.pyplot as plt
import sys


sys.path.append("/Users/nickhattrup/Documents/research/MolDyML/src/random/")
from dictionary import atoms_dict

from ase import units
from ase.md import Langevin
import torch




def get_lj_parameters(atom1, atom2):
    # Define the Lennard-Jones parameters for each atom type
    lj_parameters = {
    "H": {"epsilon": 0.00146, "sigma": 2.886},
    "C": {"epsilon": 0.00284, "sigma": 3.851},
    "O": {"epsilon": 0.0104, "sigma": 3.500},
    }


    epsilon1 = lj_parameters[atom1.symbol]['epsilon']
    sigma1 = lj_parameters[atom1.symbol]['sigma']
    epsilon2 = lj_parameters[atom2.symbol]['epsilon']
    sigma2 = lj_parameters[atom2.symbol]['sigma']

    # Combine the parameters using the Lorentz-Berthelot mixing rules
    epsilon = np.sqrt(epsilon1 * epsilon2)
    sigma = 0.5 * (sigma1 + sigma2)

    return epsilon, sigma









class CustomLennardJones(LennardJones):
    def __init__(self, molecule_indices, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.molecule_indices = molecule_indices


    def calculate(self, atoms=None, properties=('energy', 'forces'), *args, **kwargs):
        super().calculate(atoms, properties, *args, **kwargs)

        positions = self.atoms.get_positions()
        distances = squareform(pdist(positions))

        # Loop over all atom pairs and set intramolecular interaction forces and energies to zero
        energy = 0
        for i in range(len(self.atoms) - 1):
            for j in range(i + 1, len(self.atoms)):
                # If both atoms belong to the same molecule, skip the interaction
                if not any(set((i, j)).issubset(mol_indices) for mol_indices in self.molecule_indices):
                    epsilon, sigma = get_lj_parameters(self.atoms[i], self.atoms[j])
                    d = distances[i, j]
                    d6 = d**6
                    d12 = d6 * d6
                    lj_energy = 4 * epsilon * (sigma**12 / d12 - sigma**6 / d6)
                    energy += lj_energy
                    f = 48 * epsilon * (sigma**12 / d12 - 0.5 * sigma**6 / d6) / d
                    force = (positions[i] - positions[j]) * f / d
                    self.results['forces'][i] += force
                    self.results['forces'][j] -= force

        self.results['energy'] = energy





class CustomxTBLennardJones(CustomLennardJones):
    def __init__(self, molecule_indices, *args, **kwargs):
        self.XTB = XTB()
        CustomLennardJones.__init__(self, molecule_indices, *args, **kwargs)

    def calculate(self, atoms=None, properties=('energy', 'forces'), system_changes=None):
        # Calculate intramolecular energy and forces using EMT
        self.XTB.calculate(atoms, properties, system_changes)
        intramolecular_energy = self.XTB.results['energy']
        intramolecular_forces = self.XTB.results['forces']

        # Calculate intermolecular energy and forces using CustomLennardJones
        CustomLennardJones.calculate(self, atoms, properties, system_changes)  # Pass 'atoms' here
        intermolecular_energy = self.results['energy']
        intermolecular_forces = self.results['forces']

        # Combine intramolecular and intermolecular energies and forces
        self.results['energy'] = intramolecular_energy + intermolecular_energy
        self.results['forces'] = intramolecular_forces + intermolecular_forces

        # Store intramolecular and intermolecular energies separately
        self.results['intramolecular_energy'] = intramolecular_energy
        self.results['intermolecular_energy'] = intermolecular_energy



class PolyMLCalculators(CustomLennardJones):
    def __init__(self, molecule_indices, mol_types, ml_calculators, *args, **kwargs):
        self.mol_types = mol_types
        self.ml_calculators = ml_calculators
        CustomLennardJones.__init__(self, molecule_indices, *args, **kwargs)

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=None):
        # Initialize energy and forces
        ml_energy = 0.0
        ml_forces = np.zeros_like(self.atoms.get_positions())

        # Calculate intramolecular energy and forces for each unique molecule type
        for mol_type, ml_calculator in self.ml_calculators.items():
            # Get the indices of the molecules of the current type
            mol_indices_list = [mol_indices for mol_indices, mol_type_iter in zip(self.molecule_indices, self.mol_types) if mol_type_iter == mol_type]

            for mol_indices in mol_indices_list:
                # Get the atoms belonging to the current molecule
                mol_atoms = atoms[mol_indices]

                # Calculate energy and forces using the ML calculator
                ml_calculator.calculate(mol_atoms, properties, system_changes)
                ml_energy += ml_calculator.results["energy"]
                ml_forces[mol_indices] = ml_calculator.results["forces"]

        # Calculate intermolecular energy and forces using CustomLennardJones
        CustomLennardJones.calculate(self, atoms, properties, system_changes)

        # Combine intramolecular and intermolecular energies and forces
        self.results["energy"] = ml_energy + self.results["intermolecular_energy"]
        self.results["forces"] = ml_forces + self.results["intermolecular_forces"]

        # Store intramolecular and intermolecular energies separately
        self.results["intramolecular_energy"] = ml_energy
        self.results["intermolecular_energy"] = self.results["intermolecular_energy"]








'''
# Create two methanol molecules
methanol1 = read("/Users/nickhattrup/Documents/research/MolDyML/geom/methanol/methanol_one.xyz")
methanol2 = read("/Users/nickhattrup/Documents/research/MolDyML/geom/methanol/methanol_two.xyz")



# Combine the two methanol molecules
system = methanol1 + methanol2

# Define the molecule_indices for the CustomLennardJones calculator
molecule_indices = [
    [0, 1, 2, 3, 4, 5],  # Indices of the first methanol molecule
    [6, 7, 8, 9, 10, 11]   # Indices of the second methanol molecule
]
'''




'''
# Create two toluene molecules
toluene1 = read("./toluene_one.xyz")
toluene2 = read("./toluene_two.xyz")
toluene3 = read("./toluene_three.xyz")

# Combine the two toluene molecules
system = toluene1 + toluene2 + toluene3

# Define the molecule_indices for the CustomLennardJones calculator
molecule_indices = [
    [i for i in  range(15)],  # Indices of the first toluene molecule
    [i for i in  range(15, 30)],   # Indices of the second toluene molecule
    [i for i in  range(30, 45)]   # Indices of the third toluene molecule

]


'''


# Polymer system 
System = read("./../../System.xyz")
molecule_indices = [
    [i for i in  range(196)],  # Indices of the first polymer molecule
    [i for i in  range(196, 392)],   # Indices of the second polymer molecule

]




device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
Poly_calc = NequIPCalculator.from_deployed_model(
        model_path='/Users/nickhattrup/Documents/research/MolDyML/deployed_models/3_monomer/3_monomer_1.pth',
        device=device,
        species_to_type_name = atoms_dict,
        energy_units_to_eV=27.2144,
        length_units_to_A=1
    )

ml_calculators = {'Poly': Poly_calc}

mol_types = ['Poly', 'Poly']


# Set up the custom Lennard-Jones calculator
calculator = CustomxTBLennardJones(molecule_indices)

# Assign the calculator to the system
System.set_calculator(calculator)

# Define the simulation parameters
temperature = 300  # Kelvin
friction = 0.01    # Friction coefficient
time_step = .5  # Time step in fs
n_steps = 10000    # Number of steps (10 ps * 1000 fs/ps)

# Set initial velocities
MaxwellBoltzmannDistribution(System, temperature * units.kB)
# Initialize the Langevin dynamics
dyn = Langevin(System, time_step * units.fs, temperature * units.kB, friction)

# Define a function to print progress during the simulation
def print_status():
    step = dyn.get_time() // dyn.dt
    print(f"Step: {step}, Time: {dyn.get_time() * 1e-3:.2f} ps, Energy: {System.get_potential_energy():.4f} eV")

# Attach the print_status function to the dynamics
dyn.attach(print_status, interval=1000)


# Define a function to save the geometry every 10 steps
def save_geometry():
    with open("trajectory.xyz", "a") as f:
        # Write the number of atoms
        f.write(f"{len(System)}\n")
        # Write a comment line (optional, you can leave it empty)
        f.write(f"Step: {dyn.get_time() // dyn.dt}, Time: {dyn.get_time() * 1e-3:.2f} ps\n")
        # Write the atomic positions
        for atom in System:
            f.write(f"{atom.symbol} {atom.x} {atom.y} {atom.z}\n")


# Attach the save_geometry function to the dynamics
dyn.attach(save_geometry, interval=100)





# Run the simulation
dyn.run(n_steps)