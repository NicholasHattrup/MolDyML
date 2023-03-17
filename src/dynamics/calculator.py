from ase import Atoms
from ase.calculators.calculator import Calculator
from nequip.ase import NequIPCalculator
from xtb.ase.calculator import XTB
from ase.neighborlist import NeighborList
from nequip.ase import nequip_calculator
from ase.build import attach
from ase.calculators.calculator import Calculator, all_properties
from ase.calculators.calculator import PropertyNotImplementedError
import numpy as np



# Class for creating a multi-calculator ASE atoms object using ASE built-in calculators or custom calclulators (i.e. MLIPs!)


class poly_NequIP(Calculator):
    """Composite calculator that combines the results of two calculators."""

    def __init__(self, calc1, calc2):
        Calculator.__init__(self)
        self.calc1 = calc1 #Machine learning or QM calculator 
        self.calc2 = calc2 #Long range calculator 

    def calculate(self, molecules=None, mask=None, properties=all_properties):
        for molecule in molecules:
            self.calc1.calculate(molecule, properties)
        if mask is None:
            mask = np.True_(shape=())
        
        

        '''
        self.calc1.calculate(atoms, properties)
        self.calc2.calculate(atoms, properties)
        self.results = {}
        for prop in properties:
            if prop in self.calc1.results and prop in self.calc2.results:
                self.results[prop] = self.calc1.results[prop] + self.calc2.results[prop]
            elif prop in self.calc1.results:
                self.results[prop] = self.calc1.results[prop]
            elif prop in self.calc2.results:
                self.results[prop] = self.calc2.results[prop]
            else:
                raise PropertyNotImplementedError(prop)


        '''



class PolyCalculator(Calculator):
    def __init__(self, ml_model, epsilon, sigma, n_atoms_molecule_1, **kwargs):
        super().__init__(**kwargs)
        self.ml_model = ml_model
        self.epsilon = epsilon
        self.sigma = sigma
        self.n_atoms_molecule_1 = n_atoms_molecule_1

    def lj_force(self, r, epsilon, sigma):
        q = (sigma / r) ** 6
        return -48 * epsilon * q * (q - 0.5) / r

    def calculate(self, atoms=None, properties=None, system_changes=None):
        super().calculate(atoms, properties, system_changes)

        positions = atoms.get_positions()
        n_atoms = len(positions)

        # Calculate intramolecular potential using the ML model
        intramolecular_energy = 0
        forces = np.zeros((n_atoms, 3))

        for molecule_atoms in [positions[:self.n_atoms_molecule_1], positions[self.n_atoms_molecule_1:]]:
            intramolecular_energy += self.ml_model.predict_energy(molecule_atoms)
            forces += self.ml_model.predict_forces(molecule_atoms)

        # Calculate intermolecular potential using the LJ potential
        intermolecular_energy = 0
        for i in range(self.n_atoms_molecule_1):
            for j in range(self.n_atoms_molecule_1, n_atoms):
                r_vec = positions[i] - positions[j]
                r = np.linalg.norm(r_vec)
                intermolecular_energy += lj_potential(r, self.epsilon, self.sigma)

                f_ij = self.lj_force(r, self.epsilon, self.sigma) * r_vec / r
                forces[i] += f_ij
                forces[j] -= f_ij

        total_energy = intramolecular_energy + intermolecular_energy

        self.results["energy"] = total_energy
        self.results["forces"] = forces
    
        
