from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.calculators.lj import LennardJones
from xtb.ase.calculator import XTB
from ase.neighborlist import NeighborList
from nequip.ase import nequip_calculator
from ase.build import attach



# Class for creating a multi-calculator ASE atoms object using ASE built-in calculators or custom calclulators (i.e. MLIPs!)


class MyCalculator(Calculator):
    """Calculator that combines QM/ML Interatomic Potential and Long-Range closed-form potentials for a system of multiple molecules."""
    
    def __init__(self, molecules=None, lj_params=None, lj_cutoff=None):
        """Initialize the calculator."""
        Calculator.__init__(self)
        self.molecules = molecules
        self.lj_params = lj_params
        self.lj_cutoff = lj_cutoff
    
    def calculate(self, atoms=None, properties=['energy'], system_changes='all'):
        """Calculate the energy of the system."""
        # Set up DFTB calculator for each molecule
        for i, mol in enumerate(self.molecules):
            xtb_calc = XTB(mol=mol)
            mol.set_calculator(xtb_calc)
        
        # Set up Lennard-Jones calculator for atoms between molecules
        lj_atoms = self.get_intermolecular_atoms(self.molecules, self.lj_cutoff)
        lj_calc = LennardJones(**self.lj_params)
        lj_atoms.set_calculator(lj_calc)
        
        # Calculate the energy of the system
        energy = 0.0
        for mol in self.molecules:
            energy += mol.get_potential_energy()
        energy += lj_atoms.get_potential_energy()
        
        atoms.set_calculator(self)
        atoms.set_potential_energy(energy)
    
    def get_potential_energy(self, atoms=None, force_consistent=False):
        """Get the potential energy of the system."""
        return atoms.get_potential_energy()
    
    def get_intermolecular_atoms(self, molecules, cutoff):
        """Get the atoms between different molecules."""
        nl = NeighborList([cutoff / 2] * len(molecules), self_interaction=False, bothways=True)
        positions = [mol.get_positions() for mol in molecules]
        nl.update(positions)
        lj_atoms = []
        for i, _ in enumerate(molecules[:-1]):
            for j, _ in enumerate(molecules[i+1:]):
                for offset, d in nl.get_neighbors(i, j+i+1):
                    if d < cutoff:
                        lj_atoms.append(self.atoms[offset])
        return self.atoms[lj_atoms]
    
        
