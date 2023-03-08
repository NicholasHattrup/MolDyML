from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.build import attach



# Class for creating a multi-calculator ASE atoms object using ASE built-in calculators or custom calclulators (i.e. MLIPs!)


class CombinedSystem:
    def __init__(self, atoms_list, calculators_list, epsilon=0.1, sigma=3.0, cutoff=10.0):
        if not isinstance(atoms_list, list):
            atoms_list = [atoms_list]
        if not isinstance(calculators_list, list):
            calculators_list = [calculators_list]

        self.atoms = atoms_list
        # Check that the number of atoms objects and calculators match
        assert len(atoms_list) == len(calculators_list)

        # Set the calculators for each atoms object
        for i, atoms in enumerate(atoms_list):
            if calculators_list[i] is not None:
                atoms.set_calculator(calculators_list[i])

        # Combine the atoms objects into a single system
        self.system = atoms_list[0]
        for atoms in atoms_list[1:]:
            self.system += atoms

        # Add non-bonded interactions 
        if len(atoms_list) > 1:
            num_mols=len(atoms_list)
            self.lj = LennardJones(epsilon=epsilon, sigma=sigma, cutoff=cutoff)
            for i in range(num_mols):
                for j in range(i+1,num_mols):
                    for atom_one in atoms_list[i]:
                        for atom_two in atoms_list[j]:
                            self.system.add_nonbonded(self.lj, atoms_one, atoms_two)
                    

    def get_energy(self):
        # Calculate the energy of the combined system
        return self.system.get_potential_energy()
    

    def build(self): # Given the atoms objects, place them in a simulation box without overlap
        system = self.atoms[0]
        for i in range(1, len(self.atoms)):
            system = attach.attach_randomly(system, self.atoms[i], distance=5)
    
        
