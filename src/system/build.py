import ase
from ase import units
from ase.calculators.emt import EMT
from ase.calculators.mmff import MMFF94
from ase.data import bond_lengths
from ase.optimize import BFGS



class Polymer:
    
    def __init__(self, xyz_file):
        self.xyz=xyz_file
        self.monomer = ase.io.read(xyz_file.xyz)
        self.atoms=len(self.monomer.get_atoms())
    
    def generate(self,index,repeats=2):
         self.polymer=self.monomer
         self.units=repeats
         for n in range(1,repeats):
            atom1 = self.polymer.get_atoms()[index+(n-1)*self.atoms]
            atom2 = self.monomer.get_atoms()[index]
            bond_length = bond_lengths[(atom1, atom2)]
            self.polymer[self.polymer.get_number_of_atoms()-1].symbol = 'C'
            self.polymer[self.polymer.get_number_of_atoms()-1].position = atom2.position + bond_length * np.array([1, 0, 0])
            

    def write_xyz(self, filename=None):
        if filename is None:
            filename=self.xyz.replace('.xyz',f'{self.units}.xyz')
        mb.formats.XYZWriter(self.polymer.to_xyz()).write_file(filename)





