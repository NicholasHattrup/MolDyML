from ase.io import read, write
from ase import units
from ase.build.attach import attach


class System:
    def __init__(self, xyz_files):
        # Define a collection of seperate molecule objects as well as an overall Atoms object for the system 
        self.mols=[]
        if not isinstance(xyz_files,list):
            xyz_files=[xyz_files]
        mol=read(xyz_files[0])
        atoms=mol
        self.mols.append(mol)
        for i in range(1,len(xyz_files)):
            mol=read(xyz_files[i])
            atoms=attach(atoms,mol, distance=10)
            self.mols.append(mol)
        self.atoms=atoms
    
            

    def write_xyz(self, filename=None):
        if filename is None:
            filename='System.xyz'
        write(filename, self.atoms)

        





