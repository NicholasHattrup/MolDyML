import ase 
from src.train.initial_configuation import Initial
from src.train import md
from src.system.build import System
from src.dynamics.calculator import CombinedSystem
# General script where I will build a python package around the workflow I hope to devise below
xyz_file = '/Users/nickhattrup/Documents/research/MolDyML/geom/3_monomer/3_monomer_conf0_siman.xyz'



'''
configs=Initial(xyz_file=xyz_file,method='mmff',configs=12)
configs.generate()
print(configs.configs)
'''




#md.concatenate_xyz_files(directories='/Users/nickhattrup/Documents/research/MolDyML/geometries/training_data/dynamics/MD/3_monomer',substring='forces')
#md.combine_pos_force_xyz(pos_file='/Users/nickhattrup/Documents/research/MolDyML/geometries/training_data/dynamics/MD/3_monomer/positions.xyz',
#                         force_file='/Users/nickhattrup/Documents/research/MolDyML/geometries/training_data/dynamics/MD/3_monomer/forces.xyz')

#md.get_unique_frames('/Users/nickhattrup/Documents/research/MolDyML/geometries/training_data/dynamics/MD/3_monomer/positions.xyz')






if __name__ == '__main__':
    polymer_lst = [ase.io.read(xyz_file)] * 2 
    box = System(xyz_files=polymer_lst)
    box.write_xyz()
    #calc_lst = 
    #system=CombinedSystem() # Initialize a combined system class




