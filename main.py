import ase
import torch
from nequip.ase import NequIPCalculator
from src.train.initial_configuation import Initial
from src.train import md
from src.random import random
from src.random.dictionary import atoms_dict
from src.dynamics.calculator import CombinedSystem
# General script where I will build a python package around the workflow I hope to devise below
xyz_file = '/home/nhattrup/MolDyML/geom/3_monomer/3_monomer_conf0_siman.xyz'



'''
configs=Initial(xyz_file=xyz_file,method='mmff',configs=12)
configs.generate()
print(configs.configs)
'''




#md.concatenate_xyz_files(directories='/Users/nickhattrup/Documents/research/MolDyML/geometries/training_data/dynamics/MD/3_monomer',substring='forces')
#md.combine_pos_force_xyz(pos_file='/Users/nickhattrup/Documents/research/MolDyML/geometries/training_data/dynamics/MD/3_monomer/positions.xyz',
#                         force_file='/Users/nickhattrup/Documents/research/MolDyML/geometries/training_data/dynamics/MD/3_monomer/forces.xyz')


'''
_, unique_indexes=md.get_unique_frames('/home/nhattrup/MolDyML/data/3_monomer/positions.xyz')

random.extract_configurations('data/3_monomer/combined.xyz', unique_indexes)
'''

if __name__ == '__main__':


    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    calculator = NequIPCalculator.from_deployed_model(
        model_path='/home/nhattrup/MolDyML/deployed_models/3_monomer/3_monomer_1.pth',
        device=device,
        species_to_type_name = atoms_dict,
        energy_units_to_eV=27.2144,
        length_units_to_A=1,
    )




    polymer_lst = [ase.io.read(xyz_file), ase.io.read(xyz_file)]
    calc_lst = [calculator, calculator] 
    system=CombinedSystem(atoms_list=polymer_lst, calculators_list=calc_lst) # Initialize a combined system class
    print(system.get_energy())


