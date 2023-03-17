import ase
import torch
from nequip.ase import NequIPCalculator
from ase.calculators.lj import LennardJones
from src.train.initial_configuation import Initial
from src.train import md
from src.system.build import System
from src.dynamics.calculator import poly_NequIP
from src.random.dictionary import atoms_dict
import numpy as np
from ase.calculators.calculator import Calculator
import sys
# General script where I will build a python package around the workflow I hope to devise below
xyz_file = sys.argv[1]



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
    polymer_lst = [xyz_file] * 2
    box = System(xyz_files=polymer_lst)
    #box.write_xyz()
    

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    ML_calc = NequIPCalculator.from_deployed_model(
        model_path='/Users/nickhattrup/Documents/research/MolDyML/geom/toluene/model/toluene_model.pth',
        device=device,
        species_to_type_name = atoms_dict,
        energy_units_to_eV=27.2144,
        length_units_to_A=1
    )


    LR_calc = LennardJones()



    #polymer_lst = [ase.io.read(xyz_file), ase.io.read(xyz_file)]
    #calc_lst = [calculator, calculator]
    #calculator=poly_NequIP(calc1=ML_calc, calc2=LR_calc)
    
    print(box.atoms)
    print(box.mols)


