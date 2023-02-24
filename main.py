from src.train.initial_configuation import Initial


# General script where I will build a python package around the workflow I hope to devise below
xyz_file = '/Users/nhattrup/Documents/research/Fluxional_MD/nequip/molecules/golder_polymer/3_units/3_monomer.xyz'

configs=Initial(xyz_file=xyz_file,method='mmff',configs=50)
configs.generate()
print(configs.configs)