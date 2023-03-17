import autode as ade

mol = ade.Molecule(smiles='CC1=CC=CC=C1')
mol.optimise(method=ade.methods.XTB())
