import autode as ade



mol = ade.Molecule(smiles='C1=CC2C3C2C=CC1C=C3')
mol.optimise(method=ade.methods.XTB())

