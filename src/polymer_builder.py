import autode as ade



monomer = ade.Molecule(smiles='CCCCCCC1=CC(C2=C[C@H]3[C@@H]4[C@H]3C(C5=C(CCCCCC)C=C(C6=C[C@H](C=C7)C=C(C8=C(CCCCCC)C=CC(CCCCCC)=C8)[C@@H]9[C@H]7[C@H]69)C(CCCCCC)=C5)=C[C@H]2C=C4)=C(CCCCCC)C=C1C%10=C[C@@H](C=C%11)C=C[C@H]%12[C@@H]%11[C@H]%12%10')

monomer.optimise(method=ade.methods.XTB())
monomer.print_xyz_file()