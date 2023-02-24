import mbuild as mb



class Polymer:
    
    def __init__(self, xyz_file):
        self.xyz=xyz_file
        self.monomer=mb.load(self.xyz)
    
    def generate(self,index,repeats=2):
        self.polymer=mb.Polymer(self.monomer,repeats)
        self.units=repeats

    def write_xyz(self, filename=None):
        if filename is None:
            filename=self.xyz.replace('.xyz',f'{self.units}.xyz')
        mb.formats.XYZWriter(self.polymer.to_xyz()).write_file(filename)





