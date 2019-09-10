import atom

class Residue:
    def __init__(self, name, resid, atoms):
        self.name = name # Residue type - prefer single letter code
        self.resid = resid # ID number for residue in structure
        self.atoms = atoms # what is atoms? List of atom positions?
        self.centroid = None # or COM?

    def get_atoms(self):
        return self.atoms

    def get_resid(self):
        return self.resid

    def get_COM(self):
        if self.name == "GLY" or self.name == "G": # Make flexible with single letter or triple letter codes. 
            if self.atoms[1].get_name() == "CA": return self.atoms[1].get_coords()
            else: print("No CA atom found for GLY molecule at {}".format(self.resid))          
        COM = [0.0, 0.0, 0.0]
        mass_sum = 0
        for atm in self.atoms:
            COM[0] += atm.get_mass()*atm.get_coords()[0]
            COM[1] += atm.get_mass()*atm.get_coords()[1]
            COM[2] += atm.get_mass()*atm.get_coords()[2]
            mass_sum += atm.get_mass()
        COM[0] /= float(mass_sum)
        COM[1] /= float(mass_sum)
        COM[2] /= float(mass_sum)
        self.centroid = tuple(COM) # Set centroid attribute. 
        return tuple(COM)

    def get_centroid(self):
        if self.centroid == None:
            self.centroid = self.get_COM
        return self.centroid

    def add_atom(self, atm):
        self.atoms.append(atm)

    def update_COM(self):
        self.centroid = self.get_COM()