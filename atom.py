from mendeleev import element
import re

element_mass = {} # why is this needed?
class Atom:
    def __init__(self, symbol, name, atomid, coords):
        self.symbol = symbol.capitalize()
        self.name = name
        self.atomid = atomid
        self.coords = coords # coords a list? or can be a numpy array?
        if element_mass.get(self.symbol) ==  None:
            element_mass[self.symbol] = element(self.symbol).atomic_weight
            # why would element mass be zero, but atomic_weight be fine?
        self.atomic_mass = element_mass[self.symbol]
        m = re.search(r'CA$|C$|N$|O$|', self.name)
        if m.group(0) != "":
            self.mcsc = 'mc'
        else:
            self.mcsc = 'sc' 
        
    def get_mcsc(self):
        return self.mcsc
    
    def get_symbol(self):
        return self.symbol

    def get_atomid(self):
        return self.atomid

    def get_name(self):
        return self.name

    def get_coords(self):
        return self.coords

    def get_mass(self):
        return self.atomic_mass