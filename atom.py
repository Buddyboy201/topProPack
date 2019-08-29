from mendeleev import element


element_mass = {}
class Atom:
    def __init__(self, symbol, name, atomid, coords):
        self.symbol = symbol.capitalize()
        self.name = name
        self.atoid = atomid
        self.coords = coords
        if element_mass.get(self.symbol) ==  None:
            element_mass[self.symbol] = element(self.symbol).atomic_weight
        self.atomic_mass = element_mass[self.symbol]
        
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

    


    
        

