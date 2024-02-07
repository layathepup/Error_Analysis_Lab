import string
from molmass import Formula


class Solution:
    solute_id: string
    solute_mass: float
    volume: float
    solute_molmass: float

    class Concentration:
        mass: float
        volume: float

        def __init__(self, solution):
            self.mass = solution.solute_mass / solution.volume
            self.mol = self.mass / solution.solute_molmass

    concentration: Concentration

    def __init__(self, solute_id, solute_mass, volume):
        # need generic setattr to avoid calling overriden setattr before initialization
        object.__setattr__(self, 'solute_id', solute_id)
        object.__setattr__(self, 'solute_mass', solute_mass)
        object.__setattr__(self, 'volume', volume)
        object.__setattr__(self, 'solute_molmass', Formula(self.solute_id).mass)
        object.__setattr__(self, 'solute_mass', solute_mass)
        object.__setattr__(self, 'concentration', self.Concentration(self))

    def solute_id(self):
        return self.solute_id

    def __setattr__(self, key, value):
        super().__setattr__(key, value)
        self.concentration = self.Concentration(self)
        return value
