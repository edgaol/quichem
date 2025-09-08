from __future__ import annotations
from pint import Quantity, UnitRegistry
from sympy import Matrix, lcm
from sys import maxsize
from typing import Dict, List
import periodictable as pt
import re

u = UnitRegistry()
Q_ = u.Quantity

def f_(formula: str) -> str:
    subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

    def replacer(match):
        element = match.group(1)
        number = match.group(2)
        if number:
            return element + number.translate(subscript_map)
        else:
            return element

    return re.sub(r"([A-Z][a-z]*)(\d*)", replacer, formula)

def formula_parser(formula: str) -> Dict[str, int]:
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula)
    elements = {}
    for e, n in matches:
        elements[e] = 1 if n == "" else int(n)
    return elements

def molecular_weight(elements: Dict[str, int]) -> Quantity:
    mw = 0
    for e, n in elements.items():
        mw += n * getattr(pt, e).mass
    mw *= u.g / u.mol
    return mw


class Molecule():
    def __init__(self, formula: str, quantity: Quantity=None, precision=3):
        self.formula = formula
        self.elements = formula_parser(formula)
        self.mw = molecular_weight(self.elements)
        if quantity:
            if quantity.units == u.moles:
                self.mol = quantity
            elif quantity.units == u.gram:
                self.mol = quantity / self.mw
            else:
                raise TypeError(f"Unit {quantity.unit} not implemented for Molecule's input")
        self.precision = precision

    def format_quantity(self, q: Quantity):
        return f"{q:.{self.precision}e}"
    
    @property
    def grams(self):
        return self.mol * self.mw

    @property
    def atoms(self):
        return (self.mol * pt.constants.avogadro_number).magnitude
    
    def describe(self):
        def fill(header: str, q: Quantity):
            return f"| {header}" + (20-len(header))*" " + ":    " + self.format_quantity(q) + (12-len(str(q.units)))*" " + "|\n"
        
        out = ""
        out += str(self) + " " + (27+len(str(self.mw)))*"-" + "\n"
        out += fill("Molecular weight", self.mw)
        out += fill("Moles", self.mol)
        out += fill("Mass", self.grams)
        out += (30+len(str(self.mw)))*"-" + "\n"
        return out

    @property
    def centesimal(self) -> Dict[str, float]:
        centesimals = {}
        for (e, count) in self.elements.items():
            val = count * getattr(pt, e).mass * u.g / u.mole / self.mw
            val = round(val, self.precision)
            centesimals[e] = val.magnitude

        return centesimals

    def __str__(self):
        return f_(self.formula)

    def __add__(self, mol: Molecule):
        return [self, mol]


class Reaction():
    def __init__(self, reactants: List[str|Molecule]=[], products: List[str|Molecule]=[]):
        self.r = {}
        self.p = {}
        for r in reactants:
            if not isinstance(r, Molecule):
                r = Molecule(r)
            self.r[r.formula] = r
        for p in products:
            if not isinstance(p, Molecule):
                p = Molecule(p)
            self.p[p.formula] = p
        self.stoichiometry = {}
        self.limitant = None

    def balance(self, inplace=False):
        elements_table = []
        elements = []
        elements_blacklist = set()
        for mol in self.r.values():
            for e in mol.elements:
                if e not in elements_blacklist:
                    elements.append(e)
                    elements_blacklist.add(e)
        for mol in self.p.values():
            for e in mol.elements:
                if e not in elements_blacklist:
                    elements.append(e)
                    elements_blacklist(e)
        for e in elements:
            row = []
            rl = len(self.r)
            for i, mol in enumerate(list(self.r.values()) + list(self.p.values())):
                n = mol.elements.get(e)
                if n:
                    row.append(n if i < rl else -n)
                else:
                    row.append(0)
            elements_table.append(row)
        a = Matrix(elements_table)
        ns = a.nullspace()
        if ns:
            coeffs = ns[0]
            multiplier = lcm([term.q for term in coeffs])
            integers = [int(term * multiplier) for term in coeffs]
            formulas = list(self.r.keys()) + list(self.p.keys())
            self.stoichiometry = dict(zip(formulas, integers))
        else:
            raise ValueError("The reaction could not be balanced.")

    def find_limitant(self):
        if not self.stoichiometry:
            self.balance()
        moles = maxsize * u.mol
        for r in self.r.values():
            if r.mol < moles:
                moles = r.mol / self.stoichiometry[r.formula]
                self.limitant = r

    def __str__(self):
        out = ""
        if self.stoichiometry:
            rl = len(self.r)
            pl = len(self.p)
            species = list(self.stoichiometry.keys())
            coefficients = list(self.stoichiometry.values())
            for i in range(rl):
                out += f"{coefficients[i]} {f_(species[i])} + "
            out = out[:-2] + "-> "
            for i in range(rl, rl+pl):
                out += f"{coefficients[i]} {f_(species[i])} + "
        else:
            for r in self.r.keys():
                out += f"{f_(r)} + "
            out = out[:-2] + "-> "
            for p in self.p.keys():
                out += f"{f_(p)} + "
        return out[:-3]
    
    def __getitem__(self, index: str):
        mol = self.r.get(index)
        if not mol:
            mol = self.p.get(index)
        return mol

    @property
    def reactants(self):
        return list(map(str, self.r))

    @property
    def products(self):
        return list(map(str, self.p))
    
