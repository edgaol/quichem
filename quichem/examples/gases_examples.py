import os
import periodictable as pt
import sys
from pint import Quantity as Q_

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

from core import Molecule
from gases import law, u

def ex_4_1():
    v0 = 2.00 * u.L
    p0 = 2.00 * u.atm
    p = 3.00 * u.atm
    v = law(p0=p0, v0=v0, p=p)
    print(f"4.1) The new volume is {v}.")

def ex_4_2():
    a = Q_(45, u.degC)
    b = Q_(-74, u.degC)
    c = Q_(273, u.degC)
    print(f"4.2) The temperatures are a={a.to('kelvin')}, b={b.to('kelvin')} and c={c.to('kelvin')}.")

def ex_4_3():
    v0 = 2 * u.L
    t0 = 273.15 * u.K
    t = 373.15 * u.K
    v = law(v0=v0, t0=t0, t=t)
    print(f"4.3) The new volume is {v}.")

def ex_4_4():
    v0 = 1 * u.L
    p0 = 2 * u.atm
    t0 = 283.15 * u.K
    t = 373.15 * u.K
    p = 1 * u.atm
    v = law(v0=v0, p0=p0, t0=t0, t=t, p=p)
    print(f"4.4) The new volume is {v}.")

def ex_4_5():
    n = 2
    p = 1.5
    t = 400.15
    v = law(n=n, p=p, t=t)
    print(f"4.5) The new volume is {v}.")

def ex_4_6():
    co2 = Molecule("CO2", 22 * u.g)
    p = 0.9
    t = 283.15
    v = law(p=p, t=t, n=co2.mol)
    print(f"4.6) The volume is {v}.")

def ex_4_7():
    v0 = 1
    p0 = 1
    t0 = 283.15
    t = 373.15
    p = law(p0=p0, t0=t0, t=t)
    print(f"4.7) The new pressure is {p}.")

def ex_4_8():
    # Oxygen gas volume approximately 21% of air's
    t0 = 298.15
    p0 = 0.93
    v = law(t0=t0, p0=p0)
    vo2 = 0.21 * v
    print(f"4.8) The volume will be of {vo2}.")

def ex_4_9():
    # Reaction: 4 TNG(l) -> 6 N2(g) + 12 CO2(g) + 10 H2O(g) + 1 O2(g)
    # Since the TNG's volume is negligible, the initial mol count (n) for the air occupying the recipient needs to be calculated.
    # Given stoichometry, there will produced a total of 2.9 mol of gas by the 0.4 inital TNG. There is also 0.4 mol of air, so 3.3 mol total gas.
    p0 = 1
    v0 = 10
    t0 = 298.15
    t = 573.15
    # _n0 = law(p=1, v=v0, t=t0)
    n = 3.3
    p = law(n=n, t=t, v=v0)
    print(f"4.9) The new pressure would be of {p}.")

def ex_4_10():
    # Reaction: KClO3 -> KCl + 3/2 O2
    kclo3 = Molecule("KClO3", 3 * u.g)
    o2 = Molecule("O2", kclo3.mol * 3/2)
    p = 0.93
    t = 298.15
    v = law(n=o2.mol, t=t, p=p)
    print(f"4.10) The volume is {v}.")

if __name__ == "__main__":
    ex_4_1()
    ex_4_2()
    ex_4_3()
    ex_4_4()
    ex_4_5()
    ex_4_6()
    ex_4_7()
    ex_4_8()
    ex_4_9()
    ex_4_10()
