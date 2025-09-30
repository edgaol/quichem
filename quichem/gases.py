from core import u
from numbers import Number
from pint import Quantity
from sympy import solve
from sympy.abc import x as X
from typing import List
import re

def law(
    p0: Quantity | Number=None,
    v0: Quantity | Number=None,
    t0: Quantity | Number=None,
    n0: Quantity | Number=1.0 * u.mole,
    a0: Quantity | Number=0 * u.L**2*u.atm/u.mole**2,
    b0: Quantity | Number=0 * u.L/u.mole,
    p: Quantity | Number=None,
    v: Quantity | Number=None,
    t: Quantity | Number=None,
    n: Quantity | Number=1.0 * u.mole,
    a: Quantity | Number=0 * u.L**2*u.atm/u.mole**2,
    b: Quantity | Number=0 * u.L/u.mole,
    r: Quantity | Number=0.08206 * u.atm*u.L/(u.mole*u.K)
) -> (Quantity, str):
    data = {}
    i = 0
    non_empty_params = []
    for x in ("p0", "v0", "t0", "p", "v", "t", "a0", "b0", "n0", "a", "b", "n", "r"):
        val = locals()[x]
        if i < 6:
            i += 1
            if val:
                non_empty_params.append(x)
        if not isinstance(val, Quantity):
            if val:
                if re.match(x[0], r"p"):
                    data[x] = val * u.atm
                elif re.match(x[0], r"v"):
                    data[x] = val * u.L
                elif re.match(x[0], r"t"):
                    data[x] = val * u.K
                elif re.match(x, r"n"):
                    data[x] = val * u.mole
                elif re.match(x, r"r"):
                    data[x] = val * u.atm*u.L/(u.mole*u.K)
            else:
                if re.match(x[0], r"a"):
                    data[x] = val * u.L**2*u.atm/u.mole**2
                elif re.match(x[0], r"b"):
                    data[x] = val * u.L/u.mole
                else:
                    data[x] = val
        else:
            data[x] = val
    np = len(non_empty_params)
    if np == 2:
        use_0, use_1 = False, False
        try:
            use_0 = all(x[1] == "0" for x in non_empty_params)
        except:
            pass
        try:
            use_1 = all(len(x) == 1 for x in non_empty_params)
        except:
            pass
        if use_0:
            if use_1:
                raise ValueError(
                    "To calculate PV = nRT, provide 2 among \"p0, v0 and t0\" or \"p, v and t\".")
            else:
                return _f(
                    data["a0"], data["b0"], data["n0"], data["r"],
                    data["p0"], data["v0"], data["t0"])
        else:
            if use_1:
                return _f(
                    data["a"], data["b"], data["n"], data["r"],
                    data["p"], data["v"], data["t"])
            else:
                raise ValueError("Something went wrong.")
    elif np == 3:
        return _g( 
            data["p0"], data["v0"], data["t0"], data["p"], data["v"],
            data["t"])
    elif np == 5:
        return _h(
            data["a0"], data["b0"], data["n0"], data["a"], data["b"], 
            data["n"], data["p0"], data["v0"], data["t0"], data["p"], 
            data["v"], data["t"])
    else:
        raise ValueError("Unmatched amount of variables.")

def _f(
    a: Quantity, b: Quantity, n: Quantity, r: Quantity, 
    p: Quantity=None, v: Quantity=None, t: Quantity=None
) -> Quantity:
    """
    PVT⁻¹ = nR
    """
    if not p:
        return n * r * t / (v-n*b) - n**2*a/v**2
    elif not v:
        sol = solve(-n.m**3*a.m*b.m*X**(-2) + n.m**2*a.m*X**(-1) + p.m*X - p.m*n.m*b.m -n.m*r.m*t.m)
        units = n.units * r.units * t.units / p.units
        return Quantity(_parse_solution(sol), units)
    elif not t:
        return (p+n**2*a/v**2) * (v-n*b) / (n * r)
    else:
        raise ValueError("Something went wrong.")
    
def _g(
    p0: Quantity=None, v0: Quantity=None, t0: Quantity=None,
    p: Quantity=None, v: Quantity=None, t: Quantity=None
) -> Quantity:
    """
    VoTo⁻¹ = VT⁻¹ or
    PoTo⁻¹ = PT⁻¹ or
    PoVo = PV
    """
    if not (p0 or p):
        if not v0:
            return v * t0 / t
        elif not v:
            return v0 * t / t0
        elif not t0:
            return v0 * t / v
        elif not t:
            return v * t0 / v0
    elif not (v0 or v):
        if not p0:
            return p * t0 / t
        elif not p:
            return p0 * t / t0
        elif not t0:
            return p0 * t / p
        elif not t:
            return p * t0 / p0
    elif not(t or t0):
        if not p0:
            return p * v / v0
        elif not p:
            return p0 * v0 / v
        elif not v0:
            return p * v / p0
        elif not v:
            return p0 * v0 / p
    else:
        raise ValueError("Something went wrong.")

def _h(
    a0: Quantity, b0: Quantity, n0: Quantity, a: Quantity, b: Quantity,
    n: Quantity, p0: Quantity=None, v0: Quantity=None, t0: Quantity=None,
    p: Quantity=None, v: Quantity=None, t: Quantity=None
) -> Quantity:
    """
    PoVoTo⁻¹ = PVT⁻¹
    """
    if not (t0 and t):
        if not p0:
            return (p-n**2*a/v**2) * (v-n**b) / ((v0-n0*b0)) + n0**2*a0/v0**2
        elif not p:
            return (p0-n0**a0/v0**2) * (v0-n0**b0) / ((v-n*b)) + n**2*a/v**2
        elif not v0:
            return p * v * t0 / (t * p0)
        elif not v:
            return p0 * v0 * t / (t0 * p)
    if not p0:
        return (p-n**2*a/v**2) * (v-n**b) * t0 / ((v0-n0*b0)*t) + n0**2*a0/v0**2
    elif not p:
        return (p0-n0**a0/v0**2) * (v0-n0**b0) * t / ((v-n*b)*t0) + n**2*a/v**2
    elif not v0:
        return p * v * t0 / (t * p0)
    elif not v:
        return p0 * v0 * t / (t0 * p)
    elif not t0:
        return (p0-n0**2*a0/v0**2) * (v0-n0*b0) * t / ((p-n**2*a/v**2) * (v-n*b))
    elif not t:
        return (p-n**2*a/v**2) * (v-n*b) * t / ((p0-n0**2*a0/v0**2) * (v0-n0*b0))
    else:
        raise ValueError("Something went wrong.")

def _parse_solution(sol: List[float]) -> float:
    return [s for s in sol if s.is_real][0]
