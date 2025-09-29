from core import u
from numbers import Number
from pint import Quantity, UnitRegistry
from typing import Any
import inspect
import re

def ideal_gases(
    p0: Quantity | Number=None,
    v0: Quantity | Number=None,
    t0: Quantity | Number=None,
    p: Quantity | Number=None,
    v: Quantity | Number=None,
    t: Quantity | Number=None,
    n: Quantity | Number=1.0 * u.mole,
    r: Quantity | Number=0.08206 * u.atm*u.L/(u.mole*u.K)
) -> (Quantity, str):
    sig = inspect.signature(ideal_gases)
    params = sig.parameters.keys()
    params = list(params)
    data = {}
    i = 0
    non_empty_params = []
    for x in params:
        val = locals()[x]
        if i < 6:
            i += 1
            if val:
                non_empty_params.append(x)
        if val and not isinstance(val, Quantity):
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
            data[x] = val
    e_ = _number_exists
    np = len(non_empty_params)
    if np == 2:
        use_0, use_1 = False, False
        try:
            use_0 = all(x[1] == "0" for x in non_empty_params)
        except Exception as e:
            print(e)
        try:
            use_1 = all(len(x) == 1 for x in non_empty_params)
        except Exception as e:
            print(e)
        if use_0:
            if use_1:
                raise ValueError("To calculate PV = nRT, provide 2 among \"p0, v0 and t0\" or \"p, v and t\".")
            else:
                return _f(data["n"], data["r"], data["p0"], data["v0"], data["t0"])
        else:
            if use_1:
                return _f(data["n"], data["r"], data["p"], data["v"], data["t"])
            else:
                raise ValueError("Something went wrong.")
    elif np == 3:
        return _g(data["p0"], data["v0"], data["t0"], data["p"], data["v"], data["t"])
    elif np == 5:
        return _h(data["p0"], data["v0"], data["t0"], data["p"], data["v"], data["t"])
    else:
        raise ValueError("Unmatched amount of variables.")

def _number_exists(x: Any) -> Any | int:
    return x if x else 1

def _f(
    n: Quantity, r: Quantity, p: Quantity=None, v: Quantity=None, 
    t: Quantity=None
) -> Quantity:
    """
    PVT⁻¹ = nR
    """
    if not p:
        return n * r * t / v
    elif not v:
        return n * r * t / p
    elif not t:
        return p * v / (n * r)
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
    p0: Quantity=None, v0: Quantity=None, t0: Quantity=None, 
    p: Quantity=None, v: Quantity=None, t: Quantity=None    
) -> Quantity:
    """
    PoVoTo⁻¹ = PVT⁻¹
    """
    if not p0:
        return p * v * t0 / (t * v0)
    elif not p:
        return p0 * v0 * t / (t0 * v)
    elif not v0:
        return p * v * t0 / (t * p0)
    elif not v:
        return p0 * v0 * t / (t0 * p)
    elif not t0:
        return p0 * v0 * t / (p * v)
    elif not t:
        return p * v * t0 / (p0 * v0)
    else:
        raise ValueError("Something went wrong.")
    