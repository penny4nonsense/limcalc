"""
Code generation from limcalc Expr JSON to NumPy/Numba/Python expressions.

Takes the JSON representation of an Expr (as returned by the limcalc CLI)
and produces a string of Python code suitable for use with NumPy or Numba.

Example:
    expr = {"tag": "Mul", "left": {"tag": "Const", "value": 2.0},
            "right": {"tag": "Cos", "arg": {"tag": "Var", "name": "x"}}}
    to_numpy(expr)  # "2.0 * np.cos(x)"
"""

from __future__ import annotations
from typing import Any


# Precedence levels for parenthesization (mirrors Pretty.hs)
_PREC_ADD  = 0
_PREC_MUL  = 1
_PREC_NEG  = 2
_PREC_POW  = 3
_PREC_ATOM = 4


def to_numpy(expr: dict) -> str:
    """Convert a limcalc Expr JSON object to a NumPy expression string.

    The result is a valid Python expression using numpy (imported as np)
    and scipy.special for special functions.

    Args:
        expr: A dict representing a limcalc Expr in JSON form.

    Returns:
        A Python expression string.

    Example:
        >>> to_numpy({"tag": "Sin", "arg": {"tag": "Var", "name": "x"}})
        'np.sin(x)'
    """
    return _gen(expr, _PREC_ADD)


def to_python(expr: dict) -> str:
    """Convert a limcalc Expr JSON to a plain Python expression (no numpy).

    Uses math module functions instead of numpy. Useful for scalar evaluation.

    Args:
        expr: A dict representing a limcalc Expr in JSON form.

    Returns:
        A Python expression string using the math module.
    """
    return _gen(expr, _PREC_ADD, use_numpy=False)


def to_lambda(expr: dict, vars: list[str]) -> Any:
    """Convert a limcalc Expr JSON to a callable NumPy lambda.

    Args:
        expr: A dict representing a limcalc Expr in JSON form.
        vars: List of variable names (e.g. ["x"] or ["x", "y"]).

    Returns:
        A callable that accepts numpy arrays and returns a numpy array.

    Example:
        >>> import numpy as np
        >>> f = to_lambda({"tag": "Sin", "arg": {"tag": "Var", "name": "x"}}, ["x"])
        >>> f(np.array([0.0, 1.0]))
        array([0.        , 0.84147098])
    """
    import numpy as np
    import scipy.special
    code = to_numpy(expr)
    args = ", ".join(vars)
    fn = eval(f"lambda {args}: {code}",
              {"np": np, "scipy": scipy.special})
    return fn


def to_numba(expr: dict, vars: list[str]) -> Any:
    """Convert a limcalc Expr JSON to a Numba JIT-compiled function.

    Requires numba to be installed. The function is compiled on first call.

    Args:
        expr: A dict representing a limcalc Expr in JSON form.
        vars: List of variable names.

    Returns:
        A Numba JIT-compiled callable.
    """
    try:
        import numba
    except ImportError:
        raise ImportError("numba is required for to_numba(). "
                          "Install it with: pip install numba")
    fn = to_lambda(expr, vars)
    return numba.jit(fn, nopython=True)


def _paren(s: str, condition: bool) -> str:
    return f"({s})" if condition else s


def _gen(expr: dict, prec: int, use_numpy: bool = True) -> str:
    """Recursive code generation with precedence-based parenthesization."""
    tag = expr["tag"]
    np_prefix = "np" if use_numpy else "math"

    if tag == "Const":
        v = expr["value"]
        # Render integers without decimal point
        if v == int(v):
            return str(int(v))
        return repr(v)

    elif tag == "Var":
        return expr["name"]

    elif tag == "Pi":
        return f"{np_prefix}.pi"

    elif tag == "E":
        return f"{np_prefix}.e"

    elif tag == "I":
        return "1j"

    elif tag == "Add":
        left  = _gen(expr["left"],  _PREC_ADD, use_numpy)
        right = _gen(expr["right"], _PREC_ADD, use_numpy)
        s = f"{left} + {right}"
        return _paren(s, prec > _PREC_ADD)

    elif tag == "Sub":
        left  = _gen(expr["left"],  _PREC_ADD, use_numpy)
        right = _gen(expr["right"], _PREC_MUL, use_numpy)  # right-assoc
        s = f"{left} - {right}"
        return _paren(s, prec > _PREC_ADD)

    elif tag == "Mul":
        left  = _gen(expr["left"],  _PREC_MUL, use_numpy)
        right = _gen(expr["right"], _PREC_MUL, use_numpy)
        s = f"{left} * {right}"
        return _paren(s, prec > _PREC_MUL)

    elif tag == "Div":
        left  = _gen(expr["left"],  _PREC_MUL, use_numpy)
        right = _gen(expr["right"], _PREC_NEG, use_numpy)  # right-assoc
        s = f"{left} / {right}"
        return _paren(s, prec > _PREC_MUL)

    elif tag == "Pow":
        base = _gen(expr["base"], _PREC_POW, use_numpy)
        exp  = _gen(expr["exp"],  _PREC_POW, use_numpy)
        s = f"{base} ** {exp}"
        return _paren(s, prec > _PREC_POW)

    elif tag == "Neg":
        arg = _gen(expr["arg"], _PREC_NEG, use_numpy)
        return _paren(f"-{arg}", prec > _PREC_NEG)

    elif tag == "Abs":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        return f"{np_prefix}.abs({arg})"

    elif tag == "Exp":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        return f"{np_prefix}.exp({arg})"

    elif tag == "Log":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        return f"{np_prefix}.log({arg})"

    elif tag == "Sin":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        return f"{np_prefix}.sin({arg})"

    elif tag == "Cos":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        return f"{np_prefix}.cos({arg})"

    elif tag == "Arcsin":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        return f"{np_prefix}.arcsin({arg})"

    elif tag == "Arccos":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        return f"{np_prefix}.arccos({arg})"

    elif tag == "Arctan":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        return f"{np_prefix}.arctan({arg})"

    # Special functions via scipy.special
    elif tag == "Erf":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        if use_numpy:
            return f"scipy.erf({arg})"
        else:
            return f"math.erf({arg})"

    elif tag == "Si":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        if use_numpy:
            return f"scipy.sici({arg})[0]"
        else:
            return f"scipy.special.sici({arg})[0]"

    elif tag == "Ci":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        if use_numpy:
            return f"scipy.sici({arg})[1]"
        else:
            return f"scipy.special.sici({arg})[1]"

    elif tag == "Ei":
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        if use_numpy:
            return f"scipy.expi({arg})"
        else:
            return f"scipy.special.expi({arg})"

    elif tag == "Li":
        # li(x) = Ei(log(x))
        arg = _gen(expr["arg"], _PREC_ADD, use_numpy)
        if use_numpy:
            return f"scipy.expi(np.log({arg}))"
        else:
            return f"scipy.special.expi(math.log({arg}))"

    else:
        raise ValueError(f"Unknown Expr tag: {tag}")
