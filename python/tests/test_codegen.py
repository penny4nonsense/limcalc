"""
Tests for limcalc.codegen — NumPy/Numba code generation from Expr JSON.

Run with: python -m pytest test_codegen.py -v
"""

import pytest
import numpy as np
from limcalc.codegen import to_numpy, to_python, to_lambda


# ---------------------------------------------------------------------------
# Atoms
# ---------------------------------------------------------------------------

def test_const_integer():
    assert to_numpy({"tag": "Const", "value": 2.0}) == "2"

def test_const_float():
    assert to_numpy({"tag": "Const", "value": 0.5}) == "0.5"

def test_const_negative():
    assert to_numpy({"tag": "Const", "value": -1.0}) == "-1"

def test_var():
    assert to_numpy({"tag": "Var", "name": "x"}) == "x"

def test_var_multichar():
    assert to_numpy({"tag": "Var", "name": "theta"}) == "theta"

def test_pi():
    assert to_numpy({"tag": "Pi"}) == "np.pi"

def test_e():
    assert to_numpy({"tag": "E"}) == "np.e"

def test_i():
    assert to_numpy({"tag": "I"}) == "1j"

# ---------------------------------------------------------------------------
# Arithmetic
# ---------------------------------------------------------------------------

def test_add():
    expr = {"tag": "Add",
            "left":  {"tag": "Var", "name": "x"},
            "right": {"tag": "Var", "name": "y"}}
    assert to_numpy(expr) == "x + y"

def test_sub():
    expr = {"tag": "Sub",
            "left":  {"tag": "Var", "name": "x"},
            "right": {"tag": "Var", "name": "y"}}
    assert to_numpy(expr) == "x - y"

def test_mul():
    expr = {"tag": "Mul",
            "left":  {"tag": "Const", "value": 2.0},
            "right": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "2 * x"

def test_div():
    expr = {"tag": "Div",
            "left":  {"tag": "Var", "name": "x"},
            "right": {"tag": "Var", "name": "y"}}
    assert to_numpy(expr) == "x / y"

def test_pow():
    expr = {"tag": "Pow",
            "base": {"tag": "Var", "name": "x"},
            "exp":  {"tag": "Const", "value": 2.0}}
    assert to_numpy(expr) == "x ** 2"

def test_neg():
    expr = {"tag": "Neg", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "-x"

# ---------------------------------------------------------------------------
# Parenthesization
# ---------------------------------------------------------------------------

def test_add_inside_mul_gets_parens():
    expr = {"tag": "Mul",
            "left":  {"tag": "Add",
                      "left":  {"tag": "Var", "name": "x"},
                      "right": {"tag": "Var", "name": "y"}},
            "right": {"tag": "Var", "name": "z"}}
    assert to_numpy(expr) == "(x + y) * z"

def test_sub_right_assoc():
    expr = {"tag": "Sub",
            "left":  {"tag": "Var", "name": "a"},
            "right": {"tag": "Sub",
                      "left":  {"tag": "Var", "name": "b"},
                      "right": {"tag": "Var", "name": "c"}}}
    assert to_numpy(expr) == "a - (b - c)"

def test_div_right_assoc():
    expr = {"tag": "Div",
            "left":  {"tag": "Var", "name": "a"},
            "right": {"tag": "Div",
                      "left":  {"tag": "Var", "name": "b"},
                      "right": {"tag": "Var", "name": "c"}}}
    assert to_numpy(expr) == "a / (b / c)"

def test_mul_no_parens():
    expr = {"tag": "Mul",
            "left":  {"tag": "Mul",
                      "left":  {"tag": "Var", "name": "x"},
                      "right": {"tag": "Var", "name": "y"}},
            "right": {"tag": "Var", "name": "z"}}
    assert to_numpy(expr) == "x * y * z"

# ---------------------------------------------------------------------------
# Elementary functions
# ---------------------------------------------------------------------------

def test_sin():
    expr = {"tag": "Sin", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "np.sin(x)"

def test_cos():
    expr = {"tag": "Cos", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "np.cos(x)"

def test_exp():
    expr = {"tag": "Exp", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "np.exp(x)"

def test_log():
    expr = {"tag": "Log", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "np.log(x)"

def test_abs():
    expr = {"tag": "Abs", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "np.abs(x)"

def test_arcsin():
    expr = {"tag": "Arcsin", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "np.arcsin(x)"

def test_arccos():
    expr = {"tag": "Arccos", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "np.arccos(x)"

def test_arctan():
    expr = {"tag": "Arctan", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "np.arctan(x)"

# ---------------------------------------------------------------------------
# Special functions
# ---------------------------------------------------------------------------

def test_erf():
    expr = {"tag": "Erf", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "scipy.erf(x)"

def test_si():
    expr = {"tag": "Si", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "scipy.sici(x)[0]"

def test_ci():
    expr = {"tag": "Ci", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "scipy.sici(x)[1]"

def test_ei():
    expr = {"tag": "Ei", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "scipy.expi(x)"

def test_li():
    expr = {"tag": "Li", "arg": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "scipy.expi(np.log(x))"

# ---------------------------------------------------------------------------
# Nested expressions
# ---------------------------------------------------------------------------

def test_sin_x_squared():
    # sin(x^2)
    expr = {"tag": "Sin",
            "arg": {"tag": "Pow",
                    "base": {"tag": "Var", "name": "x"},
                    "exp":  {"tag": "Const", "value": 2.0}}}
    assert to_numpy(expr) == "np.sin(x ** 2)"

def test_diff_sin_x_squared():
    # d/dx sin(x^2) = 2*x*cos(x^2)
    expr = {"tag": "Mul",
            "left": {"tag": "Mul",
                     "left":  {"tag": "Const", "value": 2.0},
                     "right": {"tag": "Cos",
                               "arg": {"tag": "Pow",
                                       "base": {"tag": "Var", "name": "x"},
                                       "exp":  {"tag": "Const", "value": 2.0}}}},
            "right": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "2 * np.cos(x ** 2) * x"

def test_rational_expression():
    # sin(x)/x
    expr = {"tag": "Div",
            "left":  {"tag": "Sin", "arg": {"tag": "Var", "name": "x"}},
            "right": {"tag": "Var", "name": "x"}}
    assert to_numpy(expr) == "np.sin(x) / x"

# ---------------------------------------------------------------------------
# to_python (no numpy)
# ---------------------------------------------------------------------------

def test_to_python_sin():
    expr = {"tag": "Sin", "arg": {"tag": "Var", "name": "x"}}
    assert to_python(expr) == "math.sin(x)"

def test_to_python_pi():
    assert to_python({"tag": "Pi"}) == "math.pi"

# ---------------------------------------------------------------------------
# to_lambda — numerical correctness
# ---------------------------------------------------------------------------

def test_lambda_sin():
    expr = {"tag": "Sin", "arg": {"tag": "Var", "name": "x"}}
    f = to_lambda(expr, ["x"])
    x = np.array([0.0, np.pi/2, np.pi])
    np.testing.assert_allclose(f(x), np.sin(x))

def test_lambda_cos():
    expr = {"tag": "Cos", "arg": {"tag": "Var", "name": "x"}}
    f = to_lambda(expr, ["x"])
    x = np.linspace(0, 2*np.pi, 100)
    np.testing.assert_allclose(f(x), np.cos(x))

def test_lambda_exp():
    expr = {"tag": "Exp", "arg": {"tag": "Var", "name": "x"}}
    f = to_lambda(expr, ["x"])
    x = np.array([-1.0, 0.0, 1.0, 2.0])
    np.testing.assert_allclose(f(x), np.exp(x))

def test_lambda_polynomial():
    # x^2 + 2*x + 1 = (x+1)^2
    expr = {"tag": "Add",
            "left": {"tag": "Add",
                     "left":  {"tag": "Pow",
                               "base": {"tag": "Var", "name": "x"},
                               "exp":  {"tag": "Const", "value": 2.0}},
                     "right": {"tag": "Mul",
                               "left":  {"tag": "Const", "value": 2.0},
                               "right": {"tag": "Var", "name": "x"}}},
            "right": {"tag": "Const", "value": 1.0}}
    f = to_lambda(expr, ["x"])
    x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    np.testing.assert_allclose(f(x), (x + 1)**2)

def test_lambda_multivariate():
    # x^2 + y^2
    expr = {"tag": "Add",
            "left":  {"tag": "Pow",
                      "base": {"tag": "Var", "name": "x"},
                      "exp":  {"tag": "Const", "value": 2.0}},
            "right": {"tag": "Pow",
                      "base": {"tag": "Var", "name": "y"},
                      "exp":  {"tag": "Const", "value": 2.0}}}
    f = to_lambda(expr, ["x", "y"])
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 5.0, 6.0])
    np.testing.assert_allclose(f(x, y), x**2 + y**2)

def test_lambda_erf():
    expr = {"tag": "Erf", "arg": {"tag": "Var", "name": "x"}}
    f = to_lambda(expr, ["x"])
    from scipy.special import erf
    x = np.array([0.0, 0.5, 1.0, 2.0])
    np.testing.assert_allclose(f(x), erf(x))

def test_lambda_large_array():
    # Performance check: should handle 1M points without issue
    expr = {"tag": "Sin", "arg": {"tag": "Var", "name": "x"}}
    f = to_lambda(expr, ["x"])
    x = np.linspace(0, 2*np.pi, 1_000_000)
    result = f(x)
    assert result.shape == (1_000_000,)
    np.testing.assert_allclose(result, np.sin(x))
