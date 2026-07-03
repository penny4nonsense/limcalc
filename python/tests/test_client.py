"""
Tests for limcalc.client — subprocess wrapper around limcalc-cli.

Requires the limcalc-cli binary to be available, either via the
LIMCALC_CLI environment variable or on the system PATH. Tests are
skipped automatically if no binary is found.

Run with: uv run pytest tests/ -v
"""

import os
import pytest
from limcalc.client import LimCalc, LimCalcError, _find_cli


# ---------------------------------------------------------------------------
# Skip all tests if CLI binary is not available
# ---------------------------------------------------------------------------

def _cli_available() -> bool:
    try:
        _find_cli()
        return True
    except FileNotFoundError:
        return False


skip_if_no_cli = pytest.mark.skipif(
    not _cli_available(),
    reason="limcalc-cli binary not found. Set LIMCALC_CLI env var."
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def lc():
    """A single LimCalc client shared across all tests in this module."""
    client = LimCalc()
    yield client
    client.close()


# ---------------------------------------------------------------------------
# Helper expressions
# ---------------------------------------------------------------------------

def var(name):
    return {"tag": "Var", "name": name}

def const(v):
    return {"tag": "Const", "value": float(v)}

def sin(arg):
    return {"tag": "Sin", "arg": arg}

def cos(arg):
    return {"tag": "Cos", "arg": arg}

def exp(arg):
    return {"tag": "Exp", "arg": arg}

def log(arg):
    return {"tag": "Log", "arg": arg}

def add(a, b):
    return {"tag": "Add", "left": a, "right": b}

def mul(a, b):
    return {"tag": "Mul", "left": a, "right": b}

def div(a, b):
    return {"tag": "Div", "left": a, "right": b}

def pow_(base, exp_):
    return {"tag": "Pow", "base": base, "exp": exp_}

def neg(arg):
    return {"tag": "Neg", "arg": arg}


# ---------------------------------------------------------------------------
# diff tests
# ---------------------------------------------------------------------------

@skip_if_no_cli
def test_diff_sin(lc):
    result = lc.diff(sin(var("x")), "x")
    assert result["tag"] == "Cos"
    assert result["arg"]["tag"] == "Var"
    assert result["arg"]["name"] == "x"

@skip_if_no_cli
def test_diff_cos(lc):
    result = lc.diff(cos(var("x")), "x")
    assert result["tag"] == "Neg"
    assert result["arg"]["tag"] == "Sin"

@skip_if_no_cli
def test_diff_exp(lc):
    result = lc.diff(exp(var("x")), "x")
    assert result["tag"] == "Exp"

@skip_if_no_cli
def test_diff_polynomial(lc):
    # d/dx x^2 = 2x
    result = lc.diff(pow_(var("x"), const(2)), "x")
    assert result["tag"] == "Mul"

@skip_if_no_cli
def test_diff_constant_is_zero(lc):
    result = lc.diff(const(5), "x")
    assert result["tag"] == "Const"
    assert result["value"] == 0.0

@skip_if_no_cli
def test_diff_var_is_one(lc):
    result = lc.diff(var("x"), "x")
    assert result["tag"] == "Const"
    assert result["value"] == 1.0


# ---------------------------------------------------------------------------
# integrate tests
# ---------------------------------------------------------------------------

@skip_if_no_cli
def test_integrate_const(lc):
    # int 1 dx = x
    result = lc.integrate(const(1), "x")
    assert result["tag"] == "Var"
    assert result["name"] == "x"

@skip_if_no_cli
def test_integrate_exp(lc):
    # int exp(x) dx = exp(x)
    result = lc.integrate(exp(var("x")), "x")
    assert result["tag"] == "Exp"

@skip_if_no_cli
def test_integrate_sin(lc):
    # int sin(x) dx = -cos(x)
    result = lc.integrate(sin(var("x")), "x")
    assert result["tag"] == "Neg"
    assert result["arg"]["tag"] == "Cos"

@skip_if_no_cli
def test_integrate_nonelementary_raises(lc):
    # int exp(x^2) dx is non-elementary
    with pytest.raises(LimCalcError, match="NonElementary"):
        lc.integrate(exp(pow_(var("x"), const(2))), "x")

@skip_if_no_cli
def test_integrate_erf(lc):
    # int exp(-x^2) dx = (sqrt(pi)/2) * erf(x)
    expr = exp(neg(pow_(var("x"), const(2))))
    result = lc.integrate(expr, "x")
    assert result["tag"] == "Mul"


# ---------------------------------------------------------------------------
# limit tests
# ---------------------------------------------------------------------------

@skip_if_no_cli
def test_limit_sinc(lc):
    # lim_{x->0} sin(x)/x = 1
    expr = div(sin(var("x")), var("x"))
    result = lc.limit(expr, "x", 0.0)
    assert abs(result - 1.0) < 1e-6

@skip_if_no_cli
def test_limit_polynomial(lc):
    # lim_{x->2} x^2 = 4
    expr = pow_(var("x"), const(2))
    result = lc.limit(expr, "x", 2.0)
    assert abs(result - 4.0) < 1e-6

@skip_if_no_cli
def test_limit_exp(lc):
    # lim_{x->0} exp(x) = 1
    result = lc.limit(exp(var("x")), "x", 0.0)
    assert abs(result - 1.0) < 1e-6

@skip_if_no_cli
def test_limit_pole_raises(lc):
    # lim_{x->0} 1/x is a pole
    expr = div(const(1), var("x"))
    with pytest.raises(LimCalcError, match="Pole"):
        lc.limit(expr, "x", 0.0)


# ---------------------------------------------------------------------------
# simplify tests
# ---------------------------------------------------------------------------

@skip_if_no_cli
def test_simplify_zero_plus_x(lc):
    expr = add(const(0), var("x"))
    result = lc.simplify(expr)
    assert result["tag"] == "Var"
    assert result["name"] == "x"

@skip_if_no_cli
def test_simplify_x_times_one(lc):
    expr = mul(var("x"), const(1))
    result = lc.simplify(expr)
    assert result["tag"] == "Var"

@skip_if_no_cli
def test_simplify_double_neg(lc):
    expr = neg(neg(var("x")))
    result = lc.simplify(expr)
    assert result["tag"] == "Var"


# ---------------------------------------------------------------------------
# pretty tests
# ---------------------------------------------------------------------------

@skip_if_no_cli
def test_pretty_sin(lc):
    result = lc.pretty(sin(var("x")))
    assert result == "sin(x)"

@skip_if_no_cli
def test_pretty_cos(lc):
    result = lc.pretty(cos(var("x")))
    assert result == "cos(x)"

@skip_if_no_cli
def test_pretty_polynomial(lc):
    # x^2
    result = lc.pretty(pow_(var("x"), const(2)))
    assert result == "x^2"


# ---------------------------------------------------------------------------
# gradient tests
# ---------------------------------------------------------------------------

@skip_if_no_cli
def test_gradient_x_squared_plus_y_squared(lc):
    # grad(x^2 + y^2) = [2x, 2y]
    expr = add(pow_(var("x"), const(2)), pow_(var("y"), const(2)))
    result = lc.gradient(expr, ["x", "y"])
    assert len(result) == 2
    # Both components should be Mul nodes (2*x and 2*y)
    assert result[0]["tag"] == "Mul"
    assert result[1]["tag"] == "Mul"

@skip_if_no_cli
def test_gradient_constant_is_zero(lc):
    result = lc.gradient(const(5), ["x", "y"])
    assert all(r["tag"] == "Const" and r["value"] == 0.0 for r in result)


# ---------------------------------------------------------------------------
# Context manager
# ---------------------------------------------------------------------------

@skip_if_no_cli
def test_context_manager():
    with LimCalc() as client:
        result = client.diff(sin(var("x")), "x")
        assert result["tag"] == "Cos"


# ---------------------------------------------------------------------------
# Persistent subprocess — verify reuse across calls
# ---------------------------------------------------------------------------

@skip_if_no_cli
def test_subprocess_reused(lc):
    pid1 = lc._proc.pid if lc._proc else None
    lc.diff(sin(var("x")), "x")
    lc.diff(cos(var("x")), "x")
    pid2 = lc._proc.pid if lc._proc else None
    assert pid1 == pid2
