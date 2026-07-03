"""
limcalc — Limit-based symbolic calculus engine.

A Python interface to the limcalc Haskell engine, providing symbolic
differentiation, integration, and limit computation via log-Puiseux
series expansion. Results can be converted to NumPy or Numba code
for fast numerical evaluation.

Basic usage:

    import limcalc as lc

    # Differentiate sin(x^2)
    f = lc.diff("sin(x^2)", "x")
    print(f.pretty())      # "2*x*cos(x^2)"
    print(f.to_numpy())    # "2 * x * np.cos(x ** 2)"

    # Integrate exp(-x^2)
    g = lc.integrate("exp(-x^2)", "x")
    print(g.pretty())      # "(sqrt(π)/2)*erf(x)"

    # Compute a limit
    print(lc.limit("sin(x)/x", "x", 0))  # 1.0

    # Gradient of a multivariate expression
    grad = lc.gradient("x^2 + y^2", ["x", "y"])

Note: String expression input requires the LaTeX parser (coming soon).
Until then, use Expr JSON dicts directly via LimCalc() client.
"""

from __future__ import annotations

from limcalc.client import LimCalc, LimCalcError
from limcalc.codegen import to_numpy, to_python, to_lambda, to_numba

__version__ = "0.1.0"
__all__ = [
    "Expr",
    "diff",
    "partial_diff",
    "integrate",
    "limit",
    "simplify",
    "gradient",
    "LimCalc",
    "LimCalcError",
    "to_numpy",
    "to_python",
    "to_lambda",
    "to_numba",
]

# Module-level client — created lazily on first use
_client: LimCalc | None = None


def _get_client() -> LimCalc:
    """Return the module-level LimCalc client, creating it if needed."""
    global _client
    if _client is None or (_client._proc and _client._proc.poll() is not None):
        _client = LimCalc()
    return _client


class Expr:
    """A symbolic expression result from limcalc.

    Wraps an Expr JSON dict and provides convenient methods for
    pretty-printing and code generation.

    Attributes:
        _expr: The raw Expr JSON dict.

    Example:
        f = diff_expr({"tag": "Sin", "arg": {"tag": "Var", "name": "x"}}, "x")
        print(f.pretty())    # "cos(x)"
        print(f.to_numpy())  # "np.cos(x)"
    """

    def __init__(self, expr: dict):
        self._expr = expr

    def pretty(self) -> str:
        """Return a human-readable string representation.

        Example:
            Expr({"tag": "Cos", "arg": {"tag": "Var", "name": "x"}}).pretty()
            # "cos(x)"
        """
        return _get_client().pretty(self._expr)

    def to_numpy(self) -> str:
        """Return a NumPy expression string.

        Example:
            Expr(...).to_numpy()  # "np.cos(x)"
        """
        return to_numpy(self._expr)

    def to_python(self) -> str:
        """Return a plain Python expression string (using math module).

        Example:
            Expr(...).to_python()  # "math.cos(x)"
        """
        return to_python(self._expr)

    def to_lambda(self, vars: list[str]):
        """Return a callable NumPy function.

        Args:
            vars: List of variable names in the expression.

        Returns:
            A callable that accepts numpy arrays.

        Example:
            import numpy as np
            f = Expr(...).to_lambda(["x"])
            f(np.linspace(0, 1, 100))
        """
        return to_lambda(self._expr, vars)

    def to_numba(self, vars: list[str]):
        """Return a Numba JIT-compiled function.

        Args:
            vars: List of variable names in the expression.

        Returns:
            A Numba JIT-compiled callable.

        Example:
            f = Expr(...).to_numba(["x"])
            f(np.linspace(0, 1, 1000000))
        """
        return to_numba(self._expr, vars)

    def simplify(self) -> "Expr":
        """Return a simplified version of this expression."""
        return Expr(_get_client().simplify(self._expr))

    def json(self) -> dict:
        """Return the raw Expr JSON dict."""
        return self._expr

    def __repr__(self) -> str:
        return f"Expr({self.pretty()})"


def diff_expr(expr: dict, var: str) -> Expr:
    """Differentiate an Expr JSON dict with respect to var.

    Args:
        expr: An Expr JSON dict.
        var:  Variable name.

    Returns:
        An Expr wrapping the derivative.
    """
    return Expr(_get_client().diff(expr, var))


def partial_diff_expr(expr: dict, var: str) -> Expr:
    """Partial derivative of an Expr JSON dict with respect to var.

    Args:
        expr: An Expr JSON dict.
        var:  Variable name.

    Returns:
        An Expr wrapping the partial derivative.
    """
    return Expr(_get_client().partial_diff(expr, var))


def integrate_expr(expr: dict, var: str) -> Expr:
    """Integrate an Expr JSON dict with respect to var.

    Args:
        expr: An Expr JSON dict.
        var:  Variable name.

    Returns:
        An Expr wrapping the antiderivative.

    Raises:
        LimCalcError: if the integral is non-elementary.
    """
    return Expr(_get_client().integrate(expr, var))


def limit_expr(expr: dict, var: str, x0: float) -> float:
    """Compute lim_{var -> x0} expr.

    Args:
        expr: An Expr JSON dict.
        var:  Variable name.
        x0:   Limit point.

    Returns:
        The limit value as a float.

    Raises:
        LimCalcError: if the limit is a pole or does not exist.
    """
    return _get_client().limit(expr, var, x0)


def gradient_expr(expr: dict, vars: list[str]) -> list[Expr]:
    """Compute the gradient of an Expr JSON dict.

    Args:
        expr: An Expr JSON dict (scalar-valued).
        vars: List of variable names.

    Returns:
        A list of Expr objects [df/dx1, df/dx2, ...].
    """
    results = _get_client().gradient(expr, vars)
    return [Expr(r) for r in results]


# ---------------------------------------------------------------------------
# String-input API (requires LaTeX parser — coming soon)
# ---------------------------------------------------------------------------

def _parse(expr_str: str) -> dict:
    """Parse a string expression to Expr JSON.

    Currently raises NotImplementedError — the LaTeX parser is not
    yet implemented. Use the _expr variants with JSON dicts directly.
    """
    raise NotImplementedError(
        "String expression parsing is not yet implemented. "
        "Use diff_expr(), integrate_expr(), etc. with Expr JSON dicts directly, "
        "or set the LIMCALC_CLI environment variable and use LimCalc() directly."
    )


def diff(expr, var: str) -> Expr:
    """Differentiate expr with respect to var.

    Args:
        expr: An Expr JSON dict or a string expression (coming soon).
        var:  Variable name.

    Returns:
        An Expr wrapping the derivative.
    """
    if isinstance(expr, str):
        return diff_expr(_parse(expr), var)
    return diff_expr(expr, var)


def partial_diff(expr, var: str) -> Expr:
    """Partial derivative of expr with respect to var.

    Args:
        expr: An Expr JSON dict or a string expression (coming soon).
        var:  Variable name.

    Returns:
        An Expr wrapping the partial derivative.
    """
    if isinstance(expr, str):
        return partial_diff_expr(_parse(expr), var)
    return partial_diff_expr(expr, var)


def integrate(expr, var: str) -> Expr:
    """Integrate expr with respect to var.

    Args:
        expr: An Expr JSON dict or a string expression (coming soon).
        var:  Variable name.

    Returns:
        An Expr wrapping the antiderivative.

    Raises:
        LimCalcError: if the integral is non-elementary.
    """
    if isinstance(expr, str):
        return integrate_expr(_parse(expr), var)
    return integrate_expr(expr, var)


def limit(expr, var: str, x0: float) -> float:
    """Compute lim_{var -> x0} expr.

    Args:
        expr: An Expr JSON dict or a string expression (coming soon).
        var:  Variable name.
        x0:   Limit point.

    Returns:
        The limit value as a float.

    Raises:
        LimCalcError: if the limit is a pole or does not exist.
    """
    if isinstance(expr, str):
        return limit_expr(_parse(expr), var, x0)
    return limit_expr(expr, var, x0)


def gradient(expr, vars: list[str]) -> list[Expr]:
    """Compute the gradient of expr.

    Args:
        expr: An Expr JSON dict or a string expression (coming soon).
        vars: List of variable names.

    Returns:
        A list of Expr objects.
    """
    if isinstance(expr, str):
        return gradient_expr(_parse(expr), vars)
    return gradient_expr(expr, vars)


def simplify(expr) -> Expr:
    """Simplify expr algebraically.

    Args:
        expr: An Expr JSON dict.

    Returns:
        A simplified Expr.
    """
    if isinstance(expr, dict):
        return Expr(_get_client().simplify(expr))
    raise TypeError(f"Expected dict, got {type(expr)}")
