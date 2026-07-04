"""
LaTeX math parser for limcalc.

Converts LaTeX math strings to Expr JSON dicts or numeric values,
dispatching to the limcalc Haskell engine as needed.

= Architecture

1. Preprocessor: replaces LaTeX function commands (\sin, \exp, etc.)
   with placeholder tokens (SIN_, EXP_, etc.) to avoid lexer ambiguity.

2. Preprocessor: detects derivative/integral/limit operators with regex
   and dispatches to the limcalc engine directly.

3. Lark grammar: parses pure mathematical expressions using placeholder
   tokens for function names.

4. Transformer: converts the parse tree to Expr JSON dicts.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

from lark import Lark, Transformer, Token
from lark.exceptions import UnexpectedInput, UnexpectedEOF

from limcalc.client import LimCalc, LimCalcError

# ---------------------------------------------------------------------------
# Grammar loading
# ---------------------------------------------------------------------------

_GRAMMAR_PATH = Path(__file__).parent / "grammar.lark"

_parser = Lark(
    _GRAMMAR_PATH.read_text(),
    parser="earley",
    ambiguity="resolve",
)

# ---------------------------------------------------------------------------
# Module-level client (lazy)
# ---------------------------------------------------------------------------

_client: LimCalc | None = None


def _get_client() -> LimCalc:
    global _client
    if _client is None or (_client._proc and _client._proc.poll() is not None):
        _client = LimCalc()
    return _client


# ---------------------------------------------------------------------------
# Expr JSON constructors
# ---------------------------------------------------------------------------

def _var(name):   return {"tag": "Var", "name": name}
def _const(v):    return {"tag": "Const", "value": float(v)}
def _pi():        return {"tag": "Pi"}
def _e():         return {"tag": "E"}
def _i():         return {"tag": "I"}
def _add(a, b):   return {"tag": "Add", "left": a, "right": b}
def _sub(a, b):   return {"tag": "Sub", "left": a, "right": b}
def _mul(a, b):   return {"tag": "Mul", "left": a, "right": b}
def _div(a, b):   return {"tag": "Div", "left": a, "right": b}
def _pow(b, e):   return {"tag": "Pow", "base": b, "exp": e}
def _neg(a):      return {"tag": "Neg", "arg": a}
def _abs(a):      return {"tag": "Abs", "arg": a}
def _fn(tag, a):  return {"tag": tag, "arg": a}
def _sqrt(a):     return _pow(a, _const(0.5))
def _tan(a):      return _div(_fn("Sin", a), _fn("Cos", a))
def _sec(a):      return _div(_const(1), _fn("Cos", a))
def _csc(a):      return _div(_const(1), _fn("Sin", a))
def _cot(a):      return _div(_fn("Cos", a), _fn("Sin", a))
def _sinc(a):     return _div(_fn("Sin", a), a)
def _sinh(a):     return _div(_sub(_fn("Exp", a), _fn("Exp", _neg(a))), _const(2))
def _cosh(a):     return _div(_add(_fn("Exp", a), _fn("Exp", _neg(a))), _const(2))
def _tanh(a):     return _div(_sinh(a), _cosh(a))
def _arcsinh(a):  return _fn("Log", _add(a, _sqrt(_add(_pow(a, _const(2)), _const(1)))))
def _arccosh(a):  return _fn("Log", _add(a, _sqrt(_sub(_pow(a, _const(2)), _const(1)))))
def _arctanh(a):  return _mul(_const(0.5), _fn("Log", _div(_add(_const(1), a), _sub(_const(1), a))))


# ---------------------------------------------------------------------------
# Function command preprocessor
# ---------------------------------------------------------------------------

# Map LaTeX function commands to placeholder tokens.
# Sorted longest-first so \arcsin is replaced before \sin.
_FUNC_CMDS = [
    (r"\operatorname{erf}",  "@F00@"),
    (r"\operatorname{sinc}", "@F01@"),
    (r"\operatorname{Si}",   "@F02@"),
    (r"\operatorname{Ci}",   "@F03@"),
    (r"\operatorname{Ei}",   "@F04@"),
    (r"\operatorname{li}",   "@F05@"),
    (r"\mathrm{erf}",        "@F00@"),
    (r"\mathrm{Si}",         "@F02@"),
    (r"\mathrm{Ci}",         "@F03@"),
    (r"\mathrm{Ei}",         "@F04@"),
    (r"\mathrm{li}",         "@F05@"),
    (r"\text{erf}",          "@F00@"),
    (r"\text{sinc}",         "@F01@"),
    (r"\text{Si}",           "@F02@"),
    (r"\text{Ci}",           "@F03@"),
    (r"\text{Ei}",           "@F04@"),
    (r"\text{li}",           "@F05@"),
    (r"\arcsinh",            "@F06@"),
    (r"\arccosh",            "@F07@"),
    (r"\arctanh",            "@F08@"),
    (r"\arcsin",             "@F09@"),
    (r"\arccos",             "@F10@"),
    (r"\arctan",             "@F11@"),
    (r"\asin",               "@F09@"),
    (r"\acos",               "@F10@"),
    (r"\atan",               "@F11@"),
    (r"\sinh",               "@F12@"),
    (r"\cosh",               "@F13@"),
    (r"\tanh",               "@F14@"),
    (r"\sqrt",               "@F15@"),
    (r"\sin",                "@F16@"),
    (r"\cos",                "@F17@"),
    (r"\tan",                "@F18@"),
    (r"\sec",                "@F19@"),
    (r"\csc",                "@F20@"),
    (r"\cot",                "@F21@"),
    (r"\exp",                "@F22@"),
    (r"\log",                "@F23@"),
    (r"\ln",                 "@F23@"),
]

# Dispatch table: placeholder -> builder
_FN_DISPATCH = {
    "@F00@": lambda a: _fn("Erf", a),
    "@F01@": lambda a: _sinc(a),
    "@F02@": lambda a: _fn("Si", a),
    "@F03@": lambda a: _fn("Ci", a),
    "@F04@": lambda a: _fn("Ei", a),
    "@F05@": lambda a: _fn("Li", a),
    "@F06@": lambda a: _arcsinh(a),
    "@F07@": lambda a: _arccosh(a),
    "@F08@": lambda a: _arctanh(a),
    "@F09@": lambda a: _fn("Arcsin", a),
    "@F10@": lambda a: _fn("Arccos", a),
    "@F11@": lambda a: _fn("Arctan", a),
    "@F12@": lambda a: _sinh(a),
    "@F13@": lambda a: _cosh(a),
    "@F14@": lambda a: _tanh(a),
    "@F15@": lambda a: _sqrt(a),
    "@F16@": lambda a: _fn("Sin", a),
    "@F17@": lambda a: _fn("Cos", a),
    "@F18@": lambda a: _tan(a),
    "@F19@": lambda a: _sec(a),
    "@F20@": lambda a: _csc(a),
    "@F21@": lambda a: _cot(a),
    "@F22@": lambda a: _fn("Exp", a),
    "@F23@": lambda a: _fn("Log", a),
}


def _protect_funcs(s: str) -> str:
    """Replace LaTeX function commands with placeholder tokens."""
    for cmd, placeholder in _FUNC_CMDS:
        s = s.replace(cmd, placeholder)
    return s


# ---------------------------------------------------------------------------
# Greek name mapping
# ---------------------------------------------------------------------------

_GREEK = {
    r"\alpha": "alpha", r"\beta": "beta", r"\gamma": "gamma",
    r"\delta": "delta", r"\epsilon": "epsilon", r"\varepsilon": "varepsilon",
    r"\zeta": "zeta", r"\eta": "eta", r"\theta": "theta",
    r"\vartheta": "vartheta", r"\iota": "iota", r"\kappa": "kappa",
    r"\lambda": "lambda", r"\mu": "mu", r"\nu": "nu", r"\xi": "xi",
    r"\rho": "rho", r"\varrho": "varrho", r"\sigma": "sigma",
    r"\varsigma": "varsigma", r"\tau": "tau", r"\upsilon": "upsilon",
    r"\phi": "phi", r"\varphi": "varphi", r"\chi": "chi",
    r"\psi": "psi", r"\omega": "omega",
    r"\Gamma": "Gamma", r"\Delta": "Delta", r"\Theta": "Theta",
    r"\Lambda": "Lambda", r"\Xi": "Xi", r"\Pi": "Pi",
    r"\Sigma": "Sigma", r"\Upsilon": "Upsilon", r"\Phi": "Phi",
    r"\Psi": "Psi", r"\Omega": "Omega",
}


def _greek(s):
    return _GREEK.get(s, s.lstrip("\\"))


def _varname_str(expr):
    if isinstance(expr, dict) and expr.get("tag") == "Var":
        return expr["name"]
    return str(expr)


def _expr_to_float(expr):
    if isinstance(expr, (int, float)):
        return float(expr)
    if isinstance(expr, dict):
        tag = expr["tag"]
        if tag == "Const":  return expr["value"]
        if tag == "Pi":     return math.pi
        if tag == "E":      return math.e
        if tag == "Neg":    return -_expr_to_float(expr["arg"])
        if tag == "Add":    return _expr_to_float(expr["left"]) + _expr_to_float(expr["right"])
        if tag == "Sub":    return _expr_to_float(expr["left"]) - _expr_to_float(expr["right"])
        if tag == "Mul":    return _expr_to_float(expr["left"]) * _expr_to_float(expr["right"])
        if tag == "Div":    return _expr_to_float(expr["left"]) / _expr_to_float(expr["right"])
        if tag == "Pow":    return _expr_to_float(expr["base"]) ** _expr_to_float(expr["exp"])
    raise ValueError(f"Cannot evaluate to float: {expr}")


# ---------------------------------------------------------------------------
# Preprocessor — detect derivative/integral/limit before parsing
# ---------------------------------------------------------------------------

_VAR_PAT = r"([a-zA-Z]|\\(?:alpha|beta|gamma|delta|epsilon|varepsilon|zeta|eta|theta|vartheta|iota|kappa|lambda|mu|nu|xi|rho|varrho|sigma|varsigma|tau|upsilon|phi|varphi|chi|psi|omega|Gamma|Delta|Theta|Lambda|Xi|Pi|Sigma|Upsilon|Phi|Psi|Omega))"

_DIFF_PATTERNS = [
    re.compile(
        r"^\\frac\{(?:d|\\partial)\}\{(?:d|\\partial)\s*" + _VAR_PAT + r"\}\s*(.+)$",
        re.DOTALL
    ),
    re.compile(
        r"^(?:d|\\partial)\s*/\s*(?:d|\\partial)\s*" + _VAR_PAT + r"\s+(.+)$",
        re.DOTALL
    ),
]

_INT_INDEF = re.compile(
    r"^\\int\s+(.+?)\s+d" + _VAR_PAT + r"\s*$",
    re.DOTALL
)
_INT_DEF = re.compile(
    r"^\\int_\{?(.+?)\}?\^\{?(.+?)\}?\s+(.+?)\s+d" + _VAR_PAT + r"\s*$",
    re.DOTALL
)
_LIM_PATTERN = re.compile(
    r"^\\lim_\{" + _VAR_PAT + r"\s*\\to\s*(.+?)(\^\+|\^-)?\}\s*(.+)$",
    re.DOTALL
)


def _parse_expr(s: str) -> dict:
    """Parse a pure expression string to Expr JSON."""
    s = s.strip()
    s = _protect_funcs(s)
    try:
        tree = _parser.parse(s)
    except (UnexpectedInput, UnexpectedEOF) as e:
        raise ParseError(f"Could not parse: {s!r}\n{e}") from e
    return LatexTransformer().transform(tree)


def _detect_and_dispatch(latex: str):
    """Detect derivative/integral/limit and dispatch to engine."""
    s = latex.strip()

    # Derivative
    for pat in _DIFF_PATTERNS:
        m = pat.match(s)
        if m:
            var_token = m.group(1)
            expr_str  = m.group(2)
            var_name  = _greek(var_token) if var_token.startswith("\\") else var_token
            expr      = _parse_expr(expr_str)
            return _get_client().diff(expr, var_name)

    # Definite integral
    m = _INT_DEF.match(s)
    if m:
        lower_str, upper_str, integrand_str, var_token = m.groups()
        var_name  = _greek(var_token) if var_token.startswith("\\") else var_token
        integrand = _parse_expr(integrand_str)
        antideriv = _get_client().integrate(integrand, var_name)
        upper_val = _get_client().limit(antideriv, var_name, _expr_to_float(_parse_expr(upper_str)))
        lower_val = _get_client().limit(antideriv, var_name, _expr_to_float(_parse_expr(lower_str)))
        return upper_val - lower_val

    # Indefinite integral
    m = _INT_INDEF.match(s)
    if m:
        integrand_str, var_token = m.groups()
        if not integrand_str.strip():
            raise ParseError("Missing integrand")
        var_name  = _greek(var_token) if var_token.startswith("\\") else var_token
        integrand = _parse_expr(integrand_str)
        return _get_client().integrate(integrand, var_name)

    # Limit
    m = _LIM_PATTERN.match(s)
    if m:
        var_token, point_str, direction, expr_str = m.groups()
        var_name = _greek(var_token) if var_token.startswith("\\") else var_token
        x0       = _expr_to_float(_parse_expr(point_str))
        expr     = _parse_expr(expr_str)
        return _get_client().limit(expr, var_name, x0)

    return None


# ---------------------------------------------------------------------------
# Transformer
# ---------------------------------------------------------------------------

class LatexTransformer(Transformer):

    # Passthrough rules
    def start(self, c):         return c[0]
    def expr(self, c):          return c[0]
    def term(self, c):          return c[0]
    def factor(self, c):        return c[0]
    def implicit_mul_pass(self, c): return c[0]
    def power(self, c):         return c[0]
    def postfix(self, c):       return c[0]
    def atom(self, c):          return c[0]
    def func_app(self, c):      return c[0]
    def eval_at(self, c):       return c[0]

    # Arithmetic
    def add(self, c):           return _add(c[0], c[1])
    def sub(self, c):           return _sub(c[0], c[1])
    def mul(self, c):           return _mul(c[0], c[1])
    def div(self, c):           return _div(c[0], c[1])
    def frac(self, c):          return _div(c[0], c[1])
    def neg(self, c):           return _neg(c[0])
    def pos(self, c):           return c[0]

    def implicit_mul(self, c):
        result = c[0]
        for arg in c[1:]:
            result = _mul(result, arg)
        return result

    # Powers
    def pow_braced(self, c):
        base, exp = c[0], c[1]
        if base == _e(): return _fn("Exp", exp)
        return _pow(base, exp)

    def pow_simple(self, c):
        base, exp = c[0], c[1]
        if base == _e(): return _fn("Exp", exp)
        return _pow(base, exp)

    # Grouping — find the expr child by type
    def grouped_lr(self, c):
        for child in c:
            if isinstance(child, dict): return child
        return c[2]

    def grouped_paren(self, c):   return c[0]
    def grouped_bracket(self, c): return c[0]
    def brace(self, c):           return c[0]

    def abs_val(self, c):
        for child in c:
            if isinstance(child, dict): return _abs(child)
        return _abs(c[1])

    def norm(self, c):
        for child in c:
            if isinstance(child, dict): return _abs(child)
        return _abs(c[1])

    # Evaluation notation
    def eval_at_braced(self, c):
        expr, varname, value = c[0], c[2], c[3]
        return _get_client().limit(expr, _varname_str(varname), _expr_to_float(value))

    def eval_at_inline(self, c):
        expr, varname, value = c[0], c[2], c[3]
        return _get_client().limit(expr, _varname_str(varname), _expr_to_float(value))

    def eval_bounds(self, c):
        return c[0]

    # Function application — FUNC_NAME token dispatches via _FN_DISPATCH
    def func_braced(self, c):
        fn_name = str(c[0])
        return _FN_DISPATCH[fn_name](c[1])

    def func_paren(self, c):
        fn_name = str(c[0])
        return _FN_DISPATCH[fn_name](c[1])

    def func_grouped(self, c):
        fn_name = str(c[0])
        return _FN_DISPATCH[fn_name](c[1])

    def func_atom(self, c):
        fn_name = str(c[0])
        return _FN_DISPATCH[fn_name](c[1])

    # Atoms
    def number(self, c):
        return _const(float(str(c[0])))

    def variable(self, c):
        return c[0]

    def constant(self, c):
        s = str(c[0])
        if s == r"\pi":                 return _pi()
        if s in ("e", r"\e"):           return _e()
        if s in ("i", r"\i"):           return _i()
        if s == r"\infty":              return _const(float("inf"))
        return _const(float(s))

    # Variable names
    def varname_simple(self, c):
        s = str(c[0])
        return _var(_greek(s) if s.startswith("\\") else s)

    def varname_sub(self, c):
        return _var(f"{_tok(c[0])}_{_tok(c[1])}")

    def varname_sub_num(self, c):
        return _var(f"{_tok(c[0])}_{c[1]}")

    def varname_body(self, c):
        return "".join(str(t) for t in c)

    def varname(self, c):
        return c[0]


def _tok(t) -> str:
    if isinstance(t, dict) and t.get("tag") == "Var":
        return t["name"]
    s = str(t)
    return _greek(s) if s.startswith("\\") else s


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

class ParseError(Exception):
    """Raised when the LaTeX input cannot be parsed."""
    pass


def parse(latex: str) -> Any:
    """Parse a LaTeX math string and evaluate it.

    Returns an Expr JSON dict for symbolic expressions, or a float
    for limits and evaluated expressions.

    Args:
        latex: A LaTeX math string.

    Returns:
        Expr JSON dict or float.

    Raises:
        ParseError: if the input cannot be parsed.
        LimCalcError: if the engine returns an error.
    """
    if not latex or not latex.strip():
        raise ParseError("Empty input")

    latex = latex.strip()

    # Validate: d/dx with no expression should raise
    for pat in _DIFF_PATTERNS:
        m = pat.match(latex)
        if m and not m.group(2).strip():
            raise ParseError("Missing expression after derivative operator")

    # Validate: \int with no dx should raise
    if latex.startswith(r"\int") and not re.search(r"\sd[a-zA-Z\\]", latex):
        raise ParseError("Missing integration variable (expected 'dx', 'dt', etc.)")

    try:
        result = _detect_and_dispatch(latex)
        if result is not None:
            return result
    except LimCalcError:
        raise
    except ParseError:
        raise
    except Exception as e:
        raise ParseError(f"Error in dispatch: {e}") from e

    return _parse_expr(latex)
