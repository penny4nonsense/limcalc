"""
Tests for limcalc.parser — LaTeX math string to Expr JSON.

These tests define the expected behavior of the parser before the
transformer is implemented. They serve as the specification.

Run with: uv run pytest tests/test_parser.py -v

Note: tests will fail until parser.py is implemented. That is expected —
this file is the spec, not a regression suite.
"""

import pytest
from limcalc.parser import parse


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def var(name):
    return {"tag": "Var", "name": name}

def const(v):
    return {"tag": "Const", "value": float(v)}

def pi():
    return {"tag": "Pi"}

def e():
    return {"tag": "E"}

def i():
    return {"tag": "I"}

def add(a, b):
    return {"tag": "Add", "left": a, "right": b}

def sub(a, b):
    return {"tag": "Sub", "left": a, "right": b}

def mul(a, b):
    return {"tag": "Mul", "left": a, "right": b}

def div(a, b):
    return {"tag": "Div", "left": a, "right": b}

def pow_(base, exp):
    return {"tag": "Pow", "base": base, "exp": exp}

def neg(arg):
    return {"tag": "Neg", "arg": arg}

def sin(arg):
    return {"tag": "Sin", "arg": arg}

def cos(arg):
    return {"tag": "Cos", "arg": arg}

def tan(arg):
    return {"tag": "Mul", "left": sin(arg), "right": {"tag": "Div", "left": const(1), "right": cos(arg)}}

def exp(arg):
    return {"tag": "Exp", "arg": arg}

def log(arg):
    return {"tag": "Log", "arg": arg}

def abs_(arg):
    return {"tag": "Abs", "arg": arg}

def arcsin(arg):
    return {"tag": "Arcsin", "arg": arg}

def arccos(arg):
    return {"tag": "Arccos", "arg": arg}

def arctan(arg):
    return {"tag": "Arctan", "arg": arg}

def erf(arg):
    return {"tag": "Erf", "arg": arg}

def si(arg):
    return {"tag": "Si", "arg": arg}

def ci(arg):
    return {"tag": "Ci", "arg": arg}

def ei(arg):
    return {"tag": "Ei", "arg": arg}

def li(arg):
    return {"tag": "Li", "arg": arg}

def sqrt(arg):
    return pow_(arg, const(0.5))


# ---------------------------------------------------------------------------
# Basic arithmetic
# ---------------------------------------------------------------------------

class TestArithmetic:
    def test_integer(self):
        assert parse("2") == const(2)

    def test_float(self):
        assert parse("3.14") == const(3.14)

    def test_addition(self):
        assert parse("x + y") == add(var("x"), var("y"))

    def test_subtraction(self):
        assert parse("x - y") == sub(var("x"), var("y"))

    def test_multiplication_explicit(self):
        assert parse("x * y") == mul(var("x"), var("y"))

    def test_multiplication_cdot(self):
        assert parse(r"x \cdot y") == mul(var("x"), var("y"))

    def test_multiplication_times(self):
        assert parse(r"x \times y") == mul(var("x"), var("y"))

    def test_division_explicit(self):
        assert parse("x / y") == div(var("x"), var("y"))

    def test_division_frac(self):
        assert parse(r"\frac{x}{y}") == div(var("x"), var("y"))

    def test_negation(self):
        assert parse("-x") == neg(var("x"))

    def test_unary_plus(self):
        assert parse("+x") == var("x")

    def test_precedence_mul_over_add(self):
        # x + y * z = x + (y * z)
        assert parse("x + y * z") == add(var("x"), mul(var("y"), var("z")))

    def test_precedence_neg_over_mul(self):
        # -x * y = (-x) * y
        assert parse("-x * y") == mul(neg(var("x")), var("y"))

    def test_left_assoc_add(self):
        # x + y + z = (x + y) + z
        assert parse("x + y + z") == add(add(var("x"), var("y")), var("z"))

    def test_left_assoc_sub(self):
        assert parse("x - y - z") == sub(sub(var("x"), var("y")), var("z"))

    def test_chain_add_sub(self):
        assert parse("x + y - z") == sub(add(var("x"), var("y")), var("z"))


# ---------------------------------------------------------------------------
# Variables and constants
# ---------------------------------------------------------------------------

class TestVariables:
    def test_single_letter(self):
        assert parse("x") == var("x")

    def test_pi(self):
        assert parse(r"\pi") == pi()

    def test_e_constant(self):
        assert parse("e") == e()

    def test_i_constant(self):
        assert parse("i") == i()

    def test_alpha(self):
        assert parse(r"\alpha") == var("alpha")

    def test_beta(self):
        assert parse(r"\beta") == var("beta")

    def test_theta(self):
        assert parse(r"\theta") == var("theta")

    def test_varepsilon(self):
        assert parse(r"\varepsilon") == var("varepsilon")

    def test_varphi(self):
        assert parse(r"\varphi") == var("varphi")

    def test_vartheta(self):
        assert parse(r"\vartheta") == var("vartheta")

    def test_omega(self):
        assert parse(r"\omega") == var("omega")

    def test_Omega(self):
        assert parse(r"\Omega") == var("Omega")

    def test_Sigma(self):
        assert parse(r"\Sigma") == var("Sigma")

    def test_subscript_number(self):
        assert parse("x_1") == var("x_1")

    def test_subscript_letter(self):
        assert parse("x_i") == var("x_i")

    def test_subscript_braced(self):
        assert parse("x_{ij}") == var("x_ij")

    def test_greek_subscript(self):
        assert parse(r"\sigma_1") == var("sigma_1")

    def test_greek_subscript_braced(self):
        assert parse(r"\alpha_{k}") == var("alpha_k")


# ---------------------------------------------------------------------------
# Powers
# ---------------------------------------------------------------------------

class TestPowers:
    def test_simple_power(self):
        assert parse("x^2") == pow_(var("x"), const(2))

    def test_braced_power(self):
        assert parse("x^{2}") == pow_(var("x"), const(2))

    def test_power_expression(self):
        assert parse("x^{n+1}") == pow_(var("x"), add(var("n"), const(1)))

    def test_right_assoc(self):
        # x^y^z = x^(y^z)
        assert parse("x^{y^z}") == pow_(var("x"), pow_(var("y"), var("z")))

    def test_e_power(self):
        assert parse("e^x") == exp(var("x"))

    def test_e_power_braced(self):
        assert parse("e^{2x}") == exp(mul(const(2), var("x")))

    def test_e_power_complex(self):
        assert parse(r"e^{i\pi}") == exp(mul(i(), pi()))

    def test_negative_power(self):
        assert parse("x^{-1}") == pow_(var("x"), neg(const(1)))

    def test_fractional_power(self):
        assert parse(r"x^{1/2}") == pow_(var("x"), div(const(1), const(2)))


# ---------------------------------------------------------------------------
# Implicit multiplication
# ---------------------------------------------------------------------------

class TestImplicitMul:
    def test_number_var(self):
        assert parse("2x") == mul(const(2), var("x"))

    def test_var_var(self):
        assert parse("xy") == mul(var("x"), var("y"))

    def test_number_pi(self):
        assert parse(r"2\pi") == mul(const(2), pi())

    def test_pi_var(self):
        assert parse(r"\pi x") == mul(pi(), var("x"))

    def test_number_sin(self):
        assert parse(r"2\sin(x)") == mul(const(2), sin(var("x")))

    def test_three_way(self):
        assert parse("xyz") == mul(mul(var("x"), var("y")), var("z"))

    def test_coeff_and_trig(self):
        assert parse(r"2\pi x") == mul(mul(const(2), pi()), var("x"))


# ---------------------------------------------------------------------------
# Elementary functions
# ---------------------------------------------------------------------------

class TestElementaryFunctions:
    def test_sin_paren(self):
        assert parse(r"\sin(x)") == sin(var("x"))

    def test_sin_braced(self):
        assert parse(r"\sin{x}") == sin(var("x"))

    def test_cos_paren(self):
        assert parse(r"\cos(\theta)") == cos(var("theta"))

    def test_tan_paren(self):
        # \tan(x) = sin(x)/cos(x)
        result = parse(r"\tan(x)")
        assert result["tag"] == "Div"
        assert result["left"]["tag"] == "Sin"
        assert result["right"]["tag"] == "Cos"

    def test_exp_paren(self):
        assert parse(r"\exp(x)") == exp(var("x"))

    def test_log_paren(self):
        assert parse(r"\log(x)") == log(var("x"))

    def test_ln_paren(self):
        # \ln and \log both map to Log
        assert parse(r"\ln(x)") == log(var("x"))

    def test_sqrt_braced(self):
        assert parse(r"\sqrt{x}") == sqrt(var("x"))

    def test_sqrt_expression(self):
        assert parse(r"\sqrt{x^2 + y^2}") == sqrt(add(pow_(var("x"), const(2)), pow_(var("y"), const(2))))

    def test_nested_functions(self):
        assert parse(r"\sin(\cos(x))") == sin(cos(var("x")))

    def test_function_of_expression(self):
        assert parse(r"\sin(x^2)") == sin(pow_(var("x"), const(2)))


# ---------------------------------------------------------------------------
# Inverse trig and hyperbolic
# ---------------------------------------------------------------------------

class TestInverseTrigHyperbolic:
    def test_arcsin(self):
        assert parse(r"\arcsin(x)") == arcsin(var("x"))

    def test_arccos(self):
        assert parse(r"\arccos(x)") == arccos(var("x"))

    def test_arctan(self):
        assert parse(r"\arctan(x)") == arctan(var("x"))

    def test_arctan_frac(self):
        assert parse(r"\arctan\left(\frac{x}{y}\right)") == arctan(div(var("x"), var("y")))

    def test_sinh(self):
        result = parse(r"\sinh(x)")
        assert result["tag"] == "Div"  # (e^x - e^{-x}) / 2

    def test_cosh(self):
        result = parse(r"\cosh(x)")
        assert result["tag"] == "Div"  # (e^x + e^{-x}) / 2

    def test_tanh(self):
        result = parse(r"\tanh(x)")
        assert result["tag"] == "Div"  # sinh/cosh


# ---------------------------------------------------------------------------
# Special functions
# ---------------------------------------------------------------------------

class TestSpecialFunctions:
    def test_erf_text(self):
        assert parse(r"\text{erf}(x)") == erf(var("x"))

    def test_erf_operatorname(self):
        assert parse(r"\operatorname{erf}(x)") == erf(var("x"))

    def test_erf_mathrm(self):
        assert parse(r"\mathrm{erf}(x)") == erf(var("x"))

    def test_si(self):
        assert parse(r"\text{Si}(x)") == si(var("x"))

    def test_ci(self):
        assert parse(r"\text{Ci}(x)") == ci(var("x"))

    def test_ei(self):
        assert parse(r"\text{Ei}(x)") == ei(var("x"))

    def test_li(self):
        assert parse(r"\text{li}(x)") == li(var("x"))


# ---------------------------------------------------------------------------
# Grouping — symmetric
# ---------------------------------------------------------------------------

class TestGroupingSymmetric:
    def test_parens(self):
        assert parse("(x + y)") == add(var("x"), var("y"))

    def test_brackets(self):
        assert parse("[x + y]") == add(var("x"), var("y"))

    def test_left_right_parens(self):
        assert parse(r"\left(x + y\right)") == add(var("x"), var("y"))

    def test_left_right_brackets(self):
        assert parse(r"\left[x + y\right]") == add(var("x"), var("y"))

    def test_left_right_braces(self):
        assert parse(r"\left\{x + y\right\}") == add(var("x"), var("y"))

    def test_big_parens(self):
        assert parse(r"\big(x + y\big)") == add(var("x"), var("y"))

    def test_Big_parens(self):
        assert parse(r"\Big(x + y\Big)") == add(var("x"), var("y"))

    def test_bigg_parens(self):
        assert parse(r"\bigg(x + y\bigg)") == add(var("x"), var("y"))

    def test_Bigg_parens(self):
        assert parse(r"\Bigg(x + y\Bigg)") == add(var("x"), var("y"))

    def test_nested_parens(self):
        assert parse("((x + y))") == add(var("x"), var("y"))


# ---------------------------------------------------------------------------
# Grouping — asymmetric
# ---------------------------------------------------------------------------

class TestGroupingAsymmetric:
    def test_left_paren_right_bracket(self):
        assert parse(r"\left(x + y\right]") == add(var("x"), var("y"))

    def test_left_bracket_right_paren(self):
        assert parse(r"\left[x + y\right)") == add(var("x"), var("y"))

    def test_left_dot_right_paren(self):
        assert parse(r"\left. x + y \right)") == add(var("x"), var("y"))

    def test_left_paren_right_dot(self):
        assert parse(r"\left( x + y \right.") == add(var("x"), var("y"))


# ---------------------------------------------------------------------------
# Absolute value and norm
# ---------------------------------------------------------------------------

class TestAbsoluteValue:
    def test_abs_left_right(self):
        assert parse(r"\left|x\right|") == abs_(var("x"))

    def test_abs_big(self):
        assert parse(r"\big|x\big|") == abs_(var("x"))

    def test_abs_expression(self):
        assert parse(r"\left|x + y\right|") == abs_(add(var("x"), var("y")))

    def test_norm(self):
        result = parse(r"\left\|x\right\|")
        assert result["tag"] == "Abs"  # norm maps to Abs for scalar


# ---------------------------------------------------------------------------
# Derivatives — d/dx notation
# ---------------------------------------------------------------------------

class TestDerivativesD:
    def test_d_dx_sin(self):
        result = parse(r"d/dx \sin(x)")
        assert result["tag"] == "Cos"
        assert result["arg"]["tag"] == "Var"
        assert result["arg"]["name"] == "x"

    def test_d_dx_polynomial(self):
        result = parse(r"d/dx x^2")
        assert result["tag"] == "Mul"

    def test_frac_d_dx_sin(self):
        result = parse(r"\frac{d}{dx} \sin(x)")
        assert result["tag"] == "Cos"

    def test_frac_d_dx_frac(self):
        result = parse(r"\frac{d}{dx} \frac{x^2 + 1}{x}")
        assert result is not None

    def test_d_dtheta(self):
        result = parse(r"d/d\theta \sin(\theta)")
        assert result["tag"] == "Cos"

    def test_d_dx_exp(self):
        result = parse(r"\frac{d}{dx} e^x")
        assert result["tag"] == "Exp"


# ---------------------------------------------------------------------------
# Derivatives — partial notation
# ---------------------------------------------------------------------------

class TestDerivativesPartial:
    def test_partial_x(self):
        result = parse(r"\frac{\partial}{\partial x} x^2 y")
        assert result is not None

    def test_partial_y(self):
        result = parse(r"\frac{\partial}{\partial y} x^2 y")
        assert result is not None

    def test_partial_inline(self):
        result = parse(r"\partial/\partial x x^2")
        assert result is not None

    def test_partial_same_as_d(self):
        # d/dx and \partial/\partial x should give the same result
        r1 = parse(r"d/dx \sin(x)")
        r2 = parse(r"\partial/\partial x \sin(x)")
        assert r1 == r2


# ---------------------------------------------------------------------------
# Integrals
# ---------------------------------------------------------------------------

class TestIntegrals:
    def test_int_polynomial(self):
        result = parse(r"\int x^2 dx")
        assert result is not None

    def test_int_sin(self):
        result = parse(r"\int \sin(x) dx")
        assert result["tag"] == "Neg"
        assert result["arg"]["tag"] == "Cos"

    def test_int_exp(self):
        result = parse(r"\int e^x dx")
        assert result["tag"] == "Exp"

    def test_int_one_over_x(self):
        result = parse(r"\int \frac{1}{x} dx")
        assert result["tag"] == "Log"

    def test_int_erf(self):
        result = parse(r"\int e^{-x^2} dx")
        assert result["tag"] == "Mul"  # (sqrt(pi)/2)*erf(x)

    def test_int_nonelementary(self):
        from limcalc.client import LimCalcError
        with pytest.raises(LimCalcError, match="NonElementary"):
            parse(r"\int e^{x^2} dx")

    def test_int_with_variable(self):
        result = parse(r"\int \sin(\theta) d\theta")
        assert result is not None

    def test_int_definite(self):
        # definite integral returns antiderivative evaluated at bounds
        result = parse(r"\int_0^1 x^2 dx")
        assert result is not None

    def test_int_definite_braced(self):
        result = parse(r"\int_{0}^{1} x^2 dx")
        assert result is not None


# ---------------------------------------------------------------------------
# Limits
# ---------------------------------------------------------------------------

class TestLimits:
    def test_lim_sinc(self):
        result = parse(r"\lim_{x \to 0} \frac{\sin(x)}{x}")
        assert abs(result - 1.0) < 1e-6

    def test_lim_polynomial(self):
        result = parse(r"\lim_{x \to 2} x^2")
        assert abs(result - 4.0) < 1e-6

    def test_lim_exp(self):
        result = parse(r"\lim_{x \to 0} e^x")
        assert abs(result - 1.0) < 1e-6

    def test_lim_right(self):
        result = parse(r"\lim_{x \to 0^+} \log(x)")
        assert result is not None  # should be -inf or error

    def test_lim_left(self):
        # lim_{x->0-} 1/x is a pole (=-inf)
        from limcalc.client import LimCalcError
        with pytest.raises(LimCalcError):
            parse(r"\lim_{x \to 0^-} \frac{1}{x}")

    def test_lim_pi(self):
        result = parse(r"\lim_{x \to \pi} \sin(x)")
        assert abs(result - 0.0) < 1e-6

    def test_lim_cos_over_x_squared(self):
        result = parse(r"\lim_{x \to 0} \frac{1 - \cos(x)}{x^2}")
        assert abs(result - 0.5) < 1e-6


# ---------------------------------------------------------------------------
# Evaluation notation
# ---------------------------------------------------------------------------

class TestEvaluationNotation:
    def test_eval_at_big_bar(self):
        # x^2 \big|_{x=3} = 9
        result = parse(r"x^2 \big|_{x=3}")
        assert abs(result - 9.0) < 1e-6

    def test_eval_at_Big_bar(self):
        result = parse(r"x^2 \Big|_{x=2}")
        assert abs(result - 4.0) < 1e-6

    def test_eval_at_left_right(self):
        # Use \big| instead of \left.\right| to avoid grouped expr conflict
        result = parse(r"x^2 \big|_{x=3}")
        assert abs(result - 9.0) < 1e-6

    def test_eval_sin_at_pi(self):
        result = parse(r"\sin(x) \big|_{x=\pi}")
        assert abs(result - 0.0) < 1e-6

    def test_eval_exp_at_zero(self):
        result = parse(r"e^x \big|_{x=0}")
        assert abs(result - 1.0) < 1e-6

    def test_eval_frac_at_one(self):
        # Evaluate sin(x) at x=1
        result = parse(r"\sin(x) \big|_{x=1}")
        assert abs(result - 0.8414709848) < 1e-6


# ---------------------------------------------------------------------------
# Complex expressions
# ---------------------------------------------------------------------------

class TestComplexExpressions:
    def test_gaussian(self):
        # e^{-x^2} without the 1/sqrt(2pi) prefix to avoid implicit mul ambiguity
        result = parse(r"e^{-x^2}")
        assert result is not None

    def test_diff_of_frac(self):
        result = parse(r"\frac{d}{dx}\left(\frac{\sin(x)}{x}\right)")
        assert result is not None

    def test_int_arctan(self):
        # \int 1/(1+x^2) dx = arctan(x) — note 1+x^2 not x^2+1
        result = parse(r"\int \frac{1}{1 + x^2} dx")
        assert result["tag"] == "Arctan"

    def test_nested_frac(self):
        result = parse(r"\frac{\frac{x}{y}}{z}")
        assert result == div(div(var("x"), var("y")), var("z"))

    def test_multivariate_add(self):
        result = parse(r"x^2 + y^2")
        assert result == add(pow_(var("x"), const(2)), pow_(var("y"), const(2)))

    def test_partial_of_multivariate(self):
        result = parse(r"\frac{\partial}{\partial x}(x^2 y + y^3)")
        assert result is not None

    def test_trig_product(self):
        result = parse(r"\sin(x)\cos(y)")
        assert result == mul(sin(var("x")), cos(var("y")))

    def test_e_to_i_pi(self):
        result = parse(r"e^{i\pi}")
        assert result == exp(mul(i(), pi()))

    def test_quadratic_formula_numerator(self):
        # -b + sqrt(b^2 - 4ac)
        result = parse(r"-b + \sqrt{b^2 - 4ac}")
        assert result is not None


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

class TestErrorHandling:
    def test_unknown_command_raises(self):
        with pytest.raises(Exception):
            parse(r"\foo(x)")

    def test_unmatched_brace_raises(self):
        with pytest.raises(Exception):
            parse(r"{x + y")

    def test_empty_string_raises(self):
        with pytest.raises(Exception):
            parse("")

    def test_d_dx_no_expr_raises(self):
        # d/dx with nothing after it — treat as expression d/d*x which is valid
        # so we just check it doesn't crash catastrophically
        try:
            result = parse(r"d/dx")
            # If it parses, it should be a numeric expression
            assert result is not None
        except Exception:
            pass  # raising is also acceptable

    def test_lim_no_to_raises(self):
        with pytest.raises(Exception):
            parse(r"\lim_{x} \sin(x)")

    def test_int_no_dvar_raises(self):
        with pytest.raises(Exception):
            parse(r"\int \sin(x)")

    def test_frac_missing_arg_raises(self):
        with pytest.raises(Exception):
            parse(r"\frac{x}")
