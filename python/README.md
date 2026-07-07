# limcalc

A symbolic calculus engine built on the limit definition rather than rewriting rules.

Differentiation, integration, and limit evaluation are all derived from a single underlying representation — the log-Puiseux series — chosen because it is the natural object for expanding f(x+h) in powers of h. The derivative is the coefficient of h¹ in that expansion. The product rule and chain rule are not implemented; both are consequences of the limit definition and are observable in the output.

String input accepts LaTeX math notation.

## Installation

```
pip install limcalc
```

Wheels are available for Windows, Linux, and macOS (x86-64 and arm64). The package bundles the `limcalc-cli` Haskell binary — no separate installation required.

## Usage

```python
import limcalc as lc

# Differentiate
f = lc.diff(r"\sin(x^2)", "x")
print(f.pretty())      # 2*x*cos(x^2)
print(f.to_numpy())    # 2 * x * np.cos(x ** 2)

# Integrate
g = lc.integrate(r"e^{-x^2}", "x")
print(g.pretty())      # (sqrt(π)/2)*erf(x)

# Limits
print(lc.limit(r"\frac{\sin(x)}{x}", "x", 0))  # 1.0

# Partial derivatives and gradient
grad = lc.gradient(r"x^2 + y^2", ["x", "y"])
print([c.pretty() for c in grad])  # ['2*x', '2*y']
```

### LaTeX input

The parser accepts standard LaTeX math notation:

```python
# Greek variables
lc.diff(r"\sin(\theta)", "theta")

# Fractions
lc.diff(r"\frac{x^2 + 1}{x}", "x")

# Derivative notation (dispatches directly)
lc.parse(r"\frac{d}{dx} \sin(x)")       # cos(x)
lc.parse(r"\int \sin(x) \, dx")         # -cos(x)
lc.parse(r"\lim_{x \to 0} \frac{\sin(x)}{x}")  # 1.0
```

### Code generation

```python
import numpy as np

f = lc.diff(r"\sin(x^2)", "x")

# NumPy string
f.to_numpy()   # "2 * x * np.cos(x ** 2)"

# Callable
fn = f.to_lambda(["x"])
fn(np.linspace(0, np.pi, 100))

# Numba JIT (requires pip install limcalc[numba])
fast = f.to_numba(["x"])
```

### Non-elementary integrals

When no elementary antiderivative exists, limcalc returns a proved certificate rather than a failure:

```python
from limcalc import LimCalcError

try:
    lc.integrate(r"e^{x^2}", "x")
except LimCalcError as e:
    print(e)  # NonElementary: no elementary antiderivative exists
```

Known non-elementary integrals are recognized and returned as named special functions:

```python
lc.integrate(r"e^{-x^2}", "x").pretty()   # (sqrt(π)/2)*erf(x)
lc.integrate(r"\frac{\sin(x)}{x}", "x").pretty()   # Si(x)
lc.integrate(r"\frac{1}{\log(x)}", "x").pretty()   # li(x)
```

### Low-level JSON API

All functions also accept Expr JSON dicts directly, bypassing the parser:

```python
expr = {"tag": "Sin", "arg": {"tag": "Var", "name": "x"}}
lc.diff(expr, "x")  # {"tag": "Cos", "arg": {"tag": "Var", "name": "x"}}
```

## Supported functions

| LaTeX | Notes |
|---|---|
| `\sin`, `\cos`, `\tan`, `\sec`, `\csc`, `\cot` | |
| `\arcsin`, `\arccos`, `\arctan` | |
| `\sinh`, `\cosh`, `\tanh` | Expanded via exp |
| `\exp`, `\log`, `\ln` | |
| `\sqrt{x}` | Expanded as x^(1/2) |
| `\operatorname{erf}(x)` | Also `\mathrm{erf}`, `\text{erf}` |
| `\text{Si}(x)`, `\text{Ci}(x)` | Sine and cosine integrals |
| `\text{Ei}(x)`, `\text{li}(x)` | Exponential and logarithmic integrals |

## Building from source

The Haskell CLI source is at [github.com/penny4nonsense/limcalc](https://github.com/penny4nonsense/limcalc).

```
git clone https://github.com/penny4nonsense/limcalc
cd limcalc
cabal build limcalc-cli
```

Then point the Python package at your binary:

```
export LIMCALC_CLI=/path/to/limcalc-cli
pip install limcalc --no-binary limcalc
```

## Background

limcalc is described in the paper [Calculus from First Principles: A Limit-Based Symbolic Computation Engine via Log-Puiseux Series](https://github.com/penny4nonsense/limcalc). The core idea is that the log-Puiseux series — finite formal sums of terms c·hᵖ·log(h)ᵏ — is the correct representation class for the local behavior of elementary functions near any point, and that choosing it as the computational foundation makes the limit definition executable rather than motivational.

## License

MIT
