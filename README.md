# limcalc

A symbolic calculus engine built on the limit definition rather than rewriting rules.

Differentiation, integration, and limit evaluation are all derived from a single underlying representation — the log-Puiseux series — finite formal sums of terms `c · hᵖ · log(h)ᵏ` with `c` an algebraic number, `p` rational, and `k` a natural number. The derivative is the coefficient of `h¹` in the expansion of `f(x+h)`. The product rule and chain rule are not implemented; both are consequences of the limit definition and are observable in the output.

## Quick start

```haskell
import LimCalc.Core.Expr
import LimCalc.Differentiation.Calculus
import LimCalc.Differentiation.Limit
import LimCalc.Integration.Risch

-- Symbolic expressions
let f = Sin (Var "x")

-- Symbolic differentiation
diff f "x"
-- Right (Cos (Var "x"))

-- Gradient
gradient (Add (Pow (Var "x") (Lit 2)) (Pow (Var "y") (Lit 2))) ["x", "y"]
-- Right [Mul (Lit 2) (Var "x"), Mul (Lit 2) (Var "y")]

-- Limits
limit (Div (Sin (Var "x")) (Var "x")) "x" 0.0
-- Exists 1.0

limit (Div (Lit 1) (Var "x")) "x" 0.0
-- Pole (-1.0)

-- Integration
rischIntegrate (Sin (Var "x")) "x"
-- Elementary (Neg (Cos (Var "x")))

rischIntegrate (Exp (Neg (Pow (Var "x") (Lit 2)))) "x"
-- SpecialFunction (Mul (Div (Sqrt Pi) (Lit 2)) (Erf (Var "x")))

rischIntegrate (Exp (Pow (Var "x") (Lit 2))) "x"
-- NonElementary
```

## API

### Differentiation

```haskell
-- Symbolic derivative (algebraic path)
diff :: Expr -> String -> Either ExpandError Expr

-- Partial derivative
partialDiff :: Expr -> String -> Either ExpandError Expr

-- Gradient, Jacobian, Hessian
gradient :: Expr -> [String] -> Either ExpandError [Expr]
jacobian :: [Expr] -> [String] -> Either ExpandError [[Expr]]
hessian  :: Expr -> [String] -> Either ExpandError [[Expr]]

-- Numeric derivative at a point
derivative    :: Expr -> Map String Double -> String -> Either ExpandError Double
partial       :: Expr -> Map String Double -> String -> Either ExpandError Double
nthDerivative :: Int -> Expr -> Map String Double -> String -> Either ExpandError Double
```

### Limits

```haskell
-- One-sided and two-sided limits
limit      :: Expr -> String -> Double -> LimitResult Double
limitRight :: Expr -> String -> Double -> LimitResult Double
limitLeft  :: Expr -> String -> Double -> LimitResult Double
```

`LimitResult` is a certified result type:

```haskell
data LimitResult a
  = Exists a            -- limit exists and equals a
  | Pole Double         -- diverges; leading log-Puiseux exponent
  | DoesNotExist String -- certified non-existence with witness
  | LimitError String   -- outside the representation class
```

### Integration

```haskell
rischIntegrate :: Expr -> String -> RischResult
```

`RischResult` distinguishes proved outcomes:

```haskell
data RischResult
  = Elementary Expr       -- elementary antiderivative
  | SpecialFunction Expr  -- named special function (erf, Si, Ci, Ei, li)
  | NonElementary         -- proved: no elementary antiderivative exists
```

`NonElementary` is a mathematical certificate, not a failure. It is returned only when the Risch algorithm has determined non-existence via Liouville's theorem — never on timeout or truncation.

### Supported special functions

The following non-elementary integrals are recognized and returned as named constructors:

| Integrand | Result |
|---|---|
| `exp(-x²)` | `(√π/2) · erf(x)` |
| `sin(x)/x` | `Si(x)` |
| `cos(x)/x` | `Ci(x)` |
| `exp(x)/x` | `Ei(x)` |
| `1/log(x)` | `li(x)` |

## CLI

The package includes `limcalc-cli`, a JSON line protocol server suitable for use from other languages:

```
$ limcalc-cli
{"op":"diff","expr":{"tag":"Sin","arg":{"tag":"Var","name":"x"}},"var":"x"}
{"ok":true,"result":{"tag":"Cos","arg":{"tag":"Var","name":"x"}}}
{"op":"limit","expr":{"tag":"Div","left":{"tag":"Sin","arg":{"tag":"Var","name":"x"}},"right":{"tag":"Var","name":"x"}},"var":"x","x0":0}
{"ok":true,"result":1.0}
```

A Python interface wrapping this CLI is available as [limcalc-py](https://pypi.org/project/limcalc-py/) on PyPI.

## Background

The log-Puiseux type is strictly more expressive than pure Puiseux series. The cosine integral `Ci(x)` near `x = 0` expands as `γ + log(h) - h²/4 + ⋯` — the `log(h)` term cannot be represented as a finite sum of powers `hᵖ`. The system proves this formally (Lemma 2.1 in the paper) and uses it to motivate the representation class.

The boundary of the class is sharp: `li(x) = Ei(log x)` near `x = 0` requires `log(log(h))` terms, which lie outside the single-logarithm tower. The system correctly computes `d/dx li(x) = 1/log(x)` via the algebraic path while returning `LimitError` for the series path — an honest boundary, not an error.

## License

MIT
