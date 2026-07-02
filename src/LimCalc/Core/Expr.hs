-- | Core expression AST for limcalc.
--
-- 'Expr' is the symbolic representation of mathematical expressions
-- over which all calculus operations are defined. Differentiation,
-- integration, and limit computation all operate by recursion over
-- this tree.
--
-- The constructor set is minimal: only functions that are genuinely
-- primitive in the analytic sense are included. Derived functions
-- (e.g. 'tan'', 'sinh') are expressed in terms of primitives rather
-- than added as constructors, keeping the AST small and the
-- case-analysis in downstream modules tractable.
--
-- 'Const' currently holds a 'Double'; the intention is to replace
-- this with 'AlgNum' once the algebraic number tower is sufficiently
-- complete to handle all coefficient arithmetic symbolically.
module LimCalc.Core.Expr
  ( -- * Core AST
    Expr (..)
    -- * Derived functions
  , tan'
  ) where

-- | The core symbolic expression type.
--
-- Constructors are grouped by role:
--
-- * /Atoms/: 'Const', 'Pi', 'E', 'I', 'Var'
-- * /Arithmetic/: 'Add', 'Sub', 'Mul', 'Div', 'Pow', 'Neg'
-- * /Elementary transcendentals/: 'Exp', 'Log', 'Sin', 'Cos', 'Abs'
-- * /Inverse trig/: 'Arcsin', 'Arccos', 'Arctan'
-- * /Special functions/: 'Erf', 'Li', 'Si', 'Ci', 'Ei'
--
-- Special functions are first-class constructors rather than opaque
-- tokens because the Risch integrator needs to recognize them as
-- integration results, and the differentiation engine needs to
-- compute their derivatives symbolically via 'deriveBase'.
data Expr
  = Const Double
    -- ^ Numeric constant. Currently backed by 'Double'; intended to
    -- be replaced by 'AlgNum' for exact coefficient arithmetic.
  | Pi
    -- ^ The constant π.
  | E
    -- ^ The constant e (base of the natural logarithm).
  | I
    -- ^ The imaginary unit i, satisfying i² = −1.
  | Var String
    -- ^ A named variable. The 'String' is the variable name, matched
    -- against the point map in 'expand' and the variable name in
    -- differentiation.
  | Add Expr Expr
    -- ^ Addition: @Add f g@ represents @f + g@.
  | Sub Expr Expr
    -- ^ Subtraction: @Sub f g@ represents @f - g@.
  | Mul Expr Expr
    -- ^ Multiplication: @Mul f g@ represents @f * g@.
  | Div Expr Expr
    -- ^ Division: @Div f g@ represents @f / g@.
  | Pow Expr Expr
    -- ^ Exponentiation: @Pow f g@ represents @f ^ g@.
  | Neg Expr
    -- ^ Unary negation: @Neg f@ represents @-f@.
  | Abs Expr
    -- ^ Absolute value. Non-analytic at zeros of its argument;
    -- 'expand' returns 'NonAnalytic' when the argument has an
    -- odd-order zero at the expansion point.
  | Exp Expr
    -- ^ The exponential function: @Exp f@ represents @e^f@.
  | Log Expr
    -- ^ Natural logarithm. Undefined for non-positive arguments;
    -- 'expand' returns 'Undefined' when the argument is non-positive
    -- at the expansion point, and produces a @log(h)@ term when the
    -- argument vanishes to first order.
  | Sin Expr
    -- ^ Sine (argument in radians).
  | Cos Expr
    -- ^ Cosine (argument in radians).
  | Arcsin Expr
    -- ^ Inverse sine, with range [−π\/2, π\/2].
  | Arccos Expr
    -- ^ Inverse cosine, with range [0, π].
  | Arctan Expr
    -- ^ Inverse tangent, with range (−π\/2, π\/2).
  | Erf Expr
    -- ^ Error function: @(2\/√π) · ∫₀ˣ e^(−t²) dt@.
    -- Entire; its Puiseux expansion at any point is a pure power series.
  | Li Expr
    -- ^ Logarithmic integral: @∫₀ˣ dt\/log(t)@.
    -- Has a logarithmic singularity at x = 0 and x = 1.
    -- Expansion near x = 0 requires @log(log(h))@, outside the
    -- current log-Puiseux type; left as 'Unknown' in 'expand'.
  | Si Expr
    -- ^ Sine integral: @∫₀ˣ sin(t)\/t dt@.
    -- Entire; its Puiseux expansion at any point is a pure power series.
  | Ci Expr
    -- ^ Cosine integral: @−∫ₓ^∞ cos(t)\/t dt@.
    -- Has a @log(h)@ singularity at x = 0, handled explicitly in 'expand'.
  | Ei Expr
    -- ^ Exponential integral: @−∫₋ₓ^∞ e^(−t)\/t dt@.
    -- Has a @log(h)@ singularity at x = 0, handled explicitly in 'expand'.
  deriving (Show, Eq)

-- | Tangent, expressed as @sin(x) \/ cos(x)@.
--
-- Not a primitive constructor: all analytic behaviour of @tan@ is
-- captured by 'Sin' and 'Cos', so no special cases are needed in
-- 'expand', 'diff', or the Risch integrator.
tan' :: Expr -> Expr
tan' x = Div (Sin x) (Cos x)