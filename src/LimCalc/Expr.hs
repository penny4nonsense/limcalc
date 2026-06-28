module LimCalc.Expr where

-- | Core expression AST for limcalc.
-- All calculus operations derive from limit evaluation over this tree.
-- Const is a placeholder for AlgNum pending full algebraic number implementation.
data Expr
  = Const Double        -- ^ Numeric constant (placeholder for AlgNum)
  | Pi                  -- ^ The constant π
  | E                   -- ^ The constant e
  | I                   -- ^ The imaginary unit
  | Var String          -- ^ Variable by name
  | Add Expr Expr       -- ^ Addition
  | Sub Expr Expr       -- ^ Subtraction
  | Mul Expr Expr       -- ^ Multiplication
  | Div Expr Expr       -- ^ Division
  | Pow Expr Expr       -- ^ Power (base, exponent)
  | Neg Expr            -- ^ Negation
  | Abs Expr            -- ^ Absolute value (primitive: non-analytic over ℂ)
  | Exp Expr            -- ^ e^x
  | Log Expr            -- ^ Natural logarithm
  | Sin Expr            -- ^ Sine
  | Cos Expr            -- ^ Cosine
  deriving (Show, Eq)

-- | Tangent derived from Sin and Cos.
-- Not primitive since Sin/Cos captures all analytic behavior.
tan' :: Expr -> Expr
tan' x = Div (Sin x) (Cos x)