module LimCalc.Simplify where

import LimCalc.Expr

-- | Simplify an expression by applying algebraic identities recursively
-- until no further simplification is possible.
simplify :: Expr -> Expr
simplify expr =
  let expr' = simplifyOnce expr
  in if expr' == expr
     then expr
     else simplify expr'

-- | One pass of simplification — recurse into subexpressions first,
-- then apply identities at the top level.
simplifyOnce :: Expr -> Expr
simplifyOnce (Add f g)  = simplifyAdd (simplifyOnce f) (simplifyOnce g)
simplifyOnce (Sub f g)  = simplifySub (simplifyOnce f) (simplifyOnce g)
simplifyOnce (Mul f g)  = simplifyMul (simplifyOnce f) (simplifyOnce g)
simplifyOnce (Div f g)  = simplifyDiv (simplifyOnce f) (simplifyOnce g)
simplifyOnce (Pow f g)  = simplifyPow (simplifyOnce f) (simplifyOnce g)
simplifyOnce (Neg f)    = simplifyNeg (simplifyOnce f)
simplifyOnce (Exp f)    = Exp (simplifyOnce f)
simplifyOnce (Log f)    = Log (simplifyOnce f)
simplifyOnce (Sin f)    = Sin (simplifyOnce f)
simplifyOnce (Cos f)    = Cos (simplifyOnce f)
simplifyOnce (Abs f)    = Abs (simplifyOnce f)
simplifyOnce (Erf f)    = Erf (simplifyOnce f)
simplifyOnce (Li f)     = Li (simplifyOnce f)
simplifyOnce (Si f)     = Si (simplifyOnce f)
simplifyOnce (Ci f)     = Ci (simplifyOnce f)
simplifyOnce (Ei f)     = Ei (simplifyOnce f)
simplifyOnce e          = e

-- | Simplify addition
simplifyAdd :: Expr -> Expr -> Expr
simplifyAdd (Const 0) y          = y
simplifyAdd x (Const 0)          = x
simplifyAdd (Const a) (Const b)  = Const (a + b)
simplifyAdd x y                  = Add x y

-- | Simplify subtraction
simplifySub :: Expr -> Expr -> Expr
simplifySub x (Const 0)          = x
simplifySub (Const a) (Const b)  = Const (a - b)
simplifySub x y | x == y        = Const 0
simplifySub x y                  = Sub x y

-- | Simplify multiplication
simplifyMul :: Expr -> Expr -> Expr
-- Zero and one
simplifyMul (Const 0) _          = Const 0
simplifyMul _ (Const 0)          = Const 0
simplifyMul (Const 1) y          = y
simplifyMul x (Const 1)          = x
simplifyMul (Const (-1)) y       = Neg y
simplifyMul x (Const (-1))       = Neg x
-- Constant folding
simplifyMul (Const a) (Const b)  = Const (a * b)
-- Double negation
simplifyMul (Neg x) (Neg y)      = simplifyMul x y
-- x^n * (1/x) = x^(n-1)
simplifyMul (Pow x (Const n)) (Div (Const 1) y)
  | x == y                       = simplifyPow x (Const (n - 1))
simplifyMul (Div (Const 1) x) (Pow y (Const n))
  | x == y                       = simplifyPow y (Const (n - 1))
-- x * (1/x) = 1
simplifyMul x (Div (Const 1) y)
  | x == y                       = Const 1
simplifyMul (Div (Const 1) x) y
  | x == y                       = Const 1
-- Pull constant left: x * (c * z) = c * (x * z), but only when x
-- is not itself a constant (otherwise ping-pongs with the rule below)
simplifyMul x (Mul (Const c) z)
  | notConst x                   =
      simplifyMul (Const c) (simplifyMul x z)
-- Pull constant left: (c * y) * z = c * (y * z), but only when z
-- is not itself a constant (otherwise ping-pongs with the rule above)
simplifyMul (Mul (Const c) y) z
  | notConst z                   =
      simplifyMul (Const c) (simplifyMul y z)
-- Flatten: x * (y * z)
simplifyMul x (Mul y z) =
  let xy = simplifyMul x y
  in case xy of
       Mul _ _ ->
         let yz = simplifyMul y z
         in case yz of
              Mul _ _ -> Mul x (Mul y z)
              _       -> simplifyMul x yz
       _ -> simplifyMul xy z
-- Flatten: (x * y) * z
simplifyMul (Mul x y) z =
  let yz = simplifyMul y z
  in case yz of
       Mul _ _ ->
         let xy = simplifyMul x y
         in case xy of
              Mul _ _ -> Mul (Mul x y) z
              _       -> simplifyMul xy z
       _ -> simplifyMul x yz
-- Catch-all
simplifyMul x y                  = Mul x y

-- | Simplify division
simplifyDiv :: Expr -> Expr -> Expr
simplifyDiv (Const 0) _          = Const 0
simplifyDiv x (Const 1)          = x
simplifyDiv (Const a) (Const b)
  | b /= 0                       = Const (a / b)
simplifyDiv x y | x == y        = Const 1
-- Cancel double negation: (-a)/(-b) = a/b
simplifyDiv (Neg a) (Neg b)      = simplifyDiv a b
-- Common constant factor: (c*a)/(c*b) -> a/b
simplifyDiv (Mul (Const c1) a) (Mul (Const c2) b)
  | c1 == c2 && c1 /= 0         = simplifyDiv a b
-- (c*a)/c -> a
simplifyDiv (Mul (Const c1) a) (Const c2)
  | c1 == c2 && c1 /= 0         = a
-- c/(c*b) -> 1/b
simplifyDiv (Const c1) (Mul (Const c2) b)
  | c1 == c2 && c1 /= 0         = simplifyDiv (Const 1) b
simplifyDiv x y                  = Div x y

-- | Simplify power
simplifyPow :: Expr -> Expr -> Expr
simplifyPow _ (Const 0)          = Const 1
simplifyPow x (Const 1)          = x
simplifyPow (Const 0) _          = Const 0
simplifyPow (Const 1) _          = Const 1
simplifyPow (Const a) (Const b)  = Const (a ** b)
simplifyPow x y                  = Pow x y

-- | Simplify negation
simplifyNeg :: Expr -> Expr
simplifyNeg (Const a)            = Const (-a)
simplifyNeg (Neg x)              = x
-- Push negation inside division: -(a/b) = (-a)/b
-- Allows constant folding, e.g. -((-1)/(2x)) -> 1/(2x)
simplifyNeg (Div a b)            = simplifyDiv (simplifyNeg a) b
-- Pull negation into constant factor: -(c*x) = (-c)*x
simplifyNeg (Mul (Const c) x)   = simplifyMul (Const (-c)) x
simplifyNeg x                    = Neg x

-- | True if the expression is NOT a numeric constant.
-- Used to guard the "pull constant left" rules in simplifyMul to
-- prevent them ping-ponging with each other.
notConst :: Expr -> Bool
notConst (Const _) = False
notConst _         = True