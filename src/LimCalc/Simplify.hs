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
simplifyOnce (Arcsin f) = Arcsin (simplifyOnce f)
simplifyOnce (Arccos f) = Arccos (simplifyOnce f)
simplifyOnce (Arctan f) = Arctan (simplifyOnce f)
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

------------------------------------------------------------------------
-- Euler folding: convert e^(i*theta) expressions back to sin/cos
------------------------------------------------------------------------

-- | Fold Euler-form expressions into sin/cos wherever possible.
-- Called explicitly after Risch integration to produce human-readable
-- output, rather than as part of the simplifyOnce fixed-point loop
-- (which would make every simplification aware of trig forms).
--
-- Recognizes two patterns in Add nodes:
--   c * e^(i*theta) + c * e^(i*theta)^(-1) = 2c * cos(theta)
--   (-c*i) * e^(i*theta) + (c*i) * e^(i*theta)^(-1) = 2c * sin(theta)
foldEuler :: Expr -> Expr
foldEuler expr = simplify (foldEulerOnce expr)

foldEulerOnce :: Expr -> Expr
foldEulerOnce (Add a b) =
  case tryFoldEulerPair (foldEulerOnce a) (foldEulerOnce b) of
    Just folded -> folded
    Nothing     -> Add (foldEulerOnce a) (foldEulerOnce b)
foldEulerOnce (Sub a b) = Sub (foldEulerOnce a) (foldEulerOnce b)
foldEulerOnce (Mul a b) = Mul (foldEulerOnce a) (foldEulerOnce b)
foldEulerOnce (Div a b) = Div (foldEulerOnce a) (foldEulerOnce b)
foldEulerOnce (Neg a)   = Neg (foldEulerOnce a)
foldEulerOnce (Pow a b) = Pow (foldEulerOnce a) (foldEulerOnce b)
foldEulerOnce e         = e

-- | Try to recognize a pair of terms as an Euler identity.
-- Returns Just (sin/cos form) if recognized, Nothing otherwise.
tryFoldEulerPair :: Expr -> Expr -> Maybe Expr
tryFoldEulerPair a b = do
  (ca, ta, na) <- extractEulerTerm a
  (cb, tb, nb) <- extractEulerTerm b
  -- Terms must share the same base t = e^(i*theta)
  if ta /= tb then Nothing
  else if na == 1 && nb == -1
    -- ca*t + cb*t^(-1): check for cos pattern (ca == cb, both real)
    -- or sin pattern (ca == -cb, both imaginary)
    then tryFoldCosOrSin ca cb ta
  else if na == -1 && nb == 1
    then tryFoldCosOrSin cb ca ta
  else Nothing

-- | Extract (coefficient, base-exp-arg, power) from an Euler term.
-- Recognizes:
--   c * Exp(f)       -> (c, f, 1)
--   c * Exp(f)^n     -> (c, f, n)
--   c * I * Exp(f)   -> (c*i, f, 1)   where i = sqrt(-1) ~ (0,1) complex
-- Returns Nothing for anything not in this form.
extractEulerTerm :: Expr -> Maybe (Expr, Expr, Double)
-- Mul (Const c) (Exp f) -> (Const c, f, 1)
extractEulerTerm (Mul (Const c) (Exp f))         = Just (Const c, f, 1)
-- Mul (Const c) (Pow (Exp f) (Const n)) -> (Const c, f, n)
extractEulerTerm (Mul (Const c) (Pow (Exp f) (Const n))) = Just (Const c, f, n)
-- Mul (Const c) (Mul I (Exp f)) -> (Const c * I, f, 1)
extractEulerTerm (Mul (Const c) (Mul I (Exp f))) = Just (Mul (Const c) I, f, 1)
-- Mul (Const c) (Mul I (Pow (Exp f) (Const n)))
extractEulerTerm (Mul (Const c) (Mul I (Pow (Exp f) (Const n)))) =
  Just (Mul (Const c) I, f, n)
extractEulerTerm _ = Nothing

-- | Given coefficients ca, cb for terms t and t^(-1), try to fold
-- into cos or sin using Euler's formula.
tryFoldCosOrSin :: Expr -> Expr -> Expr -> Maybe Expr
tryFoldCosOrSin ca cb base
  -- Cosine pattern: ca == cb (same real coefficient)
  -- ca*e^(i*theta) + ca*e^(-i*theta) = 2ca*cos(theta)
  | ca == cb, Just (Const c) <- pure ca, Just theta <- extractITheta base =
      Just (simplify (Mul (Const (2 * c)) (Cos theta)))
  -- Sine pattern: ca = -c*I, cb = c*I (imaginary, opposite signs)
  -- (-c*i)*e^(i*theta) + (c*i)*e^(-i*theta) = 2c*sin(theta)
  | Just (Const c) <- extractRealFromIMul cb
  , Just (Const c') <- extractRealFromIMul ca
  , c == -c'
  , Just theta <- extractITheta base =
      Just (simplify (Mul (Const (2 * c)) (Sin theta)))
  | otherwise = Nothing

-- | Extract theta from I*theta (the argument to e^(i*theta)).
-- Returns Just theta if the expression is Mul I theta, Nothing otherwise.
extractITheta :: Expr -> Maybe Expr
extractITheta (Mul I theta)     = Just theta
extractITheta (Mul theta I)     = Just theta
extractITheta _                 = Nothing

-- | Extract the real scalar c from c*I expressions.
extractRealFromIMul :: Expr -> Maybe Expr
extractRealFromIMul (Mul (Const c) I) = Just (Const c)
extractRealFromIMul (Mul I (Const c)) = Just (Const c)
extractRealFromIMul _                 = Nothing

------------------------------------------------------------------------
-- Log-law simplification
------------------------------------------------------------------------

-- | Apply log laws to simplify expressions involving logarithms.
-- Called explicitly after Risch integration.
--   c*log(a) + d*log(b) -> various simplifications
--   c*log(a) -> log(a^c)  when c is rational
foldLogs :: Expr -> Expr
foldLogs = simplify . foldLogsOnce

foldLogsOnce :: Expr -> Expr
foldLogsOnce (Add a b) =
  case tryFoldLogPair (foldLogsOnce a) (foldLogsOnce b) of
    Just folded -> folded
    Nothing     -> Add (foldLogsOnce a) (foldLogsOnce b)
foldLogsOnce (Sub a b) =
  case tryFoldLogPair (foldLogsOnce a) (Neg (foldLogsOnce b)) of
    Just folded -> folded
    Nothing     -> Sub (foldLogsOnce a) (foldLogsOnce b)
foldLogsOnce (Mul a b) = Mul (foldLogsOnce a) (foldLogsOnce b)
foldLogsOnce (Div a b) = Div (foldLogsOnce a) (foldLogsOnce b)
foldLogsOnce (Neg a)   = Neg (foldLogsOnce a)
foldLogsOnce e         = e

-- | Extract (coefficient, log-argument) from a log term.
-- Recognizes: Log a -> (1, a), Mul (Const c) (Log a) -> (c, a),
-- Neg (Log a) -> (-1, a), Neg (Mul (Const c) (Log a)) -> (-c, a)
extractLogTerm :: Expr -> Maybe (Double, Expr)
extractLogTerm (Log a)                    = Just (1, a)
extractLogTerm (Mul (Const c) (Log a))   = Just (c, a)
extractLogTerm (Neg (Log a))             = Just (-1, a)
extractLogTerm (Neg (Mul (Const c) (Log a))) = Just (-c, a)
extractLogTerm _                          = Nothing

-- | Try to fold a pair of log terms using log laws.
tryFoldLogPair :: Expr -> Expr -> Maybe Expr
tryFoldLogPair a b = do
  (ca, la) <- extractLogTerm a
  (cb, lb) <- extractLogTerm b
  if ca == cb
    -- c*log(a) + c*log(b) = c*log(a*b)
    then Just (simplify (Mul (Const ca) (Log (Mul la lb))))
  else if ca == -cb
    -- c*log(a) - c*log(b) = c*log(a/b)
    then Just (simplify (Mul (Const ca) (Log (Div la lb))))
  else Nothing