-- | Algebraic simplification of 'Expr' trees.
--
-- Provides two distinct post-processing operations:
--
-- * 'simplify': a fixed-point algebraic simplifier applying standard
--   identities (zero\/one rules, constant folding, double negation,
--   etc.). Used throughout the codebase to normalise symbolic output.
--
-- * 'foldEuler': converts Euler-form expressions (@c·e^(iθ) + c·e^(−iθ)@)
--   to @sin\/cos@ form. Called explicitly after Risch integration to
--   produce human-readable trig output.
--
-- * 'foldLogs': applies log laws to combine @c·log(a) + c·log(b)@
--   into @c·log(a·b)@ (and similarly for differences). Called
--   explicitly after Risch integration.
--
-- 'foldEuler' and 'foldLogs' are /not/ part of the 'simplifyOnce'
-- fixed-point loop. Wiring them into 'simplify' would make every
-- simplification trig- and log-aware, introducing unwanted
-- interactions. They are instead called explicitly at integration
-- result boundaries in 'LimCalc.Risch'.
module LimCalc.Simplify
  ( -- * Algebraic simplification
    simplify
  , simplifyOnce
    -- * Operator-level simplifiers
  , simplifyAdd
  , simplifySub
  , simplifyMul
  , simplifyDiv
  , simplifyPow
  , simplifyNeg
    -- * Euler folding
  , foldEuler
  , foldEulerOnce
  , tryFoldEulerPair
  , extractEulerTerm
  , tryFoldCosOrSin
  , extractITheta
  , extractRealFromIMul
    -- * Log folding
  , foldLogs
  , foldLogsOnce
  , extractLogTerm
  , tryFoldLogPair
    -- * Utilities
  , notConst
  ) where

import LimCalc.Expr

-- | Simplify an 'Expr' by repeatedly applying algebraic identities
-- until a fixed point is reached.
--
-- The fixed-point loop calls 'simplifyOnce' until the expression no
-- longer changes. Termination is guaranteed because each successful
-- rule strictly reduces the expression (no rule increases the
-- number of constructors).
simplify :: Expr -> Expr
simplify expr =
  let expr' = simplifyOnce expr
  in if expr' == expr
     then expr
     else simplify expr'

-- | One pass of simplification: recurse into subexpressions first,
-- then apply identities at the top level.
simplifyOnce :: Expr -> Expr
simplifyOnce (Add f g)  = simplifyAdd  (simplifyOnce f) (simplifyOnce g)
simplifyOnce (Sub f g)  = simplifySub  (simplifyOnce f) (simplifyOnce g)
simplifyOnce (Mul f g)  = simplifyMul  (simplifyOnce f) (simplifyOnce g)
simplifyOnce (Div f g)  = simplifyDiv  (simplifyOnce f) (simplifyOnce g)
simplifyOnce (Pow f g)  = simplifyPow  (simplifyOnce f) (simplifyOnce g)
simplifyOnce (Neg f)    = simplifyNeg  (simplifyOnce f)
simplifyOnce (Exp f)    = Exp    (simplifyOnce f)
simplifyOnce (Log f)    = Log    (simplifyOnce f)
simplifyOnce (Sin f)    = Sin    (simplifyOnce f)
simplifyOnce (Cos f)    = Cos    (simplifyOnce f)
simplifyOnce (Arcsin f) = Arcsin (simplifyOnce f)
simplifyOnce (Arccos f) = Arccos (simplifyOnce f)
simplifyOnce (Arctan f) = Arctan (simplifyOnce f)
simplifyOnce (Abs f)    = Abs    (simplifyOnce f)
simplifyOnce (Erf f)    = Erf    (simplifyOnce f)
simplifyOnce (Li f)     = Li     (simplifyOnce f)
simplifyOnce (Si f)     = Si     (simplifyOnce f)
simplifyOnce (Ci f)     = Ci     (simplifyOnce f)
simplifyOnce (Ei f)     = Ei     (simplifyOnce f)
simplifyOnce e          = e

-- | Simplify an addition node.
simplifyAdd :: Expr -> Expr -> Expr
simplifyAdd (Const 0) y         = y
simplifyAdd x (Const 0)         = x
simplifyAdd (Const a) (Const b) = Const (a + b)
simplifyAdd x y                 = Add x y

-- | Simplify a subtraction node.
simplifySub :: Expr -> Expr -> Expr
simplifySub x (Const 0)         = x
simplifySub (Const a) (Const b) = Const (a - b)
simplifySub x y | x == y       = Const 0
simplifySub x y                 = Sub x y

-- | Simplify a multiplication node.
--
-- Handles: zero\/one annihilation, constant folding, double negation,
-- power-times-reciprocal reduction (@x^n · (1\/x) = x^(n−1)@),
-- and constant-left normalisation (to prevent ping-pong between rules).
simplifyMul :: Expr -> Expr -> Expr
simplifyMul (Const 0) _                           = Const 0
simplifyMul _ (Const 0)                           = Const 0
simplifyMul (Const 1) y                           = y
simplifyMul x (Const 1)                           = x
simplifyMul (Const (-1)) y                        = Neg y
simplifyMul x (Const (-1))                        = Neg x
simplifyMul (Const a) (Const b)                   = Const (a * b)
simplifyMul (Neg x) (Neg y)                       = simplifyMul x y
simplifyMul (Pow x (Const n)) (Div (Const 1) y)
  | x == y                                        = simplifyPow x (Const (n - 1))
simplifyMul (Div (Const 1) x) (Pow y (Const n))
  | x == y                                        = simplifyPow y (Const (n - 1))
simplifyMul x (Div (Const 1) y)
  | x == y                                        = Const 1
simplifyMul (Div (Const 1) x) y
  | x == y                                        = Const 1
simplifyMul x (Mul (Const c) z)
  | notConst x                                    =
      simplifyMul (Const c) (simplifyMul x z)
simplifyMul (Mul (Const c) y) z
  | notConst z                                    =
      simplifyMul (Const c) (simplifyMul y z)
simplifyMul x (Mul y z) =
  let xy = simplifyMul x y
  in case xy of
       Mul _ _ ->
         let yz = simplifyMul y z
         in case yz of
              Mul _ _ -> Mul x (Mul y z)
              _       -> simplifyMul x yz
       _ -> simplifyMul xy z
simplifyMul (Mul x y) z =
  let yz = simplifyMul y z
  in case yz of
       Mul _ _ ->
         let xy = simplifyMul x y
         in case xy of
              Mul _ _ -> Mul (Mul x y) z
              _       -> simplifyMul xy z
       _ -> simplifyMul x yz
simplifyMul x y                                   = Mul x y

-- | Simplify a division node.
simplifyDiv :: Expr -> Expr -> Expr
simplifyDiv (Const 0) _                           = Const 0
simplifyDiv x (Const 1)                           = x
simplifyDiv (Const a) (Const b)
  | b /= 0                                        = Const (a / b)
simplifyDiv x y | x == y                          = Const 1
simplifyDiv (Neg a) (Neg b)                       = simplifyDiv a b
simplifyDiv (Mul (Const c1) a) (Mul (Const c2) b)
  | c1 == c2 && c1 /= 0                           = simplifyDiv a b
simplifyDiv (Mul (Const c1) a) (Const c2)
  | c1 == c2 && c1 /= 0                           = a
simplifyDiv (Const c1) (Mul (Const c2) b)
  | c1 == c2 && c1 /= 0                           = simplifyDiv (Const 1) b
simplifyDiv x y                                   = Div x y

-- | Simplify a power node.
simplifyPow :: Expr -> Expr -> Expr
simplifyPow _ (Const 0)         = Const 1
simplifyPow x (Const 1)         = x
simplifyPow (Const 0) _         = Const 0
simplifyPow (Const 1) _         = Const 1
simplifyPow (Const a) (Const b) = Const (a ** b)
simplifyPow x y                 = Pow x y

-- | Simplify a negation node.
--
-- Pushes negation inside division (@−(a\/b) = (−a)\/b@) to enable
-- constant folding (e.g. @−((−1)\/(2x)) → 1\/(2x)@), and extracts
-- constant factors from negated multiplications.
simplifyNeg :: Expr -> Expr
simplifyNeg (Const a)          = Const (-a)
simplifyNeg (Neg x)            = x
simplifyNeg (Div a b)          = simplifyDiv (simplifyNeg a) b
simplifyNeg (Mul (Const c) x)  = simplifyMul (Const (-c)) x
simplifyNeg x                  = Neg x

-- | True if the expression is not a numeric constant.
--
-- Used to guard the "pull constant left" rules in 'simplifyMul',
-- preventing them from ping-ponging against each other.
notConst :: Expr -> Bool
notConst (Const _) = False
notConst _         = True

------------------------------------------------------------------------
-- Euler folding
------------------------------------------------------------------------

-- | Fold Euler-form expressions into @sin\/cos@ wherever possible,
-- then simplify.
--
-- Recognises two patterns in 'Add' nodes:
--
-- * @c · e^(iθ) + c · e^(−iθ) = 2c · cos(θ)@
-- * @(−ci) · e^(iθ) + (ci) · e^(−iθ) = 2c · sin(θ)@
--
-- Called explicitly after Risch trig integration to convert the
-- Euler-substitution output back to @sin\/cos@ form. Not part of
-- the 'simplifyOnce' fixed-point loop to avoid making every
-- simplification trig-aware.
foldEuler :: Expr -> Expr
foldEuler expr = simplify (foldEulerOnce expr)

-- | One pass of Euler folding, without the fixed-point loop.
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

-- | Try to recognise a pair of terms as an Euler identity.
-- Returns the @sin\/cos@ form if recognised, 'Nothing' otherwise.
tryFoldEulerPair :: Expr -> Expr -> Maybe Expr
tryFoldEulerPair a b = do
  (ca, ta, na) <- extractEulerTerm a
  (cb, tb, nb) <- extractEulerTerm b
  if ta /= tb then Nothing
  else if na == 1 && nb == -1
    then tryFoldCosOrSin ca cb ta
  else if na == -1 && nb == 1
    then tryFoldCosOrSin cb ca ta
  else Nothing

-- | Extract @(coefficient, base-exp-arg, power)@ from an Euler term.
--
-- Recognises:
--
-- * @c · exp(f)@       → @(c, f, 1)@
-- * @c · exp(f)^n@     → @(c, f, n)@
-- * @c · I · exp(f)@   → @(c·I, f, 1)@
-- * @c · I · exp(f)^n@ → @(c·I, f, n)@
extractEulerTerm :: Expr -> Maybe (Expr, Expr, Double)
extractEulerTerm (Mul (Const c) (Exp f))                         =
  Just (Const c, f, 1)
extractEulerTerm (Mul (Const c) (Pow (Exp f) (Const n)))         =
  Just (Const c, f, n)
extractEulerTerm (Mul (Const c) (Mul I (Exp f)))                 =
  Just (Mul (Const c) I, f, 1)
extractEulerTerm (Mul (Const c) (Mul I (Pow (Exp f) (Const n)))) =
  Just (Mul (Const c) I, f, n)
extractEulerTerm _                                                =
  Nothing

-- | Given coefficients @ca@ and @cb@ for terms @t@ and @t^(−1)@,
-- try to fold into @cos@ or @sin@ via Euler's formula.
tryFoldCosOrSin :: Expr -> Expr -> Expr -> Maybe Expr
tryFoldCosOrSin ca cb base
  | ca == cb
  , Just (Const c) <- pure ca
  , Just theta     <- extractITheta base =
      Just (simplify (Mul (Const (2 * c)) (Cos theta)))
  | Just (Const c)  <- extractRealFromIMul cb
  , Just (Const c') <- extractRealFromIMul ca
  , c == -c'
  , Just theta      <- extractITheta base =
      Just (simplify (Mul (Const (2 * c)) (Sin theta)))
  | otherwise = Nothing

-- | Extract @θ@ from @I · θ@ (the argument to @e^(iθ)@).
extractITheta :: Expr -> Maybe Expr
extractITheta (Mul I theta) = Just theta
extractITheta (Mul theta I) = Just theta
extractITheta _             = Nothing

-- | Extract the real scalar @c@ from @c · I@ expressions.
extractRealFromIMul :: Expr -> Maybe Expr
extractRealFromIMul (Mul (Const c) I) = Just (Const c)
extractRealFromIMul (Mul I (Const c)) = Just (Const c)
extractRealFromIMul _                 = Nothing

------------------------------------------------------------------------
-- Log folding
------------------------------------------------------------------------

-- | Apply log laws to combine logarithm terms, then simplify.
--
-- Recognises:
--
-- * @c · log(a) + c · log(b) = c · log(a · b)@
-- * @c · log(a) − c · log(b) = c · log(a \/ b)@
--
-- Called explicitly after Risch integration. Not part of the
-- 'simplifyOnce' fixed-point loop.
foldLogs :: Expr -> Expr
foldLogs = simplify . foldLogsOnce

-- | One pass of log folding, without the fixed-point loop.
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

-- | Extract @(coefficient, log-argument)@ from a log term.
--
-- Recognises: @log(a)@, @c · log(a)@, @−log(a)@, @−(c · log(a))@.
extractLogTerm :: Expr -> Maybe (Double, Expr)
extractLogTerm (Log a)                       = Just (1, a)
extractLogTerm (Mul (Const c) (Log a))       = Just (c, a)
extractLogTerm (Neg (Log a))                 = Just (-1, a)
extractLogTerm (Neg (Mul (Const c) (Log a))) = Just (-c, a)
extractLogTerm _                             = Nothing

-- | Try to fold a pair of log terms using log laws.
-- Returns the combined form if the coefficients match (or are
-- negatives of each other), 'Nothing' otherwise.
tryFoldLogPair :: Expr -> Expr -> Maybe Expr
tryFoldLogPair a b = do
  (ca, la) <- extractLogTerm a
  (cb, lb) <- extractLogTerm b
  if ca == cb
    then Just (simplify (Mul (Const ca) (Log (Mul la lb))))
  else if ca == -cb
    then Just (simplify (Mul (Const ca) (Log (Div la lb))))
  else Nothing