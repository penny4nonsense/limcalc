-- | Symbolic Puiseux series expansion for differentiation.
--
-- 'symExpand' expands @f(x + h)@ as a 'SymPuiseuxSeries' in @h@,
-- keeping the base variable @x@ symbolic in the coefficients. The
-- @h^1@ coefficient of the result is the symbolic derivative @f'(x)@,
-- which 'LimCalc.Calculus.diff' extracts via 'symCoeffAt'.
--
-- = Relationship to the numeric expansion path
--
-- 'LimCalc.Expand.expand' performs the same conceptual operation but
-- evaluates coefficients numerically as 'AlgNum' values. 'symExpand'
-- keeps coefficients symbolic so that the result is an 'Expr' rather
-- than a number. The Taylor series generators here (@sinSymTaylor@,
-- @expSymTaylor@, etc.) mirror their numeric counterparts in
-- 'LimCalc.Expand', but produce 'Expr' coefficients computed via
-- 'LimCalc.DiffField.deriveBase' rather than numeric evaluation.
--
-- = Limitations
--
-- Unlike the numeric path, 'symExpand' handles only pure power series
-- (no @log(h)@ terms). Functions with logarithmic singularities at
-- generic points ('Li', 'Abs') return 'Unknown'. This is acceptable
-- for differentiation, since the derivative of a smooth function at
-- a generic point is always a pure power series.
module LimCalc.SymExpand
  ( -- * Symbolic expansion
    symExpand
    -- * Taylor series generators
  , sinSymTaylor
  , cosSymTaylor
  , expSymTaylor
  , logSymTaylor
  , erfSymTaylor
  , siSymTaylor
  , ciSymTaylor
  , eiSymTaylor
    -- * Series operations
  , symEvalSeriesAt
  , symPowSeries
  , symInvertSeries
  , symGeometricSeries
  , symNormalizeW
  , symExpandPowR
  , symBinomialSeries
  , symBinomCoeffs
    -- * Utilities
  , symConstantTerm
  , depth
  ) where

import LimCalc.Expr
import LimCalc.SymPuiseux
import LimCalc.Types
import LimCalc.DiffField (deriveBase)
import LimCalc.Simplify
import Data.Ratio (numerator, denominator)

-- | Symbolically expand @f(x + h)@ as a 'SymPuiseuxSeries' in @h@.
--
-- The variable @var@ is the expansion variable; all other variables
-- remain as 'Var' nodes in the coefficients. The result is a series
-- @∑ cₙ(x) · h^n@ where each @cₙ(x)@ is an 'Expr'.
--
-- The @h^1@ coefficient is the symbolic derivative of @f@ with
-- respect to @var@, extracted by 'LimCalc.Calculus.diff' via
-- 'symCoeffAt'.
symExpand :: Expr -> String -> Either ExpandError SymPuiseuxSeries

symExpand (Const c) _ = Right $ SymPuiseuxSeries [SymPuiseuxTerm 0 (Const c)]
symExpand Pi        _ = Right $ SymPuiseuxSeries [SymPuiseuxTerm 0 Pi]
symExpand E         _ = Right $ SymPuiseuxSeries [SymPuiseuxTerm 0 E]
symExpand I         _ = Right $ SymPuiseuxSeries [SymPuiseuxTerm 0 I]

symExpand (Var name) var
  | name == var =
      Right $ SymPuiseuxSeries
        [ SymPuiseuxTerm 0 (Var name)
        , SymPuiseuxTerm 1 (Const 1)
        ]
  | otherwise =
      Right $ SymPuiseuxSeries [SymPuiseuxTerm 0 (Var name)]

symExpand (Add f g) var = do
  sf <- symExpand f var
  sg <- symExpand g var
  return $ addSymSeries sf sg

symExpand (Sub f g) var = do
  sf <- symExpand f var
  sg <- symExpand g var
  return $ addSymSeries sf (scaleSymSeries (Const (-1)) sg)

symExpand (Mul f g) var = do
  sf <- symExpand f var
  sg <- symExpand g var
  return $ mulSymSeries sf sg

symExpand (Neg f) var = do
  sf <- symExpand f var
  return $ scaleSymSeries (Const (-1)) sf

-- | @sin(f(x+h))@ via symbolic Taylor series around @f(x)@.
symExpand (Sin f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = sinSymTaylor c0
  return $ symEvalSeriesAt u s

-- | @cos(f(x+h))@ via symbolic Taylor series around @f(x)@.
symExpand (Cos f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = cosSymTaylor c0
  return $ symEvalSeriesAt u s

-- | @exp(f(x+h))@ via symbolic Taylor series around @f(x)@.
symExpand (Exp f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = expSymTaylor c0
  return $ symEvalSeriesAt u s

-- | @log(f(x+h))@ via symbolic Taylor series around @f(x)@.
symExpand (Log f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = logSymTaylor c0
  return $ symEvalSeriesAt u s

symExpand (Div f g) var = do
  sf <- symExpand f var
  sg <- symExpand g var
  return $ mulSymSeries sf (symInvertSeries sg)

symExpand (Pow f (Const r)) var     = symExpandPowR f (toRational r) var
symExpand (Pow f (Neg (Const r))) var = symExpandPowR f (toRational (-r)) var
symExpand (Pow _ _) _               = Left $ Unknown "Symbolic exponents not yet supported"

symExpand (Abs _)    _ = Left $ Unknown "Abs not yet implemented"
symExpand (Arcsin _) _ = Left $ Unknown "Arcsin expansion not yet implemented"
symExpand (Arccos _) _ = Left $ Unknown "Arccos expansion not yet implemented"
symExpand (Arctan _) _ = Left $ Unknown "Arctan expansion not yet implemented"

-- | @erf(f(x+h))@ via symbolic Taylor series. Entire function.
symExpand (Erf f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = erfSymTaylor c0
  return $ symEvalSeriesAt u s

-- | @Si(f(x+h))@ via symbolic Taylor series. Entire function.
symExpand (Si f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = siSymTaylor c0
  return $ symEvalSeriesAt u s

-- | @Ei(f(x+h))@ via symbolic Taylor series. Analytic away from @x=0@.
symExpand (Ei f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = eiSymTaylor c0
  return $ symEvalSeriesAt u s

-- | @Ci(f(x+h))@ via symbolic Taylor series. Analytic away from @x=0@.
symExpand (Ci f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = ciSymTaylor c0
  return $ symEvalSeriesAt u s

-- | @Li@ requires @log(log(x))@; not supported by the symbolic path.
symExpand (Li _) _ = Left $ Unknown
  "Li expansion not yet implemented: li(x) = Ei(log(x)), \
  \requires doubly-logarithmic Puiseux extension"

-- | Number of Taylor series terms to compute.
depth :: Int
depth = 6

-- | Extract the constant term (@h^0@ coefficient) of a symbolic series.
-- Returns @Const 0@ if the series is empty or has no constant term.
symConstantTerm :: SymPuiseuxSeries -> Expr
symConstantTerm (SymPuiseuxSeries []) = Const 0
symConstantTerm (SymPuiseuxSeries (t:_))
  | symExp t == 0 = symCoeff t
  | otherwise     = Const 0

-- | Evaluate @S = ∑ aₙ · tⁿ@ by substituting @t = u@ symbolically.
symEvalSeriesAt :: SymPuiseuxSeries -> SymPuiseuxSeries -> SymPuiseuxSeries
symEvalSeriesAt u (SymPuiseuxSeries ts) =
  foldr addSymSeries zeroSym
    [ scaleSymSeries (symCoeff t) (symPowSeries u (symExp t))
    | t <- ts
    ]

-- | Raise a symbolic series to a rational power.
-- Only handles non-negative integer powers; returns the empty series otherwise.
symPowSeries :: SymPuiseuxSeries -> Rational -> SymPuiseuxSeries
symPowSeries _ 0 = SymPuiseuxSeries [SymPuiseuxTerm 0 (Const 1)]
symPowSeries u n
  | n > 0 && denominator n == 1 =
      let k = fromIntegral (numerator n) :: Int
      in foldl mulSymSeries (SymPuiseuxSeries [SymPuiseuxTerm 0 (Const 1)]) (replicate k u)
  | otherwise = SymPuiseuxSeries []

-- | Symbolic Taylor series for @sin@ around symbolic @x@.
--
-- @sin(x + h) = sin(x) + cos(x)·h − sin(x)·h²\/2 − cos(x)·h³\/6 + …@
sinSymTaylor :: Expr -> SymPuiseuxSeries
sinSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (sinSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    sinSymCoeff n =
      let bases  = cycle [Sin x, Cos x, Neg (Sin x), Neg (Cos x)]
          facts  = scanl (*) 1 [1..] :: [Int]
          base   = bases !! n
          factor = fromIntegral (facts !! n)
      in if factor == 1 then base else Div base (Const factor)

-- | Symbolic Taylor series for @cos@ around symbolic @x@.
--
-- @cos(x + h) = cos(x) − sin(x)·h − cos(x)·h²\/2 + sin(x)·h³\/6 + …@
cosSymTaylor :: Expr -> SymPuiseuxSeries
cosSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (cosSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    cosSymCoeff n =
      let bases  = cycle [Cos x, Neg (Sin x), Neg (Cos x), Sin x]
          facts  = scanl (*) 1 [1..] :: [Int]
          base   = bases !! n
          factor = fromIntegral (facts !! n)
      in if factor == 1 then base else Div base (Const factor)

-- | Symbolic Taylor series for @exp@ around symbolic @x@.
--
-- @exp(x + h) = exp(x) + exp(x)·h + exp(x)·h²\/2 + …@
expSymTaylor :: Expr -> SymPuiseuxSeries
expSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (expSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    expSymCoeff n =
      let facts  = scanl (*) 1 [1..] :: [Int]
          factor = fromIntegral (facts !! n)
      in if factor == 1 then Exp x else Div (Exp x) (Const factor)

-- | Symbolic Taylor series for @log@ around symbolic @x@.
--
-- @log(x + h) = log(x) + h\/x − h²\/(2x²) + h³\/(3x³) − …@
logSymTaylor :: Expr -> SymPuiseuxSeries
logSymTaylor x = SymPuiseuxSeries $ take depth $
  SymPuiseuxTerm 0 (Log x) :
  [ SymPuiseuxTerm (fromIntegral n) (logSymCoeff n)
  | n <- [1..] :: [Int] ]
  where
    logSymCoeff n =
      let sign   = if even n then Const (-1) else Const 1
          factor = Const (fromIntegral n)
          xpow   = Pow x (Const (fromIntegral n))
      in Div sign (Mul factor xpow)

-- | Symbolic Taylor series for @Ei@ around symbolic @x@.
--
-- Coefficients are computed by iterating 'deriveBase' on @Ei(x)@.
-- @Ei'(x) = e^x\/x@; subsequent derivatives follow by the chain rule.
eiSymTaylor :: Expr -> SymPuiseuxSeries
eiSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (eiSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    derivs = iterate (\e -> simplify (deriveBase e)) (Ei x)
    eiSymCoeff 0 = Ei x
    eiSymCoeff n =
      let d      = derivs !! n
          facts  = scanl (*) 1 [1..] :: [Int]
          factor = fromIntegral (facts !! n)
      in if factor == 1 then d else Div d (Const factor)

-- | Symbolic Taylor series for @Ci@ around symbolic @x@.
--
-- Coefficients are computed by iterating 'deriveBase' on @Ci(x)@.
-- @Ci'(x) = cos(x)\/x@.
ciSymTaylor :: Expr -> SymPuiseuxSeries
ciSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (ciSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    derivs = iterate (\e -> simplify (deriveBase e)) (Ci x)
    ciSymCoeff 0 = Ci x
    ciSymCoeff n =
      let d      = derivs !! n
          facts  = scanl (*) 1 [1..] :: [Int]
          factor = fromIntegral (facts !! n)
      in if factor == 1 then d else Div d (Const factor)

-- | Symbolic Taylor series for @erf@ around symbolic @x@.
--
-- Coefficients are computed by iterating 'deriveBase' on @erf(x)@.
-- @erf'(x) = (2\/√π) · e^(−x²)@.
erfSymTaylor :: Expr -> SymPuiseuxSeries
erfSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (erfSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    derivs = iterate (\e -> simplify (deriveBase e)) (Erf x)
    erfSymCoeff 0 = Erf x
    erfSymCoeff n =
      let d      = derivs !! n
          facts  = scanl (*) 1 [1..] :: [Int]
          factor = fromIntegral (facts !! n)
      in if factor == 1 then d else Div d (Const factor)

-- | Symbolic Taylor series for @Si@ around symbolic @x@.
--
-- Coefficients are computed by iterating 'deriveBase' on @Si(x)@.
-- @Si'(x) = sin(x)\/x@.
siSymTaylor :: Expr -> SymPuiseuxSeries
siSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (siSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    derivs = iterate (\e -> simplify (deriveBase e)) (Si x)
    siSymCoeff 0 = Si x
    siSymCoeff n =
      let d      = derivs !! n
          facts  = scanl (*) 1 [1..] :: [Int]
          factor = fromIntegral (facts !! n)
      in if factor == 1 then d else Div d (Const factor)

-- | Invert a symbolic series: @1\/s@ via geometric series expansion.
--
-- Factors out the leading term @a · h^alpha@, then computes
-- @(1\/a) · h^(−alpha) · (1\/(1 + w))@ where @w = s\/(a·h^alpha) − 1@
-- has no constant term.
symInvertSeries :: SymPuiseuxSeries -> SymPuiseuxSeries
symInvertSeries s =
  case symLeadingTerm s of
    Nothing -> SymPuiseuxSeries []
    Just lt ->
      let alpha = symExp lt
          a     = symCoeff lt
          w     = symNormalizeW s lt alpha a
          negw  = scaleSymSeries (Const (-1)) w
          geo   = symGeometricSeries negw
          shift = negate alpha
      in shiftSymExponents shift
           (scaleSymSeries (Div (Const 1) a) geo)

-- | Geometric series @1\/(1 + w) = ∑ (−w)^n@, truncated to 'depth' terms.
symGeometricSeries :: SymPuiseuxSeries -> SymPuiseuxSeries
symGeometricSeries u =
  let upows = take (depth+1) $
                iterate (truncateSymSeries depth . mulSymSeries u)
                        (SymPuiseuxSeries [SymPuiseuxTerm 0 (Const 1)])
  in truncateSymSeries depth $ foldr addSymSeries zeroSym upows

-- | Compute @w = s\/(a · h^alpha) − 1@: normalise a series by
-- factoring out its leading term and subtracting 1.
symNormalizeW :: SymPuiseuxSeries -> SymPuiseuxTerm -> Rational -> Expr
              -> SymPuiseuxSeries
symNormalizeW (SymPuiseuxSeries ts) _lt alpha a =
  let shifted = [ SymPuiseuxTerm (symExp t - alpha) (Div (symCoeff t) a) | t <- ts ]
  in removeSymTerm 0 (SymPuiseuxSeries shifted)

-- | Symbolic power expansion @f^r@ for rational @r@.
symExpandPowR :: Expr -> Rational -> String -> Either ExpandError SymPuiseuxSeries
symExpandPowR f r var = do
  s <- symExpand f var
  let s' = truncateSymSeries depth s
  case symLeadingTerm s' of
    Nothing -> Right zeroSym
    Just lt ->
      let alpha = symExp lt
          a     = symCoeff lt
          w     = truncateSymSeries depth (symNormalizeW s' lt alpha a)
          binom = symBinomialSeries r w
          shift = alpha * r
      in Right $ truncateSymSeries depth
           (shiftSymExponents shift
             (scaleSymSeries (Pow a (Const (fromRational r))) binom))

-- | Symbolic binomial series @(1 + w)^r@, truncated to 'depth' terms.
symBinomialSeries :: Rational -> SymPuiseuxSeries -> SymPuiseuxSeries
symBinomialSeries r w =
  let bcs   = symBinomCoeffs r depth
      wpows = take (depth+1) $
                iterate (truncateSymSeries depth . mulSymSeries w)
                        (SymPuiseuxSeries [SymPuiseuxTerm 0 (Const 1)])
  in truncateSymSeries depth $ foldr addSymSeries zeroSym
       [ scaleSymSeries c wp | (c, wp) <- zip bcs wpows ]

-- | Generalised binomial coefficients as 'Expr' values:
-- @C(r, 0) = 1@, @C(r, k) = r(r−1)···(r−k+1) \/ k!@.
symBinomCoeffs :: Rational -> Int -> [Expr]
symBinomCoeffs r n = take (n+1) $ scanl step (Const 1) [0..n-1]
  where
    step acc k =
      let rk = fromRational r - fromIntegral k
          kk = fromIntegral (k+1)
      in Div (Mul acc (Const rk)) (Const kk)