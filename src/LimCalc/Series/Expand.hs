-- | Log-Puiseux series expansion engine.
--
-- 'expand' computes the local log-Puiseux series expansion of an
-- 'Expr' around a base point @x₀@, producing a 'LogPuiseuxSeries'
-- in the perturbation variable @h@ such that the series represents
-- @f(x₀ + h)@.
--
-- = Core thesis
--
-- The derivative is a limit, so compute it that way: expand
-- @f(x₀ + h)@ as a log-Puiseux series in @h@ and read off the
-- @h^1@ coefficient. This produces a system that is mathematically
-- honest, compositionally elegant, and handles edge cases (poles,
-- branch points, non-analytic points) naturally.
--
-- = Log terms
--
-- The expansion type is 'LogPuiseuxSeries', which supports terms of
-- the form @c · h^p · log(h)^k@. Log terms arise naturally from:
--
-- * @log(h)@ — the expansion of @log(x)@ near @x = 0@
-- * @Ci(h)@, @Ei(h)@ — cosine and exponential integrals near 0
--
-- The @Li@ case requires @log(log(h))@, which is outside the current
-- type; it returns 'Unknown'.
--
-- = Taylor series depth
--
-- All Taylor series are truncated to 'depth' terms. The default is
-- 8, which is sufficient for derivative extraction (the @h^1@
-- coefficient) and limit computation with good numerical stability.
module LimCalc.Series.Expand
  ( -- * Expansion
    expand
    -- * Series operations
  , composeSeries
  , evalSeriesAt
  , powSeries
  , invertSeries
  , geometricSeries
  , binomialSeries
  , normalizeW
  , expandPow
  , expandPowR
    -- * Taylor series generators
  , sinTaylor
  , cosTaylor
  , expTaylor
  , logTaylor
  , erfTaylor
  , siTaylor
  , ciTaylor
  , eiTaylor
    -- * Utilities
  , depth
  , eulerGamma
  , factorial
  , ciCoeff
  , isPositiveIntegerExp
  , isOddInteger
  , maxInvertIters
  ) where

import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.Ratio (numerator, denominator)
import LimCalc.Core.Expr
import LimCalc.Series.Puiseux
import LimCalc.Core.Types
import LimCalc.Algebra.AlgNum

-- | Expand @f(x₀ + h)@ as a log-Puiseux series in @h@.
--
-- @point@ maps each variable name to its base-point value @x₀@.
-- @var@ is the expansion variable; all other variables in @point@
-- are treated as constants. The result is a 'LogPuiseuxSeries' in
-- the perturbation @h@.
expand :: Expr -> Point -> String -> ExpandResult

-- Constants
expand (Const c) _ _ = Right $ LogPuiseuxSeries [pureTerm 0 (fromRational (toRational c))]
expand Pi        _ _ = Right $ LogPuiseuxSeries [pureTerm 0 (fromRational (toRational (pi :: Double)))]
expand E         _ _ = Right $ LogPuiseuxSeries [pureTerm 0 (fromRational (toRational (exp 1 :: Double)))]
expand I         _ _ = Right $ LogPuiseuxSeries [pureTerm 0 algI]

-- Variable
expand (Var name) point var
  | name == var =
      let x0 = Map.findWithDefault algZero name point
      in Right $ stripZeros $ LogPuiseuxSeries
           [ pureTerm 0 x0
           , pureTerm 1 algOne
           ]
  | otherwise =
      case Map.lookup name point of
        Just c  -> Right $ LogPuiseuxSeries [pureTerm 0 c]
        Nothing -> Left $ Unknown $ "Variable " ++ name ++ " not in point"

-- Arithmetic
expand (Add f g) point var = do
  sf <- expand f point var
  sg <- expand g point var
  return $ addSeries sf sg

expand (Sub f g) point var = do
  sf <- expand f point var
  sg <- expand g point var
  return $ addSeries sf (scaleSeries (negate algOne) sg)

expand (Mul f g) point var = do
  sf <- expand f point var
  sg <- expand g point var
  return $ mulSeries sf sg

expand (Neg f) point var = do
  sf <- expand f point var
  return $ scaleSeries (negate algOne) sf

expand (Div f g) point var = do
  sf <- expand f point var
  sg <- expand g point var
  case leadingTermNZ (stripZeros sg) of
    Nothing -> Left $ Singularity "Division by zero series"
    Just _  -> return $ mulSeries sf (invertSeries sg)

expand (Pow f g) point var = expandPow f g point var

expand (Exp f) point var = do
  sf <- expand f point var
  return $ composeSeries expTaylor sf

-- | Log expansion producing log-Puiseux terms.
--
-- For @f(x₀ + h)@ with leading term @a · h^α@ (@α ≥ 0@):
--
-- @log(f) = log(a · h^α · (1 + w)) = log(a) + α·log(h) + log(1 + w)@
--
-- where @w = f\/(a·h^α) − 1@ has no constant term. @log(1+w)@ is
-- computed via Taylor series. The @α·log(h)@ term is a genuine
-- log-Puiseux term with @lpLog = 1@.
--
-- For @α = 0@: reduces to an ordinary Taylor expansion of @log@
-- around @x₀@, producing no @log(h)@ term.
--
-- For @α < 0@ (pole): returns 'Undefined' — @log@ of a pole is not
-- representable as a log-Puiseux series.
expand (Log f) point var = do
  sf <- expand f point var
  let sf' = stripZeros sf
  case leadingTermNZ sf' of
    Nothing -> Left $ Undefined "log of zero series"
    Just lt ->
      let alpha = lpExp lt
          a     = lpCoeff lt
      in if alpha < 0
           then Left $ Undefined "log of a series with a pole"
           else if algToDouble a <= 0 && alpha == 0
             then Left $ Undefined "log of non-positive value"
             else do
               let logA      = algLog a
                   constPart = LogPuiseuxSeries [pureTerm 0 logA]
                   logHPart  = if alpha == 0
                                 then LogPuiseuxSeries []
                                 else LogPuiseuxSeries
                                        [logTerm 0 1 (fromRational alpha)]
                   w          = normalizeW sf' lt alpha a
                   taylorPart = removeExp 0
                                  (composeSeries logTaylor
                                    (LogPuiseuxSeries [pureTerm 0 a]
                                      `addSeries` w))
               return $ addSeries constPart
                      $ addSeries logHPart taylorPart

expand (Sin f) point var = do
  sf <- expand f point var
  return $ composeSeries sinTaylor sf

expand (Cos f) point var = do
  sf <- expand f point var
  return $ composeSeries cosTaylor sf

-- | @erf@ is entire; its expansion at any point is a pure power series
-- via Taylor.
expand (Erf f) point var = do
  sf <- expand f point var
  return $ composeSeries erfTaylor sf

-- | @Si@ is entire; its expansion at any point is a pure power series
-- via Taylor.
expand (Si f) point var = do
  sf <- expand f point var
  return $ composeSeries siTaylor sf

-- | @Ci(x)@ near @x₀ ≠ 0@: analytic, expand via Taylor.
-- At @x₀ = 0@: @Ci(h) = γ + log(h) + h²\/4 − h⁴\/96 + …@
-- The @log(h)@ term is a genuine log-Puiseux term with @lpLog = 1@.
expand (Ci f) point var = do
  sf <- expand f point var
  let x0  = constantTerm sf
      x0d = algToDouble x0
  if abs x0d > 1e-12
    then return $ composeSeries ciTaylor sf
    else do
      let gammaLogH = LogPuiseuxSeries
                        [ pureTerm 0 (fromRational (toRational eulerGamma))
                        , logTerm  0 1 algOne
                        ]
          powerTerms = LogPuiseuxSeries $ take depth
            [ pureTerm (fromIntegral (2*n))
                (fromRational (toRational (ciCoeff n)))
            | n <- [1..] :: [Int] ]
      return $ addSeries gammaLogH powerTerms

-- | @Ei(x)@ near @x₀ ≠ 0@: analytic, expand via Taylor.
-- At @x₀ = 0@: @Ei(h) = γ + log(h) + h + h²\/4 + h³\/18 + …@
-- The @log(h)@ term is a genuine log-Puiseux term with @lpLog = 1@.
expand (Ei f) point var = do
  sf <- expand f point var
  let x0  = constantTerm sf
      x0d = algToDouble x0
  if abs x0d > 1e-12
    then return $ composeSeries eiTaylor sf
    else do
      let gammaLogH = LogPuiseuxSeries
                        [ pureTerm 0 (fromRational (toRational eulerGamma))
                        , logTerm  0 1 algOne
                        ]
          powerTerms = LogPuiseuxSeries $ take depth
            [ pureTerm (fromIntegral n)
                (fromRational (toRational (1.0 / (fromIntegral n * fromIntegral (factorial n) :: Double))))
            | n <- [1..] :: [Int] ]
      return $ addSeries gammaLogH powerTerms

-- | @Li@ requires @log(log(h))@ near @x = 0@, which is outside the
-- current log-Puiseux type.
expand (Li _) _ _ = Left $ Unknown
  "Li expansion at x=0 requires log(log(h)) -- outside log-Puiseux type; \
  \expansion at x0 /= 0,1 via Taylor not yet implemented"

expand (Arcsin _) _ _ = Left $ Unknown "Arcsin expansion not yet implemented"
expand (Arccos _) _ _ = Left $ Unknown "Arccos expansion not yet implemented"
expand (Arctan _) _ _ = Left $ Unknown "Arctan expansion not yet implemented"

-- | Expand @|f|@ at a point.
--
-- Determines the sign of the leading term of @f@ to decide whether
-- @|f| = f@ or @|f| = −f@ locally. Returns 'NonAnalytic' when @f@
-- vanishes at the point with odd-order or fractional-order leading
-- term (genuine kink), and 'DomainError' when the leading coefficient
-- is non-real.
expand (Abs f) point var = do
  sf <- expand f point var
  case leadingTermNZ (stripZeros sf) of
    Nothing ->
      Left $ NonAnalytic
        "Abs has no local Puiseux expansion: f vanishes identically \
        \to all computed orders at this point"
    Just lt
      | lpExp lt > 0 && isOddInteger (lpExp lt) ->
          Left $ NonAnalytic
            "Abs: f(x0)=0 with odd-order vanishing; genuine kink"
      | lpExp lt > 0 && not (isPositiveIntegerExp (lpExp lt)) ->
          Left $ NonAnalytic
            "Abs: f(x0)=0 with fractional leading term; not analytic"
      | otherwise ->
          let a = lpCoeff lt
          in if abs (algImagDouble a) > 1e-12
               then Left $ DomainError
                      "Abs of series with non-real leading coefficient"
               else if algToDouble a > 0
                      then Right sf
                      else Right $ scaleSeries (negate algOne) sf

-- Helpers

-- | True if @r@ is a positive integer (denominator 1, numerator > 0).
isPositiveIntegerExp :: Rational -> Bool
isPositiveIntegerExp r = denominator r == 1 && numerator r > 0

-- | True if @r@ is a positive odd integer.
isOddInteger :: Rational -> Bool
isOddInteger r = isPositiveIntegerExp r && odd (numerator r)

-- | Number of Taylor series terms to compute.
depth :: Int
depth = 8

-- | Euler-Mascheroni constant, used in @Ci@ and @Ei@ expansions at 0.
eulerGamma :: Double
eulerGamma = 0.5772156649015329

-- | Factorial.
factorial :: Int -> Int
factorial n = product [1..n]

-- | @n@th coefficient in the power-series part of @Ci(h)@ at @h = 0@.
ciCoeff :: Int -> Double
ciCoeff n =
  let sign = if odd n then 1.0 else -1.0
  in sign / (fromIntegral (2*n) * fromIntegral (factorial (2*n)))

-- | Compose a Taylor series @S@ with an expansion @E@.
--
-- Extracts the constant term @c₀@ of @E@, passes it to the Taylor
-- generator @taylorAt@, then evaluates the resulting series at the
-- perturbation @u = E − c₀@. Log terms in @u@ are handled correctly
-- by 'mulSeries'.
composeSeries :: (AlgNum -> LogPuiseuxSeries AlgNum)
              -> LogPuiseuxSeries AlgNum
              -> LogPuiseuxSeries AlgNum
composeSeries taylorAt e =
  let c0 = constantTerm e
      u  = removeExp 0 e
      s  = taylorAt c0
  in evalSeriesAt u s

-- | Evaluate @S = ∑ aₙ · tⁿ@ by substituting @t = u@.
evalSeriesAt :: LogPuiseuxSeries AlgNum
             -> LogPuiseuxSeries AlgNum
             -> LogPuiseuxSeries AlgNum
evalSeriesAt u (LogPuiseuxSeries ts) =
  foldr addSeries (LogPuiseuxSeries [])
    [ scaleSeries (lpCoeff t) (powSeries u (lpExp t))
    | t <- ts
    ]

-- | Raise a series to a rational power.
-- Only handles non-negative integer powers; returns the empty series otherwise.
powSeries :: LogPuiseuxSeries AlgNum -> Rational -> LogPuiseuxSeries AlgNum
powSeries _ 0 = LogPuiseuxSeries [pureTerm 0 algOne]
powSeries u n
  | n > 0 && denominator n == 1 =
      let k = fromIntegral (numerator n) :: Int
      in foldl mulSeries (LogPuiseuxSeries [pureTerm 0 algOne]) (replicate k u)
  | otherwise = LogPuiseuxSeries []

-- | Compute @w = s\/(a · h^α) − 1@: normalise a series by factoring
-- out its leading term and subtracting 1.
normalizeW :: LogPuiseuxSeries AlgNum
           -> LogPuiseuxTerm AlgNum
           -> Rational -> AlgNum
           -> LogPuiseuxSeries AlgNum
normalizeW (LogPuiseuxSeries ts) _lt alpha a =
  let shifted = [ LogPuiseuxTerm (lpCoeff t / a) (lpExp t - alpha) (lpLog t)
                | t <- ts ]
  in removeTerm 0 0 (LogPuiseuxSeries shifted)

-- | Invert a series: @1\/s@ via geometric series expansion.
invertSeries :: LogPuiseuxSeries AlgNum -> LogPuiseuxSeries AlgNum
invertSeries s =
  case leadingTermNZ (stripZeros s) of
    Nothing -> LogPuiseuxSeries []
    Just lt ->
      let alpha = lpExp lt
          a     = lpCoeff lt
          w     = normalizeW s lt alpha a
          negw  = scaleSeries (negate algOne) w
          geo   = geometricSeries negw (fromIntegral depth)
          shift = negate alpha
      in stripZeros $ truncateToOrder (fromIntegral depth)
           (shiftExponents shift (scaleSeries (recip a) geo))

-- | Geometric series @(1 + u)^(−1) = ∑ (−u)^n@, iterated by order.
geometricSeries :: LogPuiseuxSeries AlgNum -> Rational -> LogPuiseuxSeries AlgNum
geometricSeries u targetOrder =
  case leadingTermNZ (stripZeros u) of
    Nothing -> LogPuiseuxSeries [pureTerm 0 algOne]
    Just lt ->
      let gapOrder = lpExp lt
          neededN  = if gapOrder <= 0
                       then maxInvertIters
                       else min maxInvertIters
                              (ceiling (targetOrder / gapOrder) + 1)
          upows = take (neededN + 1) $
                    iterate (truncateToOrder targetOrder . mulSeries u)
                            (LogPuiseuxSeries [pureTerm 0 algOne])
      in truncateToOrder targetOrder $ foldr addSeries (LogPuiseuxSeries []) upows

-- | Hard cap on geometric-series iterations.
maxInvertIters :: Int
maxInvertIters = 200

-- | Dispatch the 'Pow' case on whether the exponent is a numeric constant.
expandPow :: Expr -> Expr -> Point -> String -> ExpandResult
expandPow f (Const r) point var       = expandPowR f (toRational r) point var
expandPow f (Neg (Const r)) point var = expandPowR f (toRational (-r)) point var
expandPow _ _ _ _                     = Left $ Unknown "Symbolic exponents not yet supported"

-- | Expand @f^r@ for rational @r@ via the generalised binomial series.
expandPowR :: Expr -> Rational -> Point -> String -> ExpandResult
expandPowR f r point var = do
  s <- expand f point var
  let s' = stripZeros s
  case leadingTermNZ s' of
    Nothing -> Right $ LogPuiseuxSeries []
    Just lt ->
      let alpha = lpExp lt
          a     = lpCoeff lt
          w     = normalizeW s' lt alpha a
          binom = binomialSeries r w (fromIntegral depth)
          scale = a ** fromRational r
          shift = alpha * r
      in Right $ stripZeros $ truncateToOrder (shift + fromIntegral depth)
           (shiftExponents shift (scaleSeries scale binom))

-- | Generalised binomial series @(1 + w)^r@, iterated by order.
binomialSeries :: Rational -> LogPuiseuxSeries AlgNum -> Rational -> LogPuiseuxSeries AlgNum
binomialSeries r w targetOrder =
  case leadingTermNZ (stripZeros w) of
    Nothing -> LogPuiseuxSeries [pureTerm 0 algOne]
    Just lt ->
      let gapOrder = lpExp lt
          neededN  = if gapOrder <= 0
                       then maxInvertIters
                       else min maxInvertIters
                              (ceiling (targetOrder / gapOrder) + 1)
          bcs   = binomCoeffs r neededN
          wpows = take (neededN + 1) $
                    iterate (truncateToOrder targetOrder . mulSeries w)
                            (LogPuiseuxSeries [pureTerm 0 algOne])
      in truncateToOrder targetOrder $ foldr addSeries (LogPuiseuxSeries [])
           [ scaleSeries c wp | (c, wp) <- zip bcs wpows ]

-- | Generalised binomial coefficients @C(r, 0) = 1@,
-- @C(r, k) = r(r−1)···(r−k+1)\/k!@.
binomCoeffs :: Rational -> Int -> [AlgNum]
binomCoeffs r n = take (n+1) $ scanl step algOne [0..n-1]
  where
    step acc k =
      let rk = fromRational r - fromIntegral k
          kk = fromIntegral (k+1 :: Int)
      in acc * rk / kk

------------------------------------------------------------------------
-- Taylor series generators
------------------------------------------------------------------------

-- | Taylor series of @sin@ around @x₀@.
sinTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
sinTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral n) (sinCoeff n) | n <- [0..] :: [Int] ]
  where
    facts = scanl (*) 1 [1..] :: [Int]
    sinCoeff n =
      let s = algSin x0; c = algCos x0
          bases = cycle [s, c, negate s, negate c]
      in bases !! n / fromIntegral (facts !! n)

-- | Taylor series of @cos@ around @x₀@.
cosTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
cosTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral n) (cosCoeff n) | n <- [0..] :: [Int] ]
  where
    facts = scanl (*) 1 [1..] :: [Int]
    cosCoeff n =
      let s = algSin x0; c = algCos x0
          bases = cycle [c, negate s, negate c, s]
      in bases !! n / fromIntegral (facts !! n)

-- | Taylor series of @exp@ around @x₀@.
expTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
expTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral n) (algExp x0 / fromIntegral f)
  | (n, f) <- zip [0..] (scanl (*) 1 [1..] :: [Int]) ]

-- | Taylor series of @log@ around @x₀ ≠ 0@, excluding the constant
-- term @log(x₀)@. Used by 'composeSeries' in the 'Log' expansion case.
logTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
logTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral n) (logCoeff n) | n <- [1..] :: [Int] ]
  where
    logCoeff n =
      let sign = if even n then negate algOne else algOne
      in sign / (fromIntegral n * x0 ^ n)

-- | Taylor series of @erf@ around @x₀@.
erfTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
erfTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral (2*n+1)) (erfCoeff n x0) | n <- [0..] :: [Int] ]
  where
    twoOverSqrtPi = 2.0 / sqrt pi :: Double
    erfCoeff n x =
      let sign   = if even n then algOne else negate algOne
          nfact  = fromIntegral (factorial n) :: AlgNum
          twon1  = fromIntegral (2*n+1) :: AlgNum
          xpow   = x ^ (2*n)
          pre    = fromRational (toRational twoOverSqrtPi)
      in pre * sign * xpow / (nfact * twon1)

-- | Taylor series of @Si@ around @x₀@.
siTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
siTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral (2*n+1)) (siCoeff n x0) | n <- [0..] :: [Int] ]
  where
    siCoeff n x =
      let sign  = if even n then algOne else negate algOne
          twon1 = 2*n+1
          denom = fromIntegral (twon1 * factorial twon1) :: AlgNum
          xpow  = x ^ (2*n)
      in sign * xpow / denom

-- | Taylor series of @Ci@ around @x₀ ≠ 0@.
-- At @x₀ = 0@ the @log(h)@ case is handled directly in 'expand'.
ciTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
ciTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral n) (ciTaylorCoeff n x0) | n <- [0..] :: [Int] ]
  where
    ciTaylorCoeff n x =
      let xd = algToDouble x
      in if n == 0
           then fromRational (toRational (ciValue xd))
           else fromRational (toRational (ciDerivN n xd / fromIntegral (factorial n)))
    ciValue xd = eulerGamma + log (abs xd) +
      sum [ (if odd k then 1 else -1) * xd^(2*k) /
              fromIntegral (2*k * factorial (2*k))
          | k <- [1..depth] :: [Int] ]
    ciDerivN 0 xd = ciValue xd
    ciDerivN _ xd = cos xd / xd

-- | Taylor series of @Ei@ around @x₀ ≠ 0@.
-- At @x₀ = 0@ the @log(h)@ case is handled directly in 'expand'.
eiTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
eiTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral n) (eiTaylorCoeff n x0) | n <- [0..] :: [Int] ]
  where
    eiTaylorCoeff n x =
      let xd = algToDouble x
      in if n == 0
           then fromRational (toRational (eiValue xd))
           else fromRational (toRational (eiDerivN n xd / fromIntegral (factorial n)))
    eiValue xd = eulerGamma + log (abs xd) +
      sum [ xd^k / fromIntegral (k * factorial k)
          | k <- [1..depth] :: [Int] ]
    eiDerivN 0 xd = eiValue xd
    eiDerivN _ xd = exp xd / xd