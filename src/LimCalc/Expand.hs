module LimCalc.Expand where

import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.Ratio (numerator, denominator)
import LimCalc.Expr
import LimCalc.Puiseux
import LimCalc.Types
import LimCalc.AlgNum

-- | Expand f(point + h*var) as a Puiseux series in h with AlgNum coefficients
expand :: Expr -> Point -> String -> ExpandResult

-- Constants
expand (Const c) _ _ = Right $ PuiseuxSeries [PuiseuxTerm 0 (fromRational (toRational c))]
expand Pi        _ _ = Right $ PuiseuxSeries [PuiseuxTerm 0 (fromRational (toRational (pi :: Double)))]
expand E         _ _ = Right $ PuiseuxSeries [PuiseuxTerm 0 (fromRational (toRational (exp 1 :: Double)))]
expand I         _ _ = Right $ PuiseuxSeries [PuiseuxTerm 0 algI]

-- Variable
expand (Var name) point var
  | name == var =
      let x0 = Map.findWithDefault algZero name point
      in Right $ stripZeros $ PuiseuxSeries
           [ PuiseuxTerm 0 x0
           , PuiseuxTerm 1 algOne
           ]
  | otherwise =
      case Map.lookup name point of
        Just c  -> Right $ PuiseuxSeries [PuiseuxTerm 0 c]
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

expand (Log f) point var = do
  sf <- expand f point var
  let x0 = constantTerm sf
  if algToDouble x0 <= 0
    then Left $ Undefined "log of non-positive value"
    else return $ composeSeries logTaylor sf

expand (Sin f) point var = do
  sf <- expand f point var
  return $ composeSeries sinTaylor sf

expand (Cos f) point var = do
  sf <- expand f point var
  return $ composeSeries cosTaylor sf

-- | Expand |f| at the given point.
--
-- Near h=0, |f(x0+h)| is determined entirely by the sign of f's
-- leading nonzero term, since that term dominates for small h:
--   - leading coefficient > 0  =>  |f| = f locally, series unchanged
--   - leading coefficient < 0  =>  |f| = -f locally, series negated
--   - leading coefficient non-real (nonzero imaginary part) =>
--       |.| has no sensible real-analytic local expansion here
--   - f's series is identically zero (no nonzero term at all, e.g.
--       Abs (Var "x") at x0=0) => this is a genuine analytic kink.
--       No Puiseux series in h represents |h| near h=0 (it agrees
--       with neither h nor -h on both sides), so the engine reports
--       this as NonAnalytic rather than silently returning some
--       series (which would let a derivative computation downstream
--       report a wrong numeric answer instead of correctly failing).
expand (Abs f) point var = do
  sf <- expand f point var
  case leadingTermNZ (stripZeros sf) of
    Nothing ->
      Left $ NonAnalytic
        "Abs has no local Puiseux expansion: f vanishes identically \
        \to all computed orders at this point, indicating a genuine \
        \kink (e.g. |x| at x=0) where no derivative exists"
    Just lt
      -- f(x0) = 0 (leading exponent > 0, i.e. the constant term was
      -- zero and got stripped). Whether this is a genuine kink for
      -- |f| depends on whether f actually changes sign as h crosses
      -- 0, not merely on whether f(x0) = 0:
      --   - odd-order vanishing (e.g. f ~ c*h, c*h^3, ...) means f
      --     changes sign -> |f| has a true corner here, no single
      --     local Puiseux series represents it (e.g. |x| at x=0).
      --   - even-order vanishing (e.g. f ~ c*h^2) means f does NOT
      --     change sign -> |f| = +-f locally is still smooth (e.g.
      --     |x^2| at x=0 is just x^2).
      --   - fractional-order vanishing isn't a single well-defined
      --     real direction crossing in the usual sense either; treat
      --     conservatively as a kink (NonAnalytic) rather than guess.
      | pExp lt > 0 && isOddInteger (pExp lt) ->
          Left $ NonAnalytic
            "Abs has no local Puiseux expansion: f(x0) = 0 and f \
            \changes sign here (odd-order vanishing), so |f| has a \
            \genuine kink at this point (no derivative exists)"
      | pExp lt > 0 && not (isPositiveIntegerExp (pExp lt)) ->
          Left $ NonAnalytic
            "Abs has no local Puiseux expansion: f(x0) = 0 with a \
            \fractional-order leading term, so the sign behavior of \
            \f near this point is not well-defined as a single \
            \real direction; |f| is not analytic here"
      | otherwise ->
          let a = coeff lt
          in if abs (algImagDouble a) > 1e-12
               then Left $ DomainError
                      "Abs of a series with non-real leading coefficient"
               else if algToDouble a > 0
                      then Right sf
                      else Right $ scaleSeries (negate algOne) sf

-- | Local Puiseux expansion of erf/li/Si/Ci/Ei is not yet
-- implemented -- these are currently only produced as Risch
-- integration results (see Risch.hs's special-function recognition),
-- not locally expandable. Explicit Unknown errors here rather than
-- leaving the pattern incomplete (which would crash) or silently
-- producing wrong output.
expand (Erf _) _ _ = Left $ Unknown "Erf expansion not yet implemented"
expand (Li _)  _ _ = Left $ Unknown "Li expansion not yet implemented"
expand (Si _)  _ _ = Left $ Unknown "Si expansion not yet implemented"
expand (Ci _)  _ _ = Left $ Unknown "Ci expansion not yet implemented"
expand (Ei _)  _ _ = Left $ Unknown "Ei expansion not yet implemented"
expand (Arcsin _) _ _ = Left $ Unknown "Arcsin expansion not yet implemented"
expand (Arccos _) _ _ = Left $ Unknown "Arccos expansion not yet implemented"
expand (Arctan _) _ _ = Left $ Unknown "Arctan expansion not yet implemented"

-- | True if r is a positive integer (denominator 1, numerator > 0).
isPositiveIntegerExp :: Rational -> Bool
isPositiveIntegerExp r = denominator r == 1 && numerator r > 0

-- | True if r is a positive odd integer. Assumes isPositiveIntegerExp
-- r already holds where this is used.
isOddInteger :: Rational -> Bool
isOddInteger r = isPositiveIntegerExp r && odd (numerator r)

-- | How many terms
depth :: Int
depth = 8

-- | Compose a Taylor series S with expansion E
composeSeries :: (AlgNum -> PuiseuxSeries AlgNum) -> PuiseuxSeries AlgNum -> PuiseuxSeries AlgNum
composeSeries taylorAt e =
  let c0 = constantTerm e
      u  = removeTerm 0 e
      s  = taylorAt c0
  in evalSeriesAt u s

-- | Evaluate S = Σ aₙ·t^n by substituting t = u
evalSeriesAt :: PuiseuxSeries AlgNum -> PuiseuxSeries AlgNum -> PuiseuxSeries AlgNum
evalSeriesAt u (PuiseuxSeries ts) =
  foldr addSeries (PuiseuxSeries [])
    [ scaleSeries (coeff t) (powSeries u (pExp t))
    | t <- ts
    ]

-- | Raise a series to a rational power (integer powers only for now)
powSeries :: PuiseuxSeries AlgNum -> Rational -> PuiseuxSeries AlgNum
powSeries _ 0 = PuiseuxSeries [PuiseuxTerm 0 algOne]
powSeries u n
  | n > 0 && denominator n == 1 =
      let k = fromIntegral (numerator n) :: Int
      in foldl mulSeries (PuiseuxSeries [PuiseuxTerm 0 algOne]) (replicate k u)
  | otherwise = PuiseuxSeries []

-- | Remove term with given exponent
removeTerm :: Rational -> PuiseuxSeries AlgNum -> PuiseuxSeries AlgNum
removeTerm e (PuiseuxSeries ts) =
  PuiseuxSeries $ filter (\t -> pExp t /= e) ts

-- | Truncate to n terms
truncateSeries :: Int -> PuiseuxSeries AlgNum -> PuiseuxSeries AlgNum
truncateSeries n (PuiseuxSeries ts) = PuiseuxSeries (take n ts)

-- | Invert a series: 1/s
invertSeries :: PuiseuxSeries AlgNum -> PuiseuxSeries AlgNum
invertSeries s =
  case leadingTermNZ (stripZeros s) of
    Nothing -> PuiseuxSeries []
    Just lt ->
      let alpha = pExp lt
          a     = coeff lt
          w     = normalizeW s lt alpha a
          negw  = scaleSeries (negate algOne) w
          -- We want depth orders of h *past* the leading term in
          -- the final result, i.e. orders up to alpha + depth in
          -- the original series. After dividing through by h^alpha
          -- (the normalizeW step), that target becomes plain
          -- `depth` in w/h-space.
          geo   = geometricSeries negw (fromIntegral depth)
          shift = negate alpha
      in stripZeros $ truncateToOrder (fromIntegral depth)
           (shiftExponents shift (scaleSeries (recip a) geo))

-- | Geometric series (1+u)^{-1} = sum (-u)^n, truncated to reach a
-- target ORDER in h (not a fixed count of terms).
--
-- Previously this iterated a fixed number of powers of u (depth+1
-- of them), which is correct when u's own leading exponent is 1,
-- but silently under-reaches whenever u is sparse: if u's leading
-- exponent is delta (e.g. delta=3 for u ~ h^3, or delta=1/2 for
-- fractional gaps), then depth powers of u only reach order
-- depth*delta in h, while terms below that order (other than
-- multiples of delta) are simply absent with no error raised.
-- That's invisible until a downstream consumer (e.g. Div, via
-- mulSeries) needed a term at an order this silently never filled.
--
-- Iterating until the order reached exceeds the target instead
-- makes "how far this series is trustworthy" an explicit, order-
-- based invariant matching the rest of the engine (truncateSeries,
-- composeSeries, etc.), at the cost of needing more terms when u is
-- sparse. A hard iteration cap guards against runaway computation
-- if u's leading exponent is pathologically close to zero.
geometricSeries :: PuiseuxSeries AlgNum -> Rational -> PuiseuxSeries AlgNum
geometricSeries u targetOrder =
  case leadingTermNZ (stripZeros u) of
    Nothing -> PuiseuxSeries [PuiseuxTerm 0 algOne]  -- u = 0: (1+0)^-1 = 1
    Just lt ->
      let gapOrder = pExp lt
          -- Number of powers of u needed so that n*gapOrder exceeds
          -- targetOrder. Capped well above any reasonable target to
          -- avoid looping forever if gapOrder is pathologically tiny.
          neededN  = if gapOrder <= 0
                       then maxInvertIters  -- shouldn't happen: u has no zero-order term by construction (normalizeW strips it)
                       else min maxInvertIters
                              (ceiling (targetOrder / gapOrder) + 1)
          upows = take (neededN + 1) $
                    iterate (truncateToOrder targetOrder . mulSeries u)
                            (PuiseuxSeries [PuiseuxTerm 0 algOne])
      in truncateToOrder targetOrder $ foldr addSeries (PuiseuxSeries []) upows

-- | Hard cap on geometric-series iterations, guarding against
-- runaway computation when the series being inverted has a
-- pathologically small gap between its leading and next term.
maxInvertIters :: Int
maxInvertIters = 200

-- | Truncate a series to terms at or below a given order (exponent).
truncateToOrder :: Rational -> PuiseuxSeries AlgNum -> PuiseuxSeries AlgNum
truncateToOrder maxOrder (PuiseuxSeries ts) =
  PuiseuxSeries (filter (\t -> pExp t <= maxOrder) ts)

-- | Handle Pow case
expandPow :: Expr -> Expr -> Point -> String -> ExpandResult
expandPow f (Const r) point var     = expandPowR f (toRational r) point var
expandPow f (Neg (Const r)) point var = expandPowR f (toRational (-r)) point var
expandPow _ _ _ _                   = Left $ Unknown "Symbolic exponents not yet supported"

-- | Expand f^r for rational r
expandPowR :: Expr -> Rational -> Point -> String -> ExpandResult
expandPowR f r point var = do
  s <- expand f point var
  let s' = stripZeros s
  case leadingTermNZ s' of
    Nothing -> Right $ PuiseuxSeries []
    Just lt ->
      let alpha = pExp lt
          a     = coeff lt
          w     = normalizeW s' lt alpha a
          binom = binomialSeries r w (fromIntegral depth)
          scale = a ** fromRational r
          shift = alpha * r
      in Right $ stripZeros $ truncateToOrder (shift + fromIntegral depth)
           (shiftExponents shift (scaleSeries scale binom))

-- | Compute w = s/(a*h^alpha) - 1
normalizeW :: PuiseuxSeries AlgNum -> PuiseuxTerm AlgNum -> Rational -> AlgNum -> PuiseuxSeries AlgNum
normalizeW (PuiseuxSeries ts) _lt alpha a =
  let shifted = [ PuiseuxTerm (pExp t - alpha) (coeff t / a) | t <- ts ]
  in removeTerm 0 (PuiseuxSeries shifted)

-- | Binomial series (1+w)^r, truncated to reach a target ORDER in h.
--
-- Same fix as geometricSeries: a fixed count of w-powers silently
-- under-reaches when w is sparse (large or fractional gap between
-- its leading and next term). Iterate by order instead.
binomialSeries :: Rational -> PuiseuxSeries AlgNum -> Rational -> PuiseuxSeries AlgNum
binomialSeries r w targetOrder =
  case leadingTermNZ (stripZeros w) of
    Nothing -> PuiseuxSeries [PuiseuxTerm 0 algOne]  -- w = 0: (1+0)^r = 1
    Just lt ->
      let gapOrder = pExp lt
          neededN  = if gapOrder <= 0
                       then maxInvertIters
                       else min maxInvertIters
                              (ceiling (targetOrder / gapOrder) + 1)
          bcs   = binomCoeffs r neededN
          wpows = take (neededN + 1) $
                    iterate (truncateToOrder targetOrder . mulSeries w)
                            (PuiseuxSeries [PuiseuxTerm 0 algOne])
      in truncateToOrder targetOrder $ foldr addSeries (PuiseuxSeries [])
           [ scaleSeries c wp
           | (c, wp) <- zip bcs wpows
           ]

-- | Generalized binomial coefficients
binomCoeffs :: Rational -> Int -> [AlgNum]
binomCoeffs r n = take (n+1) $ scanl step algOne [0..n-1]
  where
    step acc k =
      let rk = fromRational r - fromIntegral k
          kk = fromIntegral (k+1 :: Int)
      in acc * rk / kk

-- | Shift all exponents by delta
shiftExponents :: Rational -> PuiseuxSeries AlgNum -> PuiseuxSeries AlgNum
shiftExponents delta (PuiseuxSeries ts) =
  PuiseuxSeries [ PuiseuxTerm (pExp t + delta) (coeff t) | t <- ts ]

-- | Taylor series of sin around x₀
sinTaylor :: AlgNum -> PuiseuxSeries AlgNum
sinTaylor x0 = PuiseuxSeries $ take depth
  [ PuiseuxTerm (fromIntegral n) (sinCoeff n)
  | n <- [0..] :: [Int] ]
  where
    facts  = scanl (*) 1 [1..] :: [Int]
    sinCoeff n =
      let s      = algSin x0
          c      = algCos x0
          bases  = cycle [s, c, negate s, negate c]
          base   = bases !! n
          factor = fromIntegral (facts !! n) :: AlgNum
      in base / factor

-- | Taylor series of cos around x₀
cosTaylor :: AlgNum -> PuiseuxSeries AlgNum
cosTaylor x0 = PuiseuxSeries $ take depth
  [ PuiseuxTerm (fromIntegral n) (cosCoeff n)
  | n <- [0..] :: [Int] ]
  where
    facts  = scanl (*) 1 [1..] :: [Int]
    cosCoeff n =
      let s      = algSin x0
          c      = algCos x0
          bases  = cycle [c, negate s, negate c, s]
          base   = bases !! n
          factor = fromIntegral (facts !! n) :: AlgNum
      in base / factor

-- | Taylor series of exp around x₀
expTaylor :: AlgNum -> PuiseuxSeries AlgNum
expTaylor x0 = PuiseuxSeries $ take depth
  [ PuiseuxTerm (fromIntegral n) (algExp x0 / fromIntegral f)
  | (n, f) <- zip [0..] (scanl (*) 1 [1..] :: [Int]) ]

-- | Taylor series of log around x₀
logTaylor :: AlgNum -> PuiseuxSeries AlgNum
logTaylor x0 = PuiseuxSeries $ take depth $
  PuiseuxTerm 0 (algLog x0) :
  [ PuiseuxTerm (fromIntegral n) (logCoeff n)
  | n <- [1..] :: [Int] ]
  where
    logCoeff n =
      let sign = if even n then negate algOne else algOne
          factor = fromIntegral n :: AlgNum
          xpow  = x0 ^ n
      in sign / (factor * xpow)