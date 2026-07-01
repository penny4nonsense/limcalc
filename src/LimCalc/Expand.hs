module LimCalc.Expand where

import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.Ratio (numerator, denominator)
import LimCalc.Expr
import LimCalc.Puiseux
import LimCalc.Types
import LimCalc.AlgNum

-- | Expand f(point + h*var) as a log-Puiseux series in h
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
-- For f(x0 + h) with leading term a*h^alpha (alpha >= 0):
--
--   log(f) = log(a * h^alpha * (1 + w))
--           = log(a) + alpha*log(h) + log(1 + w)
--
-- where w = f/(a*h^alpha) - 1 has no constant term.
-- log(1+w) is computed via Taylor series since w->0.
-- The alpha*log(h) term is a genuine log-Puiseux term with lpLog=1.
--
-- For alpha=0 (f(x0) != 0): reduces to ordinary Taylor expansion of
-- log around x0, no log(h) term produced.
--
-- For alpha<0 (pole): return Undefined -- log of a pole is not a
-- log-Puiseux series in any standard sense.
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

-- | erf(x) = (2/sqrt(pi)) * (x - x^3/3 + x^5/10 - x^7/42 + ...)
-- Entire function; pure power series at any point via Taylor.
expand (Erf f) point var = do
  sf <- expand f point var
  return $ composeSeries erfTaylor sf

-- | Si(x) = x - x^3/18 + x^5/600 - x^7/35280 + ...
-- Entire function; pure power series at any point via Taylor.
expand (Si f) point var = do
  sf <- expand f point var
  return $ composeSeries siTaylor sf

-- | Ci(x) near x0 /= 0: analytic, expand via Taylor.
-- At x0 = 0: Ci(h) = gamma + log(h) + h^2/4 - h^4/96 + ...
-- The log(h) term is a genuine log-Puiseux term with lpLog=1.
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

-- | Ei(x) near x0 /= 0: analytic, expand via Taylor.
-- At x0 = 0: Ei(h) = gamma + log(h) + h + h^2/4 + h^3/18 + ...
-- The log(h) term is a genuine log-Puiseux term with lpLog=1.
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

-- | Li(x) requires log(log(x)) near x=0, outside current type.
expand (Li _) _ _ = Left $ Unknown
  "Li expansion at x=0 requires log(log(h)) -- outside log-Puiseux type; \
  \expansion at x0 /= 0,1 via Taylor not yet implemented"

expand (Arcsin _) _ _ = Left $ Unknown "Arcsin expansion not yet implemented"
expand (Arccos _) _ _ = Left $ Unknown "Arccos expansion not yet implemented"
expand (Arctan _) _ _ = Left $ Unknown "Arctan expansion not yet implemented"

-- | Abs expansion (unchanged logic, updated types)
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

isPositiveIntegerExp :: Rational -> Bool
isPositiveIntegerExp r = denominator r == 1 && numerator r > 0

isOddInteger :: Rational -> Bool
isOddInteger r = isPositiveIntegerExp r && odd (numerator r)

depth :: Int
depth = 8

eulerGamma :: Double
eulerGamma = 0.5772156649015329

factorial :: Int -> Int
factorial n = product [1..n]

ciCoeff :: Int -> Double
ciCoeff n =
  let sign = if odd n then 1.0 else -1.0
  in sign / (fromIntegral (2*n) * fromIntegral (factorial (2*n)))

-- | Compose a Taylor series S with expansion E.
-- For log-Puiseux: only valid when E has no log terms at constant
-- position (i.e. the Taylor series is in a pure-power argument).
-- composeSeries strips log terms from the constant and passes only
-- the pure-power part to the Taylor generator; log terms in the
-- perturbation u are handled by mulSeries (closed under log terms).
composeSeries :: (AlgNum -> LogPuiseuxSeries AlgNum)
              -> LogPuiseuxSeries AlgNum
              -> LogPuiseuxSeries AlgNum
composeSeries taylorAt e =
  let c0 = constantTerm e
      u  = removeExp 0 e
      s  = taylorAt c0
  in evalSeriesAt u s

evalSeriesAt :: LogPuiseuxSeries AlgNum
             -> LogPuiseuxSeries AlgNum
             -> LogPuiseuxSeries AlgNum
evalSeriesAt u (LogPuiseuxSeries ts) =
  foldr addSeries (LogPuiseuxSeries [])
    [ scaleSeries (lpCoeff t) (powSeries u (lpExp t))
    | t <- ts
    ]

powSeries :: LogPuiseuxSeries AlgNum -> Rational -> LogPuiseuxSeries AlgNum
powSeries _ 0 = LogPuiseuxSeries [pureTerm 0 algOne]
powSeries u n
  | n > 0 && denominator n == 1 =
      let k = fromIntegral (numerator n) :: Int
      in foldl mulSeries (LogPuiseuxSeries [pureTerm 0 algOne]) (replicate k u)
  | otherwise = LogPuiseuxSeries []

normalizeW :: LogPuiseuxSeries AlgNum
           -> LogPuiseuxTerm AlgNum
           -> Rational -> AlgNum
           -> LogPuiseuxSeries AlgNum
normalizeW (LogPuiseuxSeries ts) _lt alpha a =
  let shifted = [ LogPuiseuxTerm (lpCoeff t / a) (lpExp t - alpha) (lpLog t)
                | t <- ts ]
  in removeTerm 0 0 (LogPuiseuxSeries shifted)

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

maxInvertIters :: Int
maxInvertIters = 200

expandPow :: Expr -> Expr -> Point -> String -> ExpandResult
expandPow f (Const r) point var       = expandPowR f (toRational r) point var
expandPow f (Neg (Const r)) point var = expandPowR f (toRational (-r)) point var
expandPow _ _ _ _                     = Left $ Unknown "Symbolic exponents not yet supported"

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

binomCoeffs :: Rational -> Int -> [AlgNum]
binomCoeffs r n = take (n+1) $ scanl step algOne [0..n-1]
  where
    step acc k =
      let rk = fromRational r - fromIntegral k
          kk = fromIntegral (k+1 :: Int)
      in acc * rk / kk

-- Taylor series generators

sinTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
sinTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral n) (sinCoeff n) | n <- [0..] :: [Int] ]
  where
    facts = scanl (*) 1 [1..] :: [Int]
    sinCoeff n =
      let s = algSin x0; c = algCos x0
          bases = cycle [s, c, negate s, negate c]
      in bases !! n / fromIntegral (facts !! n)

cosTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
cosTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral n) (cosCoeff n) | n <- [0..] :: [Int] ]
  where
    facts = scanl (*) 1 [1..] :: [Int]
    cosCoeff n =
      let s = algSin x0; c = algCos x0
          bases = cycle [c, negate s, negate c, s]
      in bases !! n / fromIntegral (facts !! n)

expTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
expTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral n) (algExp x0 / fromIntegral f)
  | (n, f) <- zip [0..] (scanl (*) 1 [1..] :: [Int]) ]

-- | Taylor series of log around x0 /= 0, excluding the constant term.
-- The constant log(x0) is handled separately in the Log expansion;
-- this provides only the h^1, h^2, ... terms for composeSeries.
logTaylor :: AlgNum -> LogPuiseuxSeries AlgNum
logTaylor x0 = LogPuiseuxSeries $ take depth
  [ pureTerm (fromIntegral n) (logCoeff n) | n <- [1..] :: [Int] ]
  where
    logCoeff n =
      let sign = if even n then negate algOne else algOne
      in sign / (fromIntegral n * x0 ^ n)

-- | erf Taylor: (2/sqrt(pi)) * sum_{n=0}^{inf} (-1)^n x^(2n+1) / (n! * (2n+1))
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

-- | Si Taylor: sum_{n=0}^{inf} (-1)^n x^(2n+1) / ((2n+1) * (2n+1)!)
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

-- | Ci Taylor around x0 /= 0 (analytic away from 0)
-- At x0=0 the log(h) case is handled directly in expand (Ci f).
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

-- | Ei Taylor around x0 /= 0
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