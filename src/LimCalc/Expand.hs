module LimCalc.Expand where

import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.Ratio (numerator, denominator)
import LimCalc.Expr
import LimCalc.Puiseux
import LimCalc.Types

-- | Expand f(point + h*var) as a Puiseux series in h.
-- point: the expansion point for all variables
-- var:   the variable we are expanding in (the h direction)
expand :: Expr -> Point -> String -> ExpandResult

-- Constants
expand (Const c) _ _ = Right $ PuiseuxSeries [PuiseuxTerm 0 c]

-- Pi and E
expand Pi _ _ = Right $ PuiseuxSeries [PuiseuxTerm 0 pi]
expand E  _ _ = Right $ PuiseuxSeries [PuiseuxTerm 0 (exp 1)]

-- Imaginary unit — placeholder
expand I  _ _ = Right $ PuiseuxSeries [PuiseuxTerm 0 0]  -- TODO: complex

-- Variable
expand (Var name) point var
  | name == var =
      -- expanding in this variable: x₀ + h
      let x0 = Map.findWithDefault 0 name point
      in Right $ stripZeros $ PuiseuxSeries
           [ PuiseuxTerm 0 x0
           , PuiseuxTerm 1 1.0
           ]
  | otherwise =
      -- not the expansion variable: treat as constant
      case Map.lookup name point of
        Just c  -> Right $ PuiseuxSeries [PuiseuxTerm 0 c]
        Nothing -> Left $ Unknown $ "Variable " ++ name ++ " not in point"

-- Addition
expand (Add f g) point var = do
  sf <- expand f point var
  sg <- expand g point var
  return $ addSeries sf sg

-- Subtraction
expand (Sub f g) point var = do
  sf <- expand f point var
  sg <- expand g point var
  return $ addSeries sf (scaleSeries (-1) sg)

-- Multiplication
expand (Mul f g) point var = do
  sf <- expand f point var
  sg <- expand g point var
  return $ mulSeries sf sg

-- Negation
expand (Neg f) point var = do
  sf <- expand f point var
  return $ scaleSeries (-1) sf

-- Division
expand (Div f g) point var = do
  sf <- expand f point var
  sg <- expand g point var
  case leadingTermNZ (stripZeros sg) of
    Nothing -> Left $ Singularity "Division by zero series"
    Just _  -> return $ mulSeries sf (invertSeries sg)

-- Power
expand (Pow f g) point var = expandPow f g point var

-- Exp
expand (Exp f) point var = do
  sf <- expand f point var
  return $ composeSeries expTaylor sf

-- Log
expand (Log f) point var = do
  sf <- expand f point var
  let c0 = constantTerm sf
  if c0 <= 0
    then Left $ Undefined $ "Log of non-positive value: " ++ show c0
    else return $ composeSeries logTaylor sf

-- Sin
expand (Sin f) point var = do
  sf <- expand f point var
  return $ composeSeries sinTaylor sf

-- Cos
expand (Cos f) point var = do
  sf <- expand f point var
  return $ composeSeries cosTaylor sf

-- Abs
expand (Abs f) point var = expandAbs f point var

-- | How many terms to compute
depth :: Int
depth = 8

-- | Compose a Taylor series S with expansion E
composeSeries :: (Double -> PuiseuxSeries) -> PuiseuxSeries -> PuiseuxSeries
composeSeries taylorAt e =
  let c0 = constantTerm e
      u  = removeTerm 0 e
      s  = taylorAt c0
  in evalSeriesAt u s

-- | Evaluate S = Σ aₙ·t^n by substituting t = u
evalSeriesAt :: PuiseuxSeries -> PuiseuxSeries -> PuiseuxSeries
evalSeriesAt u (PuiseuxSeries ts) =
  foldr addSeries zeroPuiseux
    [ scaleSeries (coeff t) (powSeries u (pExp t))
    | t <- ts
    ]

-- | Zero series
zeroPuiseux :: PuiseuxSeries
zeroPuiseux = PuiseuxSeries []

-- | Raise a series to a rational power (integer powers only for now)
powSeries :: PuiseuxSeries -> Rational -> PuiseuxSeries
powSeries _ 0 = PuiseuxSeries [PuiseuxTerm 0 1.0]
powSeries u n
  | n > 0 && denominator n == 1 =
      let k = fromIntegral (numerator n) :: Int
      in foldl mulSeries (PuiseuxSeries [PuiseuxTerm 0 1.0]) (replicate k u)
  | otherwise = PuiseuxSeries []

-- | Remove term with given exponent
removeTerm :: Rational -> PuiseuxSeries -> PuiseuxSeries
removeTerm e (PuiseuxSeries ts) =
  PuiseuxSeries $ filter (\t -> pExp t /= e) ts

-- | Get constant term
constantTerm :: PuiseuxSeries -> Double
constantTerm (PuiseuxSeries []) = 0
constantTerm (PuiseuxSeries (t:_))
  | pExp t == 0 = coeff t
  | otherwise   = 0

-- | Truncate to n terms
truncateSeries :: Int -> PuiseuxSeries -> PuiseuxSeries
truncateSeries n (PuiseuxSeries ts) = PuiseuxSeries (take n ts)

-- | Invert a series: 1/s
invertSeries :: PuiseuxSeries -> PuiseuxSeries
invertSeries s =
  case leadingTermNZ (stripZeros s) of
    Nothing -> PuiseuxSeries []
    Just lt ->
      let alpha = pExp lt
          a     = coeff lt
          w     = truncateSeries depth (normalizeW s lt alpha a)
          negw  = scaleSeries (-1) w
          geo   = geometricSeries negw
          scale = 1.0 / a
          shift = negate alpha
      in stripZeros $ truncateSeries depth
           (shiftExponents shift (scaleSeries scale geo))

-- | Geometric series 1/(1+w) = Σ (-w)^n
geometricSeries :: PuiseuxSeries -> PuiseuxSeries
geometricSeries u =
  let upows = take (depth+1) $ iterate (truncateSeries depth . mulSeries u)
                                       (PuiseuxSeries [PuiseuxTerm 0 1.0])
  in truncateSeries depth $ foldr addSeries zeroPuiseux upows

-- | Handle Pow case
expandPow :: Expr -> Expr -> Point -> String -> ExpandResult
expandPow f (Const r) point var     = expandPowR f (toRational r) point var
expandPow f (Neg (Const r)) point var = expandPowR f (toRational (-r)) point var
expandPow _ _ _ _                   = Left $ Unknown "Symbolic exponents not yet supported"

-- | Expand f^r for rational r
expandPowR :: Expr -> Rational -> Point -> String -> ExpandResult
expandPowR f r point var = do
  s <- expand f point var
  let s' = stripZeros $ truncateSeries depth s
  case leadingTermNZ s' of
    Nothing -> Right $ PuiseuxSeries []
    Just lt ->
      let alpha = pExp lt
          a     = coeff lt
          w     = truncateSeries depth (normalizeW s' lt alpha a)
          binom = binomialSeries r w
          scale = a ** fromRational r
          shift = alpha * r
      in Right $ stripZeros $ truncateSeries depth
           (shiftExponents shift (scaleSeries scale binom))

-- | Compute w = s/(a*h^alpha) - 1
normalizeW :: PuiseuxSeries -> PuiseuxTerm -> Rational -> Double -> PuiseuxSeries
normalizeW (PuiseuxSeries ts) _lt alpha a =
  let shifted = [ PuiseuxTerm (pExp t - alpha) (coeff t / a) | t <- ts ]
  in removeTerm 0 (PuiseuxSeries shifted)

-- | Binomial series (1+w)^r
binomialSeries :: Rational -> PuiseuxSeries -> PuiseuxSeries
binomialSeries r w =
  let bcs   = binomCoeffs r depth
      wpows = take (depth+1) $ iterate (truncateSeries depth . mulSeries w)
                                       (PuiseuxSeries [PuiseuxTerm 0 1.0])
  in truncateSeries depth $ foldr addSeries zeroPuiseux
       [ scaleSeries c wp
       | (c, wp) <- zip bcs wpows
       ]

-- | Generalized binomial coefficients
binomCoeffs :: Rational -> Int -> [Double]
binomCoeffs r n = take (n+1) $ scanl step 1.0 ([0..] :: [Int])
  where
    step acc k = acc * (fromRational r - fromIntegral k) / fromIntegral (k+1)

-- | Shift all exponents by delta
shiftExponents :: Rational -> PuiseuxSeries -> PuiseuxSeries
shiftExponents delta (PuiseuxSeries ts) =
  PuiseuxSeries [ PuiseuxTerm (pExp t + delta) (coeff t) | t <- ts ]

-- | Handle Abs case
expandAbs :: Expr -> Point -> String -> ExpandResult
expandAbs _ _ _ = Left $ Unknown "Abs expansion not yet implemented"

sinTaylor :: Double -> PuiseuxSeries
sinTaylor x0 = PuiseuxSeries $ take depth
  [ PuiseuxTerm (fromIntegral n) (sinCoeff n)
  | n <- [0..] :: [Int] ]
  where
    s = sin x0
    c = cos x0
    signs = cycle [s, c, -s, -c]
    facts = scanl (*) 1 [1..] :: [Int]
    sinCoeff n = (signs !! n) / fromIntegral (facts !! n)

cosTaylor :: Double -> PuiseuxSeries
cosTaylor x0 = PuiseuxSeries $ take depth
  [ PuiseuxTerm (fromIntegral n) (cosCoeff n)
  | n <- [0..] :: [Int] ]
  where
    s = sin x0
    c = cos x0
    signs = cycle [c, -s, -c, s]
    facts = scanl (*) 1 [1..] :: [Int]
    cosCoeff n = (signs !! n) / fromIntegral (facts !! n)

expTaylor :: Double -> PuiseuxSeries
expTaylor x0 = PuiseuxSeries $ take depth
  [ PuiseuxTerm (fromIntegral n) (exp x0 / fromIntegral f)
  | (n, f) <- zip ([0..] :: [Int]) facts ]
  where
    facts = scanl (*) 1 [1..] :: [Int]

logTaylor :: Double -> PuiseuxSeries
logTaylor x0 = PuiseuxSeries $ take depth $
  PuiseuxTerm 0 (log x0) :
  [ PuiseuxTerm (fromIntegral n) (logCoeff n)
  | n <- [1..] :: [Int] ]
  where
    logCoeff n =
      let sign = if even n then -1 else 1
      in sign / (fromIntegral n * x0 ^ n)