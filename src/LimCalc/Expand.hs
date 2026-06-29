module LimCalc.Expand where

import Data.Ratio (numerator, denominator)
import LimCalc.Expr
import LimCalc.Puiseux

-- | Expand f(x₀ + h) as a Puiseux series in h.
-- This is the core engine — all calculus operations derive from this.
expand :: Expr -> Double -> PuiseuxSeries

-- Constants are trivially constant series
expand (Const c) _ = PuiseuxSeries [PuiseuxTerm 0 c]

-- Pi and E as numeric constants for now
expand Pi _ = PuiseuxSeries [PuiseuxTerm 0 pi]
expand E  _ = PuiseuxSeries [PuiseuxTerm 0 (exp 1)]

-- Imaginary unit — placeholder, will need AlgNum
expand I  _ = PuiseuxSeries [PuiseuxTerm 0 0]  -- TODO: complex

-- Var x at x₀: f(x₀ + h) = x₀ + h
-- So the series is x₀*h^0 + 1*h^1
expand (Var _) x0 = stripZeros $ PuiseuxSeries
  [ PuiseuxTerm 0 x0
  , PuiseuxTerm 1 1.0
  ]

-- Addition: expand both and add series
expand (Add f g) x0 = addSeries (expand f x0) (expand g x0)

-- Subtraction
expand (Sub f g) x0 = addSeries (expand f x0) (scaleSeries (-1) (expand g x0))

-- Multiplication: expand both and multiply series
expand (Mul f g) x0 = mulSeries (expand f x0) (expand g x0)

-- Negation
expand (Neg f) x0 = scaleSeries (-1) (expand f x0)

-- Division: f/g = f * g^(-1)
expand (Div f g) x0 = mulSeries (expand f x0) (invertSeries (expand g x0))

-- Power: f^g -- the hard case, handled separately
expand (Pow f g) x0 = expandPow f g x0

-- Exp: e^f(x₀+h)
expand (Exp f) x0 = composeSeries expTaylor (expand f x0)

-- Log: ln(f(x₀+h))
expand (Log f) x0 = composeSeries logTaylor (expand f x0)

-- Sin: sin(f(x₀+h))
expand (Sin f) x0 = composeSeries sinTaylor (expand f x0)

-- Cos: cos(f(x₀+h))
expand (Cos f) x0 = composeSeries cosTaylor (expand f x0)

-- Abs: non-analytic, handle separately
expand (Abs f) x0 = expandAbs f x0

-- | How many terms to compute in series expansions
depth :: Int
depth = 8

-- | Compose a known series S(u) with expansion E(h) = c₀ + u(h).
composeSeries :: (Double -> PuiseuxSeries) -> PuiseuxSeries -> PuiseuxSeries
composeSeries taylorAt e =
  let c0 = constantTerm e
      u  = removeTerm 0 e
      s  = taylorAt c0
  in evalSeriesAt u s

-- | Evaluate a series S = Σ aₙ·t^n by substituting t = u (another series)
evalSeriesAt :: PuiseuxSeries -> PuiseuxSeries -> PuiseuxSeries
evalSeriesAt u (PuiseuxSeries ts) =
  foldr addSeries zeroPuiseux
    [ scaleSeries (coeff t) (powSeries u (pExp t))
    | t <- ts
    ]

-- | Zero series
zeroPuiseux :: PuiseuxSeries
zeroPuiseux = PuiseuxSeries []

-- | Raise a series to a rational power
-- Integer powers only for now; fractional powers (Puiseux case) are TODO
powSeries :: PuiseuxSeries -> Rational -> PuiseuxSeries
powSeries _ 0 = PuiseuxSeries [PuiseuxTerm 0 1.0]
powSeries u n
  | n > 0 && denominator n == 1 =
      let k = fromIntegral (numerator n) :: Int
      in foldl mulSeries (PuiseuxSeries [PuiseuxTerm 0 1.0]) (replicate k u)
  | otherwise = PuiseuxSeries []  -- fractional powers: TODO

-- | Remove the term with a given exponent from a series
removeTerm :: Rational -> PuiseuxSeries -> PuiseuxSeries
removeTerm e (PuiseuxSeries ts) =
  PuiseuxSeries $ filter (\t -> pExp t /= e) ts

-- | Get the constant term (h^0 coefficient) of a series
constantTerm :: PuiseuxSeries -> Double
constantTerm (PuiseuxSeries []) = 0
constantTerm (PuiseuxSeries (t:_))
  | pExp t == 0 = coeff t
  | otherwise   = 0

-- | Truncate a series to n terms
truncateSeries :: Int -> PuiseuxSeries -> PuiseuxSeries
truncateSeries n (PuiseuxSeries ts) = PuiseuxSeries (take n ts)

-- | Invert a series: compute 1/s
-- Uses geometric series: 1/(a·h^α·(1+w)) = h^(-α)/a · Σ (-w)^n
invertSeries :: PuiseuxSeries -> PuiseuxSeries
invertSeries s =
  case leadingTermNZ (stripZeros s) of
    Nothing -> PuiseuxSeries []  -- 1/0, undefined
    Just lt ->
      let alpha = pExp lt
          a     = coeff lt
          w     = truncateSeries depth (normalizeW s lt alpha a)
          negw  = scaleSeries (-1) w
          -- geometric series: Σ (-w)^n
          geo   = geometricSeries negw
          scale = 1.0 / a
          shift = negate alpha
      in stripZeros $ truncateSeries depth
           (shiftExponents shift (scaleSeries scale geo))

-- | Geometric series 1/(1-u) = Σ u^n, here we pass (-w) so it computes 1/(1+w)
geometricSeries :: PuiseuxSeries -> PuiseuxSeries
geometricSeries u =
  let upows = take (depth+1) $ iterate (truncateSeries depth . mulSeries u)
                                       (PuiseuxSeries [PuiseuxTerm 0 1.0])
  in truncateSeries depth $ foldr addSeries zeroPuiseux upows

-- | Handle Pow case: f^g
expandPow :: Expr -> Expr -> Double -> PuiseuxSeries
expandPow f (Const r) x0     = expandPowR f (toRational r) x0
expandPow f (Neg (Const r)) x0 = expandPowR f (toRational (-r)) x0
expandPow _ _ _              = PuiseuxSeries []  -- TODO: symbolic exponents

-- | Expand f^r where r is a rational number
expandPowR :: Expr -> Rational -> Double -> PuiseuxSeries
expandPowR f r x0 =
  let s = stripZeros $ truncateSeries depth (expand f x0)
  in case leadingTermNZ s of
       Nothing -> PuiseuxSeries []
       Just lt ->
         let alpha = pExp lt
             a     = coeff lt
             w     = truncateSeries depth (normalizeW s lt alpha a)
             binom = binomialSeries r w
             scale = a ** fromRational r
             shift = alpha * r
         in stripZeros $ truncateSeries depth (shiftExponents shift (scaleSeries scale binom))

-- | Compute w = s/(a*h^alpha) - 1
normalizeW :: PuiseuxSeries -> PuiseuxTerm -> Rational -> Double -> PuiseuxSeries
normalizeW (PuiseuxSeries ts) _lt alpha a =
  let shifted = [ PuiseuxTerm (pExp t - alpha) (coeff t / a) | t <- ts ]
  in removeTerm 0 (PuiseuxSeries shifted)

-- | Binomial series (1+w)^r = Σ C(r,n) * w^n
binomialSeries :: Rational -> PuiseuxSeries -> PuiseuxSeries
binomialSeries r w =
  let bcs   = binomCoeffs r depth
      wpows = take (depth+1) $ iterate (truncateSeries depth . mulSeries w)
                                       (PuiseuxSeries [PuiseuxTerm 0 1.0])
  in truncateSeries depth $ foldr addSeries zeroPuiseux
       [ scaleSeries c wp
       | (c, wp) <- zip bcs wpows
       ]

-- | Generalized binomial coefficients C(r,n) for rational r
binomCoeffs :: Rational -> Int -> [Double]
binomCoeffs r n = take (n+1) $ scanl step 1.0 [0..]
  where
    step acc k = acc * (fromRational r - fromIntegral k) / fromIntegral (k+1)

-- | Shift all exponents in a series by a rational amount
shiftExponents :: Rational -> PuiseuxSeries -> PuiseuxSeries
shiftExponents delta (PuiseuxSeries ts) =
  PuiseuxSeries [ PuiseuxTerm (pExp t + delta) (coeff t) | t <- ts ]

-- | Handle Abs case — stub for now
expandAbs :: Expr -> Double -> PuiseuxSeries
expandAbs _ _ = PuiseuxSeries []  -- TODO

-- | Taylor series of sin around x₀
sinTaylor :: Double -> PuiseuxSeries
sinTaylor x0 = PuiseuxSeries $ take depth
  [ PuiseuxTerm (fromIntegral n) (sinCoeff n)
  | n <- [0..] ]
  where
    s = sin x0
    c = cos x0
    signs = cycle [s, c, -s, -c]
    facts = scanl (*) 1 [1..]
    sinCoeff n = (signs !! n) / fromIntegral (facts !! n)

-- | Taylor series of cos around x₀
cosTaylor :: Double -> PuiseuxSeries
cosTaylor x0 = PuiseuxSeries $ take depth
  [ PuiseuxTerm (fromIntegral n) (cosCoeff n)
  | n <- [0..] ]
  where
    s = sin x0
    c = cos x0
    signs = cycle [c, -s, -c, s]
    facts = scanl (*) 1 [1..]
    cosCoeff n = (signs !! n) / fromIntegral (facts !! n)

-- | Taylor series of exp around x₀: e^x₀ · Σ tⁿ/n!
expTaylor :: Double -> PuiseuxSeries
expTaylor x0 = PuiseuxSeries $ take depth
  [ PuiseuxTerm (fromIntegral n) (exp x0 / fromIntegral f)
  | (n, f) <- zip [0..] facts ]
  where
    facts = scanl (*) 1 [1..] :: [Integer]

-- | Taylor series of log around x₀
logTaylor :: Double -> PuiseuxSeries
logTaylor x0 = PuiseuxSeries $ take depth $
  PuiseuxTerm 0 (log x0) :
  [ PuiseuxTerm (fromIntegral n) (logCoeff n)
  | n <- [1..] ]
  where
    logCoeff n =
      let sign = if even n then -1 else 1
      in sign / (fromIntegral n * x0 ^ n)