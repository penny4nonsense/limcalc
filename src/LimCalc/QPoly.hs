module LimCalc.QPoly where

import Data.Ratio

-- | Gaussian rational: a + bi where a, b ∈ ℚ
data GaussianQ = GQ
  { realQ :: Rational
  , imagQ :: Rational
  } deriving (Eq)

instance Show GaussianQ where
  show (GQ r 0) = showRat r
  show (GQ 0 i) = showRat i ++ "i"
  show (GQ r i)
    | i < 0     = showRat r ++ " - " ++ showRat (abs i) ++ "i"
    | otherwise = showRat r ++ " + " ++ showRat i ++ "i"

showRat :: Rational -> String
showRat r
  | denominator r == 1 = show (numerator r)
  | otherwise          = show (numerator r) ++ "/" ++ show (denominator r)

instance Num GaussianQ where
  (GQ a b) + (GQ c d) = GQ (a+c) (b+d)
  (GQ a b) * (GQ c d) = GQ (a*c - b*d) (a*d + b*c)
  negate (GQ a b)     = GQ (negate a) (negate b)
  abs gq              = GQ (gaussianNorm gq) 0
  signum _            = GQ 1 0
  fromInteger n       = GQ (fromInteger n) 0

instance Fractional GaussianQ where
  (GQ a b) / (GQ c d) =
    let denom = c*c + d*d
    in GQ ((a*c + b*d) / denom) ((b*c - a*d) / denom)
  fromRational r = GQ r 0

-- | Squared norm of a Gaussian rational
gaussianNorm :: GaussianQ -> Rational
gaussianNorm (GQ a b) = a*a + b*b

-- | Conjugate of a Gaussian rational
conjugateQ :: GaussianQ -> GaussianQ
conjugateQ (GQ a b) = GQ a (negate b)

-- | Polynomial over ℚ
-- Coefficients in ascending degree order
newtype QPoly = QPoly { qPolyCoef :: [Rational] }
  deriving (Eq)

instance Show QPoly where
  show (QPoly [])  = "0"
  show (QPoly cs)  = concatMap showQTerm (reverse $ zip [0..] cs)
    where
      showQTerm (0, c) = showRat c
      showQTerm (1, c) = showRat c ++ "x"
      showQTerm (n, c) = showRat c ++ "x^" ++ show (n :: Int)

-- | Degree of a QPoly
qDegree :: QPoly -> Int
qDegree (QPoly []) = -1
qDegree (QPoly cs) = length cs - 1

-- | Leading coefficient of a QPoly
qLeadingCoeff :: QPoly -> Rational
qLeadingCoeff (QPoly []) = 0
qLeadingCoeff (QPoly cs) = last cs

-- | Evaluate a QPoly at a GaussianQ point
qEval :: QPoly -> GaussianQ -> GaussianQ
qEval (QPoly cs) x =
  foldr (\c acc -> fromRational c + x * acc) 0 cs

-- | Strip trailing zeros from QPoly
qStrip :: QPoly -> QPoly
qStrip (QPoly cs) = QPoly (reverse $ dropWhile (== 0) $ reverse cs)

-- | Add two QPolys
qAddPoly :: QPoly -> QPoly -> QPoly
qAddPoly (QPoly cs1) (QPoly cs2) =
  qStrip $ QPoly $ addCoefs cs1 cs2
  where
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (x:xs) (y:ys) = (x+y) : addCoefs xs ys

-- | Multiply two QPolys
qMulPoly :: QPoly -> QPoly -> QPoly
qMulPoly (QPoly []) _            = QPoly []
qMulPoly _ (QPoly [])            = QPoly []
qMulPoly (QPoly cs1) (QPoly cs2) =
  qStrip $ QPoly $ mulCoefs cs1 cs2
  where
    mulCoefs [] _      = []
    mulCoefs (c:cs) ys =
      addCoefs (map (*c) ys) (0 : mulCoefs cs ys)
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (x:xs) (y:ys) = (x+y) : addCoefs xs ys

-- | Scale a QPoly by a rational
qScalePoly :: Rational -> QPoly -> QPoly
qScalePoly 0 _          = QPoly []
qScalePoly c (QPoly cs) = qStrip $ QPoly (map (*c) cs)

-- | Negate a QPoly
qNegPoly :: QPoly -> QPoly
qNegPoly (QPoly cs) = QPoly (map negate cs)

-- | Subtract two QPolys
qSubPoly :: QPoly -> QPoly -> QPoly
qSubPoly p q = qAddPoly p (qNegPoly q)

-- | Differentiate a QPoly
qDiffPoly :: QPoly -> QPoly
qDiffPoly (QPoly []) = QPoly []
qDiffPoly (QPoly cs) =
  qStrip $ QPoly
    [ fromIntegral n * c
    | (n, c) <- zip [1 :: Int ..] (drop 1 cs)
    ]

-- | Power of a QPoly
qPow :: QPoly -> Int -> QPoly
qPow _ 0 = QPoly [1]
qPow p n = qMulPoly p (qPow p (n-1))

-- | Num instance for QPoly, wrapping the existing qAddPoly/qMulPoly/
-- etc. This lets QPoly be used as a coefficient ring with generic
-- polynomial machinery (e.g. Poly.hs's Poly a / resultant), rather
-- than requiring a parallel, hand-rolled bivariate implementation.
instance Num QPoly where
  (+)         = qAddPoly
  (*)         = qMulPoly
  negate      = qNegPoly
  abs p       = p  -- magnitude isn't meaningful for a polynomial;
                    -- only used incidentally by generic code paths
                    -- that don't actually rely on it for QPoly
  signum _    = QPoly [1]
  fromInteger n = QPoly [fromInteger n]

-- | Fractional instance for QPoly. Division here is REAL polynomial
-- division (via qDivModPoly), discarding the remainder -- this is
-- only exact when the division is known to be exact by construction
-- (as in the subresultant/resultant algorithms that motivated this
-- instance). Using '/' on QPolys where the division isn't exact will
-- silently drop a remainder; prefer qDivModPoly directly when that
-- matters.
instance Fractional QPoly where
  p / q        = qQuotPoly p q
  fromRational r = QPoly [r]

-- | Euclidean division of QPolys: (quotient, remainder).
--
-- Mirrors Poly.divModPoly exactly. Needed because BivPoly's
-- subresultant algorithm divides by genuinely non-constant QPoly
-- leading coefficients (since BivPoly's own coefficients are QPolys
-- in x), which the previous qDivPoly stub silently failed to handle
-- -- it only worked for two degree-0 QPolys and otherwise returned
-- its numerator unchanged, corrupting the subresultant PRS recursion
-- and causing it to fail to terminate properly.
qDivModPoly :: QPoly -> QPoly -> (QPoly, QPoly)
qDivModPoly p q
  | qDegree q < 0          = error "qDivModPoly: division by zero QPoly"
  | qDegree p < qDegree q  = (QPoly [], p)
  | otherwise              = go p (QPoly [])
  where
    lc = qLeadingCoeff q
    dq = qDegree q
    go r acc
      | qDegree r < dq = (acc, r)
      | otherwise =
          let scale = qLeadingCoeff r / lc
              deg   = qDegree r - dq
              term  = QPoly (replicate deg 0 ++ [scale])
              r'    = qStrip $ qSubPoly r (qMulPoly term q)
          in go r' (qAddPoly acc term)

-- | Quotient only.
qQuotPoly :: QPoly -> QPoly -> QPoly
qQuotPoly p q = fst (qDivModPoly p q)

-- | Remainder only.
qRemPoly :: QPoly -> QPoly -> QPoly
qRemPoly p q = snd (qDivModPoly p q)