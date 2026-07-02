-- | Polynomials over ℚ and Gaussian rationals.
--
-- This module provides two types:
--
-- * 'GaussianQ' — exact Gaussian rationals @a + bi@ where @a, b ∈ ℚ@,
--   used as evaluation points in 'LimCalc.Algebra.AlgNum' and as the
--   coefficient type for isolating rectangles.
--
-- * 'QPoly' — univariate polynomials over ℚ (ascending-degree
--   coefficient lists), used as minimal polynomials in 'LimCalc.Algebra.AlgNum'
--   and as the coefficient ring for bivariate resultant computation
--   in 'LimCalc.BivPoly'.
--
-- = QPoly as a coefficient ring
--
-- 'QPoly' has 'Num' and 'Fractional' instances so that it can be
-- used as the coefficient type of 'LimCalc.Poly.Poly', enabling
-- bivariate polynomial arithmetic via the generic univariate
-- machinery. The 'Fractional' instance performs exact polynomial
-- division (via 'qDivModPoly'), which is only valid when the
-- division is known to be exact by construction; see the instance
-- documentation for the caveat.
module LimCalc.Algebra.QPoly
  ( -- * Gaussian rationals
    GaussianQ (..)
  , gaussianNorm
  , conjugateQ
  , showRat
    -- * Polynomials over ℚ
  , QPoly (..)
    -- ** Properties
  , qDegree
  , qLeadingCoeff
    -- ** Normalisation
  , qStrip
    -- ** Arithmetic
  , qAddPoly
  , qNegPoly
  , qSubPoly
  , qMulPoly
  , qScalePoly
  , qPow
    -- ** Calculus
  , qDiffPoly
  , qEval
    -- ** Division
  , qDivModPoly
  , qQuotPoly
  , qRemPoly
  ) where

import Data.Ratio

-- | An exact Gaussian rational @a + bi@ where @a, b ∈ ℚ@.
--
-- Used as evaluation points for minimal polynomials in
-- 'LimCalc.Algebra.AlgNum.algEval', and as the corner type for isolating
-- rectangles ('LimCalc.AlgNum.IsoRect').
data GaussianQ = GQ
  { realQ :: Rational
    -- ^ Real part.
  , imagQ :: Rational
    -- ^ Imaginary part.
  } deriving (Eq)

instance Show GaussianQ where
  show (GQ r 0) = showRat r
  show (GQ 0 i) = showRat i ++ "i"
  show (GQ r i)
    | i < 0     = showRat r ++ " - " ++ showRat (abs i) ++ "i"
    | otherwise = showRat r ++ " + " ++ showRat i ++ "i"

-- | Display a rational number without redundant denominators.
showRat :: Rational -> String
showRat r
  | denominator r == 1 = show (numerator r)
  | otherwise          = show (numerator r) ++ "/" ++ show (denominator r)

-- | 'Num' instance for 'GaussianQ': standard Gaussian integer arithmetic
-- extended to rationals.
instance Num GaussianQ where
  (GQ a b) + (GQ c d) = GQ (a+c) (b+d)
  (GQ a b) * (GQ c d) = GQ (a*c - b*d) (a*d + b*c)
  negate (GQ a b)     = GQ (negate a) (negate b)
  abs gq              = GQ (gaussianNorm gq) 0
  signum _            = GQ 1 0
  fromInteger n       = GQ (fromInteger n) 0

-- | 'Fractional' instance for 'GaussianQ': division by conjugate
-- multiplication.
instance Fractional GaussianQ where
  (GQ a b) / (GQ c d) =
    let denom = c*c + d*d
    in GQ ((a*c + b*d) / denom) ((b*c - a*d) / denom)
  fromRational r = GQ r 0

-- | Squared norm @|a + bi|² = a² + b²@ of a Gaussian rational.
gaussianNorm :: GaussianQ -> Rational
gaussianNorm (GQ a b) = a*a + b*b

-- | Complex conjugate @a + bi → a − bi@.
conjugateQ :: GaussianQ -> GaussianQ
conjugateQ (GQ a b) = GQ a (negate b)

-- | A univariate polynomial over ℚ, with coefficients in ascending
-- degree order. The zero polynomial is @QPoly []@. The last element
-- of a non-empty coefficient list is nonzero (maintained by 'qStrip').
newtype QPoly = QPoly { qPolyCoef :: [Rational] }
  deriving (Eq)

instance Show QPoly where
  show (QPoly [])  = "0"
  show (QPoly cs)  = concatMap showQTerm (reverse $ zip [0..] cs)
    where
      showQTerm (0, c) = showRat c
      showQTerm (1, c) = showRat c ++ "x"
      showQTerm (n, c) = showRat c ++ "x^" ++ show (n :: Int)

-- | Degree of a 'QPoly'. Returns @−1@ for the zero polynomial.
qDegree :: QPoly -> Int
qDegree (QPoly []) = -1
qDegree (QPoly cs) = length cs - 1

-- | Leading coefficient of a 'QPoly'. Returns @0@ for the zero polynomial.
qLeadingCoeff :: QPoly -> Rational
qLeadingCoeff (QPoly []) = 0
qLeadingCoeff (QPoly cs) = last cs

-- | Evaluate a 'QPoly' at a 'GaussianQ' point via Horner's method.
qEval :: QPoly -> GaussianQ -> GaussianQ
qEval (QPoly cs) x =
  foldr (\c acc -> fromRational c + x * acc) 0 cs

-- | Remove trailing zero coefficients.
qStrip :: QPoly -> QPoly
qStrip (QPoly cs) = QPoly (reverse $ dropWhile (== 0) $ reverse cs)

-- | Add two 'QPoly' values.
qAddPoly :: QPoly -> QPoly -> QPoly
qAddPoly (QPoly cs1) (QPoly cs2) =
  qStrip $ QPoly $ addCoefs cs1 cs2
  where
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (x:xs) (y:ys) = (x+y) : addCoefs xs ys

-- | Multiply two 'QPoly' values (Cauchy product).
qMulPoly :: QPoly -> QPoly -> QPoly
qMulPoly (QPoly []) _             = QPoly []
qMulPoly _ (QPoly [])             = QPoly []
qMulPoly (QPoly cs1) (QPoly cs2) =
  qStrip $ QPoly $ mulCoefs cs1 cs2
  where
    mulCoefs [] _      = []
    mulCoefs (c:cs) ys =
      addCoefs (map (*c) ys) (0 : mulCoefs cs ys)
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (x:xs) (y:ys) = (x+y) : addCoefs xs ys

-- | Scale every coefficient by a rational constant.
qScalePoly :: Rational -> QPoly -> QPoly
qScalePoly 0 _          = QPoly []
qScalePoly c (QPoly cs) = qStrip $ QPoly (map (*c) cs)

-- | Negate a 'QPoly'.
qNegPoly :: QPoly -> QPoly
qNegPoly (QPoly cs) = QPoly (map negate cs)

-- | Subtract two 'QPoly' values.
qSubPoly :: QPoly -> QPoly -> QPoly
qSubPoly p q = qAddPoly p (qNegPoly q)

-- | Formal derivative of a 'QPoly'.
qDiffPoly :: QPoly -> QPoly
qDiffPoly (QPoly []) = QPoly []
qDiffPoly (QPoly cs) =
  qStrip $ QPoly
    [ fromIntegral n * c
    | (n, c) <- zip [1 :: Int ..] (drop 1 cs)
    ]

-- | Integer power of a 'QPoly'.
qPow :: QPoly -> Int -> QPoly
qPow _ 0 = QPoly [1]
qPow p n = qMulPoly p (qPow p (n-1))

-- | 'Num' instance for 'QPoly'.
--
-- Allows 'QPoly' to be used as the coefficient type of
-- 'LimCalc.Algebra.Poly.Poly', enabling bivariate polynomial arithmetic via
-- the generic univariate machinery in 'LimCalc.Algebra.BivPoly'. The 'abs'
-- method returns the polynomial unchanged, since magnitude is not
-- meaningful for a polynomial ring; it is only invoked incidentally
-- by generic code paths that do not actually use it for 'QPoly'.
instance Num QPoly where
  (+)           = qAddPoly
  (*)           = qMulPoly
  negate        = qNegPoly
  abs p         = p
  signum _      = QPoly [1]
  fromInteger n = QPoly [fromInteger n]

-- | 'Fractional' instance for 'QPoly'.
--
-- Division is exact polynomial quotient via 'qDivModPoly', discarding
-- the remainder. This is only correct when the division is exact by
-- construction — as in the subresultant and resultant algorithms in
-- 'LimCalc.Algebra.BivPoly' that motivated this instance. Using @\/@ where
-- the division is not exact will silently drop a nonzero remainder;
-- prefer 'qDivModPoly' directly in those cases.
instance Fractional QPoly where
  p / q          = qQuotPoly p q
  fromRational r = QPoly [r]

-- | Euclidean division of 'QPoly' values: @qDivModPoly p q = (quot, rem)@
-- where @p = quot * q + rem@ and @degree rem < degree q@.
--
-- Required because 'LimCalc.Algebra.BivPoly'\'s subresultant algorithm divides
-- by non-constant 'QPoly' leading coefficients. An earlier stub only
-- handled degree-0 divisors and returned the numerator unchanged
-- otherwise, silently corrupting the subresultant pseudo-remainder
-- sequence.
qDivModPoly :: QPoly -> QPoly -> (QPoly, QPoly)
qDivModPoly p q
  | qDegree q < 0         = error "qDivModPoly: division by zero QPoly"
  | qDegree p < qDegree q = (QPoly [], p)
  | otherwise             = go p (QPoly [])
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

-- | Quotient of 'QPoly' division.
qQuotPoly :: QPoly -> QPoly -> QPoly
qQuotPoly p q = fst (qDivModPoly p q)

-- | Remainder of 'QPoly' division.
qRemPoly :: QPoly -> QPoly -> QPoly
qRemPoly p q = snd (qDivModPoly p q)