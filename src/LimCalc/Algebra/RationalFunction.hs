-- | Rational functions over an arbitrary coefficient field.
--
-- A rational function @p(x) \/ q(x)@ is represented in reduced form:
-- numerator and denominator are divided by their GCD and the
-- denominator is made monic. All arithmetic operations maintain
-- this invariant via 'ratFun'.
--
-- This module provides the arithmetic and structural operations on
-- rational functions needed by the Risch integrator
-- ('LimCalc.Integration.Risch.Primitive',
-- 'LimCalc.Integration.Risch.Exponential') and the differential field
-- tower ('LimCalc.Differentiation.DiffField'). Hermite reduction
-- ('hermiteReduce') is the key structural operation: it decomposes a
-- rational function into a "rational part" (whose integral is
-- rational) and a "logarithmic part" (a proper fraction over a
-- squarefree denominator, whose integral involves logarithms).
module LimCalc.Algebra.RationalFunction
  ( -- * Type
    RatFun (..)
    -- * Constructors
  , ratFun
  , zeroRat
  , constRat
    -- * Arithmetic
  , addRat
  , negRat
  , subRat
  , mulRat
  , invRat
  , divRat
  , diffRat
    -- * Decomposition
  , ratProperFraction
  , hermiteReduce
  , hermiteReduceProper
  , hermiteStep
    -- * Extended GCD
  , extGCD
  ) where

import LimCalc.Algebra.Poly

-- | A rational function @p(x) \/ q(x)@ in reduced form.
--
-- Invariant: @gcd(numerator, denominator) = 1@ and
-- @leadingCoeff(denominator) = 1@ (monic denominator).
-- Maintained by 'ratFun'.
data RatFun a = RatFun
  { numerator   :: Poly a
    -- ^ The numerator polynomial.
  , denominator :: Poly a
    -- ^ The denominator polynomial (monic, coprime with numerator).
  } deriving (Eq)

instance (Show a, Num a, Eq a) => Show (RatFun a) where
  show (RatFun p q)
    | degree q == 0 = show p
    | otherwise     = "(" ++ show p ++ ") / (" ++ show q ++ ")"

-- | Construct a rational function in reduced form.
--
-- Divides @p@ and @q@ by their GCD, then rescales so the denominator
-- is monic. Raises an error if @q@ is the zero polynomial.
ratFun :: (Fractional a, Eq a) => Poly a -> Poly a -> RatFun a
ratFun p q
  | degree q < 0 = error "Zero denominator"
  | otherwise    =
      let g  = gcdPoly p q
          p' = quotPoly p g
          q' = quotPoly q g
          lc = leadingCoeff q'
      in RatFun (scalePoly (1/lc) p') (scalePoly (1/lc) q')

-- | The zero rational function @0 \/ 1@.
zeroRat :: (Num a, Eq a) => String -> RatFun a
zeroRat x = RatFun (zeroPoly x) (onePoly x)

-- | A constant rational function @c \/ 1@.
constRat :: (Fractional a, Eq a) => String -> a -> RatFun a
constRat x c = RatFun (constPoly x c) (onePoly x)

-- | Add two rational functions: @p₁\/q₁ + p₂\/q₂ = (p₁q₂ + p₂q₁)\/(q₁q₂)@,
-- reduced to lowest terms.
addRat :: (Fractional a, Eq a) => RatFun a -> RatFun a -> RatFun a
addRat (RatFun p1 q1) (RatFun p2 q2) =
  ratFun (addPoly (mulPoly p1 q2) (mulPoly p2 q1)) (mulPoly q1 q2)

-- | Negate a rational function: @−(p\/q) = (−p)\/q@.
negRat :: (Num a, Eq a) => RatFun a -> RatFun a
negRat (RatFun p q) = RatFun (negPoly p) q

-- | Subtract two rational functions.
subRat :: (Fractional a, Eq a) => RatFun a -> RatFun a -> RatFun a
subRat r s = addRat r (negRat s)

-- | Multiply two rational functions: @(p₁\/q₁) · (p₂\/q₂) = (p₁p₂)\/(q₁q₂)@,
-- reduced to lowest terms.
mulRat :: (Fractional a, Eq a) => RatFun a -> RatFun a -> RatFun a
mulRat (RatFun p1 q1) (RatFun p2 q2) =
  ratFun (mulPoly p1 p2) (mulPoly q1 q2)

-- | Invert a rational function: @(p\/q)⁻¹ = q\/p@.
invRat :: (Fractional a, Eq a) => RatFun a -> RatFun a
invRat (RatFun p q) = ratFun q p

-- | Divide two rational functions.
divRat :: (Fractional a, Eq a) => RatFun a -> RatFun a -> RatFun a
divRat r s = mulRat r (invRat s)

-- | Differentiate a rational function via the quotient rule:
-- @(p\/q)' = (p'q − pq') \/ q²@, reduced to lowest terms.
diffRat :: (Fractional a, Eq a) => RatFun a -> RatFun a
diffRat (RatFun p q) =
  ratFun (subPoly (mulPoly (diffPoly p) q) (mulPoly p (diffPoly q)))
         (mulPoly q q)

-- | Split a rational function into its polynomial part and proper
-- fraction part: @f = poly + proper@ where @deg(numerator(proper)) < deg(denominator(proper))@.
ratProperFraction :: (Fractional a, Eq a) => RatFun a -> (Poly a, RatFun a)
ratProperFraction (RatFun p q)
  | degree p < degree q = (zeroPoly (polyVar p), RatFun p q)
  | otherwise           =
      let (polyQ, polyR) = divModPoly p q
      in (polyQ, RatFun polyR q)

-- | Hermite reduction of a rational function.
--
-- Decomposes @f = g + h@ where:
--
-- * @g@ is a rational function whose integral is rational (no
--   logarithmic terms); its denominator is a product of squared
--   or higher factors.
-- * @h@ is a proper fraction over a squarefree denominator; its
--   integral involves only logarithms.
--
-- The Rothstein-Trager algorithm ('LimCalc.Risch.Primitive') is then
-- applied to @h@ to complete the integration.
hermiteReduce :: (Fractional a, Eq a) => RatFun a -> (RatFun a, RatFun a)
hermiteReduce rf =
  let (polyPart, proper) = ratProperFraction rf
      polyRat = RatFun polyPart (onePoly (polyVar polyPart))
      (g, h)  = hermiteReduceProper proper
  in (addRat polyRat g, h)

-- | Hermite reduction for a proper rational function (degree of
-- numerator strictly less than degree of denominator).
--
-- Iteratively extracts the rational part by processing repeated
-- factors of the denominator via 'hermiteStep'.
hermiteReduceProper :: (Fractional a, Eq a) => RatFun a -> (RatFun a, RatFun a)
hermiteReduceProper (RatFun p q) = go p q (zeroRat (polyVar p))
  where
    go a d acc
      | degree (gcdPoly d (diffPoly d)) == 0 =
          (acc, RatFun a d)
      | otherwise =
          let d'     = diffPoly d
              v      = gcdPoly d d'
              u      = quotPoly d v
              (b, c) = hermiteStep a u v
              gContrib = RatFun (negPoly c) v
              a'     = addPoly (mulPoly u b)
                         (mulPoly (diffPoly c) (quotPoly v u))
          in go a' v (addRat acc gContrib)

-- | One step of Hermite reduction.
--
-- Given numerator @a@, squarefree part @u = d \/ gcd(d, d')@, and
-- repeated part @v = gcd(d, d')@, computes @(b, c)@ such that
-- @a \/ (u · v) = b\/u + (c\/v)'@ (the derivative of @c\/v@).
hermiteStep :: (Fractional a, Eq a) => Poly a -> Poly a -> Poly a -> (Poly a, Poly a)
hermiteStep a u v =
  let dv        = diffPoly v
      (g, s, t) = extGCD u dv
      ag        = quotPoly a g
      b         = remPoly (mulPoly s ag) v
      c         = remPoly (mulPoly t ag) v
  in (b, c)

-- | Extended Euclidean algorithm for polynomials.
--
-- Returns @(g, s, t)@ such that @g = gcd(a, b) = s·a + t·b@
-- (Bézout identity), with @g@ monic.
extGCD :: (Fractional a, Eq a) => Poly a -> Poly a -> (Poly a, Poly a, Poly a)
extGCD a b
  | degree b < 0 =
      let lc = leadingCoeff a
      in ( scalePoly (1/lc) a
         , scalePoly (1/lc) (onePoly (polyVar a))
         , zeroPoly (polyVar a)
         )
  | otherwise =
      let (q, r)    = divModPoly a b
          (g, s, t) = extGCD b r
          s'        = subPoly t (mulPoly q s)
      in (g, t, s')