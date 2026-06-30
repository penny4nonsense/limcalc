module LimCalc.RationalFunction where

import LimCalc.Poly

-- | A rational function p/q over coefficient type a
data RatFun a = RatFun
  { numerator   :: Poly a
  , denominator :: Poly a
  } deriving (Eq)

instance (Show a, Num a, Eq a) => Show (RatFun a) where
  show (RatFun p q)
    | degree q == 0 = show p
    | otherwise     = "(" ++ show p ++ ") / (" ++ show q ++ ")"

-- | Construct a rational function, reducing to lowest terms
ratFun :: (Fractional a, Eq a) => Poly a -> Poly a -> RatFun a
ratFun p q
  | degree q < 0 = error "Zero denominator"
  | otherwise    =
      let g  = gcdPoly p q
          p' = quotPoly p g
          q' = quotPoly q g
          lc = leadingCoeff q'
      in RatFun (scalePoly (1/lc) p') (scalePoly (1/lc) q')

-- | Zero rational function
zeroRat :: (Num a, Eq a) => String -> RatFun a
zeroRat x = RatFun (zeroPoly x) (onePoly x)

-- | Constant rational function
constRat :: (Fractional a, Eq a) => String -> a -> RatFun a
constRat x c = RatFun (constPoly x c) (onePoly x)

-- | Add two rational functions
addRat :: (Fractional a, Eq a) => RatFun a -> RatFun a -> RatFun a
addRat (RatFun p1 q1) (RatFun p2 q2) =
  ratFun (addPoly (mulPoly p1 q2) (mulPoly p2 q1)) (mulPoly q1 q2)

-- | Negate a rational function
negRat :: (Num a, Eq a) => RatFun a -> RatFun a
negRat (RatFun p q) = RatFun (negPoly p) q

-- | Subtract two rational functions
subRat :: (Fractional a, Eq a) => RatFun a -> RatFun a -> RatFun a
subRat r s = addRat r (negRat s)

-- | Multiply two rational functions
mulRat :: (Fractional a, Eq a) => RatFun a -> RatFun a -> RatFun a
mulRat (RatFun p1 q1) (RatFun p2 q2) =
  ratFun (mulPoly p1 p2) (mulPoly q1 q2)

-- | Invert a rational function
invRat :: (Fractional a, Eq a) => RatFun a -> RatFun a
invRat (RatFun p q) = ratFun q p

-- | Divide two rational functions
divRat :: (Fractional a, Eq a) => RatFun a -> RatFun a -> RatFun a
divRat r s = mulRat r (invRat s)

-- | Differentiate a rational function
diffRat :: (Fractional a, Eq a) => RatFun a -> RatFun a
diffRat (RatFun p q) =
  ratFun (subPoly (mulPoly (diffPoly p) q) (mulPoly p (diffPoly q)))
         (mulPoly q q)

-- | Polynomial part and proper fraction part
ratProperFraction :: (Fractional a, Eq a) => RatFun a -> (Poly a, RatFun a)
ratProperFraction (RatFun p q)
  | degree p < degree q = (zeroPoly (polyVar p), RatFun p q)
  | otherwise           =
      let (polyQ, polyR) = divModPoly p q
      in (polyQ, RatFun polyR q)

-- | Hermite reduction
hermiteReduce :: (Fractional a, Eq a) => RatFun a -> (RatFun a, RatFun a)
hermiteReduce rf =
  let (polyPart, proper) = ratProperFraction rf
      polyRat = RatFun polyPart (onePoly (polyVar polyPart))
      (g, h)  = hermiteReduceProper proper
  in (addRat polyRat g, h)

-- | Hermite reduction for proper rational functions
hermiteReduceProper :: (Fractional a, Eq a) => RatFun a -> (RatFun a, RatFun a)
hermiteReduceProper (RatFun p q) = go p q (zeroRat (polyVar p))
  where
    go a d acc
      | degree (gcdPoly d (diffPoly d)) == 0 =
          (acc, RatFun a d)
      | otherwise =
          let d'         = diffPoly d
              v          = gcdPoly d d'
              u          = quotPoly d v
              (b, c)     = hermiteStep a u v
              gContrib   = RatFun (negPoly c) v
              a'         = addPoly (mulPoly u b)
                             (mulPoly (diffPoly c) (quotPoly v u))
          in go a' v (addRat acc gContrib)

-- | One step of Hermite reduction
hermiteStep :: (Fractional a, Eq a) => Poly a -> Poly a -> Poly a -> (Poly a, Poly a)
hermiteStep a u v =
  let dv      = diffPoly v
      (g, s, t) = extGCD u dv
      ag      = quotPoly a g
      b       = remPoly (mulPoly s ag) v
      c       = remPoly (mulPoly t ag) v
  in (b, c)

-- | Extended GCD
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

-- | Partial fraction decomposition.
-- Implemented in LimCalc.Risch.Primitive (which has access to
-- Rothstein-Trager). This stub is kept as documentation only.
-- Use Risch.Primitive.partialFractions instead.
partialFractionsStub :: RatFun a -> [RatFun a]
partialFractionsStub rf = [rf]