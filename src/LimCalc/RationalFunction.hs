module LimCalc.RationalFunction where

import LimCalc.Poly

-- | A rational function p/q
-- Invariant: q is never zero, leading coeff of q is positive
data RatFun = RatFun
  { numerator   :: Poly
  , denominator :: Poly
  } deriving (Eq)

instance Show RatFun where
  show (RatFun p q)
    | degree q == 0 = show p
    | otherwise     = "(" ++ show p ++ ") / (" ++ show q ++ ")"

-- | Construct a rational function, reducing to lowest terms
ratFun :: Poly -> Poly -> RatFun
ratFun p q
  | degree q < 0 = error "Zero denominator"
  | otherwise    =
      let g  = gcdPoly p q
          p' = quotPoly p g
          q' = quotPoly q g
          lc = leadingCoeff q'
      in RatFun (scalePoly (1/lc) p') (scalePoly (1/lc) q')

-- | Zero rational function
zeroRat :: String -> RatFun
zeroRat x = RatFun (zeroPoly x) (onePoly x)

-- | Constant rational function
constRat :: String -> Double -> RatFun
constRat x c = RatFun (constPoly x c) (onePoly x)

-- | Add two rational functions
addRat :: RatFun -> RatFun -> RatFun
addRat (RatFun p1 q1) (RatFun p2 q2) =
  ratFun (addPoly (mulPoly p1 q2) (mulPoly p2 q1)) (mulPoly q1 q2)

-- | Negate a rational function
negRat :: RatFun -> RatFun
negRat (RatFun p q) = RatFun (negPoly p) q

-- | Subtract two rational functions
subRat :: RatFun -> RatFun -> RatFun
subRat r s = addRat r (negRat s)

-- | Multiply two rational functions
mulRat :: RatFun -> RatFun -> RatFun
mulRat (RatFun p1 q1) (RatFun p2 q2) =
  ratFun (mulPoly p1 p2) (mulPoly q1 q2)

-- | Invert a rational function
invRat :: RatFun -> RatFun
invRat (RatFun p q) = ratFun q p

-- | Divide two rational functions
divRat :: RatFun -> RatFun -> RatFun
divRat r s = mulRat r (invRat s)

-- | Differentiate a rational function
-- D(p/q) = (p'q - pq') / q²
diffRat :: RatFun -> RatFun
diffRat (RatFun p q) =
  ratFun (subPoly (mulPoly (diffPoly p) q) (mulPoly p (diffPoly q)))
         (mulPoly q q)

-- | Polynomial part and proper fraction part
ratProperFraction :: RatFun -> (Poly, RatFun)
ratProperFraction (RatFun p q)
  | degree p < degree q = (zeroPoly (polyVar p), RatFun p q)
  | otherwise           =
      let (quot, rem) = divModPoly p q
      in (quot, RatFun rem q)

-- | Hermite reduction
-- Given p/q, returns (g, h) such that:
--   ∫ p/q = g + ∫ h
-- where g is a rational function and h has a squarefree denominator.
-- Uses the subresultant approach to avoid factoring.
hermiteReduce :: RatFun -> (RatFun, RatFun)
hermiteReduce rf =
  let (polyPart, proper) = ratProperFraction rf
      polyRat = RatFun polyPart (onePoly (polyVar polyPart))
      (g, h)  = hermiteReduceProper proper
  in (addRat polyRat g, h)

-- | Hermite reduction for proper rational functions
-- Core algorithm: repeatedly extract square factors from denominator
hermiteReduceProper :: RatFun -> (RatFun, RatFun)
hermiteReduceProper (RatFun p q) = go p q (zeroRat (polyVar p))
  where
    go a d acc
      | degree (gcdPoly d (diffPoly d)) == 0 =
          -- d is squarefree, we're done
          (acc, RatFun a d)
      | otherwise =
          -- split d into squarefree and non-squarefree parts
          let d'  = diffPoly d
              v   = gcdPoly d d'        -- repeated factor part
              u   = quotPoly d v        -- squarefree part
              -- solve: u*b + (v'/v)*c = a/v using extended Euclidean
              -- Actually use: find b,c such that u*b - v'/(deg)*c = a
              -- Simplified: use the Hermite identity directly
              (b, c) = hermiteStep a u v
              -- g contribution: -c/v
              gContrib = RatFun (negPoly c) v
              -- remaining numerator over v
              a' = addPoly (mulPoly u b) (mulPoly (diffPoly c) (quotPoly v u))
          in go a' v (addRat acc gContrib)

-- | One step of Hermite reduction
-- Given a, u, v where d = u*v and gcd(u,v)=1:
-- Find b, c such that a = u*b + D(v)*c (mod v)
-- Returns (b, c)
hermiteStep :: Poly -> Poly -> Poly -> (Poly, Poly)
hermiteStep a u v =
  let dv    = diffPoly v
      -- Extended GCD: find s,t such that u*s + dv*t = gcd(u,dv)
      (g, s, t) = extGCD u dv
      -- Scale to match a
      -- a = g * (a/g), so b = s*(a/g), c = t*(a/g)
      ag    = quotPoly a g
      b     = remPoly (mulPoly s ag) v
      c     = remPoly (mulPoly t ag) v
  in (b, c)

-- | Extended GCD: returns (g, s, t) such that s*a + t*b = g = gcd(a,b)
extGCD :: Poly -> Poly -> (Poly, Poly, Poly)
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

-- | Partial fraction decomposition
-- Given p/q where q = q1 * q2 * ... (squarefree factors),
-- returns list of (numerator, denominator) pairs
partialFractions :: RatFun -> [RatFun]
partialFractions rf@(RatFun p q) =
  let factors = squarefree q
  in case factors of
       []  -> [rf]
       [_] -> [rf]  -- already irreducible
       _   -> splitFractions p q factors

-- | Split into partial fractions given squarefree factors
splitFractions :: Poly -> Poly -> [(Poly, Int)] -> [RatFun]
splitFractions p q factors =
  -- For now: return unsplit (stub for full implementation)
  [RatFun p q]