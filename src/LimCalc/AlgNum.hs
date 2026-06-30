module LimCalc.AlgNum where

import Data.Ratio
import Data.Complex
import Data.List (minimumBy)
import Data.Ord (comparing)
import LimCalc.QPoly
import LimCalc.BivPoly

-- | Exact rational number (alias for clarity)
type Q = Rational

-- | Isolating rectangle in ℂ
data IsoRect = IsoRect
  { lowerLeft  :: GaussianQ
  , upperRight :: GaussianQ
  } deriving (Eq)

instance Show IsoRect where
  show (IsoRect ll ur) = "[" ++ show ll ++ ", " ++ show ur ++ "]"

-- | Check if a GaussianQ is inside an IsoRect
inRect :: GaussianQ -> IsoRect -> Bool
inRect (GQ r i) (IsoRect (GQ r1 i1) (GQ r2 i2)) =
  r1 <= r && r <= r2 && i1 <= i && i <= i2

-- | A complex algebraic number
data AlgNum = AlgNum
  { algMinPoly :: QPoly
  , algIsoRect :: IsoRect
  }

instance Show AlgNum where
  show (AlgNum p rect) = "AlgNum(" ++ show p ++ ", " ++ show rect ++ ")"

-- | Construct an algebraic number from a rational
fromQ :: Rational -> AlgNum
fromQ r = AlgNum
  { algMinPoly = QPoly [negate r, 1]
  , algIsoRect = IsoRect
      (GQ (r - 1) (negate 1))
      (GQ (r + 1) 1)
  }

-- | The algebraic number 0
algZero :: AlgNum
algZero = fromQ 0

-- | The algebraic number 1
algOne :: AlgNum
algOne = fromQ 1

-- | The imaginary unit i
--
-- Previously this rectangle was [-1, 1+1i], whose midpoint is
-- 0+0.5i, NOT i -- an asymmetric, miscentered rectangle present
-- since this constant was first written. This was invisible while
-- degree >= 2 AlgNum arithmetic was broken (root selection happened
-- to not depend on algI's own rectangle being accurate), but once
-- rectSum/rectMul's midpoint-based approximate-target estimation
-- started relying on it, it caused i+i to silently pick the wrong
-- root (0 instead of 2i) of an otherwise-correctly-computed
-- resultant. Centered properly here, matching the convention used
-- by fromQ and algSqrt.
algI :: AlgNum
algI = AlgNum
  { algMinPoly = QPoly [1, 0, 1]
  , algIsoRect = IsoRect
      (GQ (negate 1) 0)
      (GQ 1 2)
  }

-- | Square root of a positive rational
algSqrt :: Rational -> AlgNum
algSqrt r = AlgNum
  { algMinPoly = QPoly [negate r, 0, 1]
  , algIsoRect = IsoRect
      (GQ 0 (negate 1))
      (GQ (r + 1) 1)
  }

-- | Evaluate AlgNum's minimal polynomial at a GaussianQ point
algEval :: AlgNum -> GaussianQ -> GaussianQ
algEval (AlgNum p _) = qEval p

-- | Refine isolating rectangle by bisection until width < epsilon
refineRect :: AlgNum -> Rational -> AlgNum
refineRect an@(AlgNum p rect) eps
  | rectWidth rect < eps = an
  | otherwise =
      let mid     = midPoint rect
          upper   = IsoRect mid (upperRight rect)
          lower   = IsoRect (lowerLeft rect) mid
          newRect = if hasRoot p upper then upper else lower
      in refineRect (AlgNum p newRect) eps

-- | Width of an isolating rectangle
rectWidth :: IsoRect -> Rational
rectWidth (IsoRect (GQ r1 _) (GQ r2 _)) = abs (r2 - r1)

-- | Midpoint of an isolating rectangle
midPoint :: IsoRect -> GaussianQ
midPoint (IsoRect ll ur) = (ll + ur) / 2

-- | Check if a polynomial has a root in a rectangle
hasRoot :: QPoly -> IsoRect -> Bool
hasRoot p rect =
  let ll  = lowerLeft rect
      ur  = upperRight rect
      -- Check real axis interval [realQ ll, realQ ur]
      -- Evaluate polynomial at real endpoints
      rLeft  = realQ (qEval p (GQ (realQ ll) 0))
      rRight = realQ (qEval p (GQ (realQ ur) 0))
      rMid   = realQ (qEval p (GQ ((realQ ll + realQ ur) / 2) 0))
  in (rLeft * rMid <= 0) || (rMid * rRight <= 0) || abs rMid < 1e-10

-- | Add two algebraic numbers
algAdd :: AlgNum -> AlgNum -> AlgNum
algAdd a b =
  let resPoly = addResultantQ (algMinPoly a) (algMinPoly b)
      rect    = rectSum (algIsoRect a) (algIsoRect b)
  in refineToRoot (AlgNum resPoly rect)

-- | Negate an algebraic number
--
-- Negating coefficients at odd indices gives a polynomial with the
-- correct root, but does NOT preserve monic-ness (e.g. negating
-- algOne's [-1, 1] gives [-1, -1], leading coeff -1, not 1).
-- addResultantQ/mulResultantQ's linear-case shortcuts assume a monic
-- input (root = -c0/1) and silently mis-extract the root otherwise.
-- Rescaling by 1/leadingCoeff after negation keeps the result monic,
-- matching the invariant the rest of the module relies on.
algNeg :: AlgNum -> AlgNum
algNeg (AlgNum p rect) =
  let negP0 = QPoly [ if even i then c else negate c
                     | (i, c) <- zip [0 :: Int ..] (qPolyCoef p) ]
      negP  = qStrip negP0
      lc    = qLeadingCoeff negP
      monicNegP = if lc == 0 then negP else qScalePoly (1/lc) negP
      -- Negating a complex rectangle requires negating BOTH the
      -- real and imaginary extents (and swapping corners so
      -- lowerLeft/upperRight remain correctly ordered), not just the
      -- real part. The previous version only negated realQ, leaving
      -- imagQ untouched -- so negate algI returned a rectangle
      -- identical to algI's own, silently corrupting every
      -- computation that negates a non-real AlgNum (e.g. subtraction
      -- via algAdd a (algNeg b)).
      negRect = IsoRect
        (GQ (negate (realQ (upperRight rect))) (negate (imagQ (upperRight rect))))
        (GQ (negate (realQ (lowerLeft rect)))  (negate (imagQ (lowerLeft rect))))
  in AlgNum monicNegP negRect

-- | Multiply two algebraic numbers
algMul :: AlgNum -> AlgNum -> AlgNum
algMul a b =
  let resPoly = mulResultantQ (algMinPoly a) (algMinPoly b)
      rect    = rectMul (algIsoRect a) (algIsoRect b)
  in refineToRoot (AlgNum resPoly rect)

-- | Invert an algebraic number
algInv :: AlgNum -> AlgNum
algInv (AlgNum p rect) =
  let cs      = qPolyCoef p
      invP    = qStrip $ QPoly (reverse cs)
      -- make monic
      lc      = qLeadingCoeff invP
      monicP  = qScalePoly (1/lc) invP
      invRect = rectInv rect
  in refineToRoot (AlgNum monicP invRect)

-- | Sum of two isolating rectangles
rectSum :: IsoRect -> IsoRect -> IsoRect
rectSum (IsoRect ll1 ur1) (IsoRect ll2 ur2) =
  IsoRect (ll1 + ll2) (ur1 + ur2)

-- | Product of two isolating rectangles.
--
-- Previously this took a bounding box of the four corner products,
-- which doesn't reliably center on the TRUE product of the two
-- represented values -- complex multiplication doesn't preserve
-- axis-aligned rectangles, so the corner bounding box can be
-- significantly off-center from the actual product (observed: for
-- algI * algI, the corner bounding box centered near 0+0.5i, nowhere
-- near the true product -1, causing refineToRoot's Durand-Kerner
-- root selection to pick the wrong root of x^2-1).
--
-- The midpoint of each rectangle is itself an approximation of the
-- algebraic number it represents (see algToDouble/algImagDouble's
-- convention), so the product of the two midpoints is the correct
-- first-order estimate of where the true product lands. This is
-- used as the rectangle's center; the corner spread still
-- contributes to sizing the rectangle's width, but no longer
-- determines its center.
rectMul :: IsoRect -> IsoRect -> IsoRect
rectMul r1 r2 =
  let mid1 = midPoint r1
      mid2 = midPoint r2
      center = mid1 * mid2
      corners = [ lowerLeft r1 * lowerLeft r2
                , lowerLeft r1 * upperRight r2
                , upperRight r1 * lowerLeft r2
                , upperRight r1 * upperRight r2
                ]
      spreadR = (maximum (map realQ corners) - minimum (map realQ corners)) / 2
      spreadI = (maximum (map imagQ corners) - minimum (map imagQ corners)) / 2
      -- Guard against a degenerate (zero-width) spread, which would
      -- otherwise produce a zero-width rectangle that excludes the
      -- center itself due to rounding.
      halfWidth d = max d (1 % 2)
  in IsoRect
       (center - GQ (halfWidth spreadR) (halfWidth spreadI))
       (center + GQ (halfWidth spreadR) (halfWidth spreadI))

-- | Inverse of an isolating rectangle
rectInv :: IsoRect -> IsoRect
rectInv rect =
  let corners = [ lowerLeft rect
                , GQ (realQ (upperRight rect)) (imagQ (lowerLeft rect))
                , upperRight rect
                , GQ (realQ (lowerLeft rect)) (imagQ (upperRight rect))
                ]
      invCorners = map (\c -> GQ 1 0 / c) corners
      minR = minimum (map realQ invCorners)
      maxR = maximum (map realQ invCorners)
      minI = minimum (map imagQ invCorners)
      maxI = maximum (map imagQ invCorners)
  in IsoRect (GQ minR minI) (GQ maxR maxI)

-- | Refine until rectangle is tight.
--
-- If the minimal polynomial is degree <= 1 after stripping, the
-- "algebraic number" is actually exactly rational: c0 + c1*x = 0
-- has the exact root x = -c0/c1, with zero imaginary part since
-- the coefficients are real rationals. In that case we construct
-- a clean isolating rectangle directly around the exact root,
-- exactly mirroring what 'fromQ' does, instead of running the
-- (real-axis-only) bisection in 'refineRect', which has no business
-- being invoked when there is nothing left to isolate.
--
-- For degree >= 2, we use Durand-Kerner to find all complex roots
-- numerically, then pick whichever root is closest to the rectangle
-- we already have (which came from rectSum/rectMul -- an approximate
-- but reliable estimate of which root we actually want, since
-- algAdd/algMul know roughly where the true sum/product should be).
refineToRoot :: AlgNum -> AlgNum
refineToRoot an@(AlgNum p _) =
  let p' = qStrip p
  in case qDegree p' of
       0 -> an  -- degenerate: no roots; leave rectangle as-is
       1 -> let c0 = qPolyCoef p' !! 0
                c1 = qPolyCoef p' !! 1
                root = negate c0 / c1
            in AlgNum p' (exactRect root)
       _ -> let approxTarget = complexMidpoint (algIsoRect an)
                roots        = durandKerner p'
            in case roots of
                 [] -> an  -- shouldn't happen for degree >= 2; fail safe
                 _  -> let chosen = nearestTo approxTarget roots
                       in AlgNum p' (rectAroundComplex chosen)

-- | Construct a clean isolating rectangle around an exact rational
-- root, matching the convention used by 'fromQ'.
exactRect :: Rational -> IsoRect
exactRect root = IsoRect
  (GQ (root - 1) (negate 1))
  (GQ (root + 1) 1)

------------------------------------------------------------------------
-- Complex root-finding (Durand-Kerner / Weierstrass iteration)
------------------------------------------------------------------------

-- | Find all complex roots of a QPoly via Durand-Kerner iteration.
-- Converts the (exact, rational) coefficients to Complex Double for
-- the iteration, since Durand-Kerner is inherently numerical. Returns
-- an empty list only for degenerate input (degree < 1).
durandKerner :: QPoly -> [Complex Double]
durandKerner p
  | n < 1     = []
  | otherwise = iterateDK cs (initialGuesses n) durandKernerMaxIters
  where
    n  = qDegree p
    cs = map fromRational (qPolyCoef p) :: [Double]
       -- ascending degree, same convention as QPoly itself

-- | Hard cap on Durand-Kerner iterations, guarding against
-- pathological non-convergence (e.g. very closely clustered or
-- repeated roots).
durandKernerMaxIters :: Int
durandKernerMaxIters = 200

-- | Convergence tolerance: stop early once every root's update step
-- is smaller than this.
durandKernerTol :: Double
durandKernerTol = 1e-12

-- | Initial guesses spread around a circle in the complex plane,
-- the standard starting point for Durand-Kerner. Using a radius
-- slightly off 1 and irrational-ish angle offsets helps avoid
-- symmetric polynomials producing degenerate initial configurations.
initialGuesses :: Int -> [Complex Double]
initialGuesses n =
  [ mkPolar 1.13 (2 * pi * fromIntegral k / fromIntegral n + 0.3)
  | k <- [0 .. n - 1]
  ]

-- | Evaluate a polynomial (ascending-degree coefficient list) at a
-- complex point via Horner's method.
evalComplexPoly :: [Double] -> Complex Double -> Complex Double
evalComplexPoly cs z = foldr (\c acc -> (c :+ 0) + z * acc) (0 :+ 0) cs

-- | One Durand-Kerner iteration step for all roots simultaneously.
--
-- The correct update is z_i <- z_i - p(z_i) / (a_n * prod_{j/=i} (z_i - z_j)),
-- where a_n is p's LEADING coefficient. Previously this divided only
-- by the product of differences, omitting a_n entirely -- correct
-- only for monic polynomials. The resultant polynomials this is
-- actually called on (from addResultantBiv/mulResultantBiv) are not
-- normalized to be monic, so omitting a_n caused the iteration to
-- diverge (observed: roots drifting to magnitude ~10^44 instead of
-- converging) rather than merely being off by a constant factor.
dkStep :: [Double] -> [Complex Double] -> [Complex Double]
dkStep cs zs =
  [ z - evalComplexPoly cs z / (leading * denom z (filter (/= z) zs'))
  | (z, zs') <- zip zs (repeat zs)
  ]
  where
    leading = (last cs) :+ 0
    denom z others = product [ z - w | w <- others ]

-- | Iterate Durand-Kerner until convergence or the iteration cap is
-- hit, whichever comes first.
iterateDK :: [Double] -> [Complex Double] -> Int -> [Complex Double]
iterateDK _  zs 0 = zs
iterateDK cs zs n =
  let zs' = dkStep cs zs
      delta = maximum [ magnitude (a - b) | (a, b) <- zip zs zs' ]
  in if delta < durandKernerTol
       then zs'
       else iterateDK cs zs' (n - 1)

-- | Midpoint of an isolating rectangle as a Complex Double, used to
-- pick the "intended" root out of all roots Durand-Kerner finds.
complexMidpoint :: IsoRect -> Complex Double
complexMidpoint rect =
  let GQ r i = midPoint rect
  in fromRational r :+ fromRational i

-- | Pick whichever candidate root is closest to a target point.
nearestTo :: Complex Double -> [Complex Double] -> Complex Double
nearestTo target = minimumBy (comparing (magnitude . subtract target))

-- | Build a fresh isolating rectangle around an approximate complex
-- root found numerically. Since Durand-Kerner gives only a Double
-- approximation (not an exact algebraic isolating interval), this
-- rectangle is sized to comfortably contain the true root assuming
-- the numerical approximation is accurate to its converged
-- tolerance, with margin for safety. This is consistent in spirit
-- with the rest of the module's already-approximate transcendental
-- functions (algSin, algExp, etc., which are likewise Double-based).
rectAroundComplex :: Complex Double -> IsoRect
rectAroundComplex z =
  let re = toRational (realPart z)
      im = toRational (imagPart z)
      margin = 1 % 1000000  -- generous relative to durandKernerTol
  in IsoRect
       (GQ (re - margin) (im - margin))
       (GQ (re + margin) (im + margin))

-- | Num instance for AlgNum
--
-- abs was previously defined as identity (abs an = an), which is
-- wrong: 'abs' must return the magnitude, not the value itself.
-- This silently broke every stripZeros/leadingTermNZ/isAlgZero-style
-- comparison for negative coefficients, since comparing 'abs x > eps'
-- via the Ord instance (which compares algToDouble values) would
-- compare the *signed* double against eps -- so any negative
-- coefficient with magnitude bigger than eps was incorrectly treated
-- as smaller than eps and stripped as if it were zero.
instance Num AlgNum where
  (+)         = algAdd
  (*)         = algMul
  negate      = algNeg
  abs an      = if algToDouble an < 0 then algNeg an else an
  signum _    = algOne
  fromInteger = fromQ . fromInteger

-- | Fractional instance for AlgNum
instance Fractional AlgNum where
  recip        = algInv
  fromRational = fromQ

-- | Approximate transcendental functions at AlgNum
-- These return approximate AlgNum values using Double arithmetic
algSin :: AlgNum -> AlgNum
algSin a = fromRational (toRational (sin (algToDouble a)))

algCos :: AlgNum -> AlgNum
algCos a = fromRational (toRational (cos (algToDouble a)))

algExp :: AlgNum -> AlgNum
algExp a = fromRational (toRational (exp (algToDouble a)))

algLog :: AlgNum -> AlgNum
algLog a = fromRational (toRational (log (algToDouble a)))

-- | Convert AlgNum to approximate Double
-- Uses midpoint of isolating rectangle
algToDouble :: AlgNum -> Double
algToDouble (AlgNum _ rect) =
  let mid = midPoint rect
  in fromRational (realQ mid)

-- | Approximate imaginary part of an AlgNum, as a Double.
-- Mirrors algToDouble's convention (midpoint of isolating rectangle)
-- but projects the imaginary component instead of the real one.
-- Needed to detect non-real leading coefficients, e.g. when deciding
-- whether Abs has a sensible real-valued local expansion.
algImagDouble :: AlgNum -> Double
algImagDouble (AlgNum _ rect) =
  let mid = midPoint rect
  in fromRational (imagQ mid)

-- | Ord instance (approximate via Double)
instance Ord AlgNum where
  compare a b = compare (algToDouble a) (algToDouble b)

-- | Real instance (approximate via Double)
instance Real AlgNum where
  toRational = toRational . algToDouble

-- | RealFrac instance (approximate)
instance RealFrac AlgNum where
  properFraction a =
    let d = algToDouble a
        n = truncate d :: Integer
    in (fromIntegral n, fromQ (toRational (d - fromIntegral n)))
  truncate = truncate . algToDouble
  round    = round    . algToDouble
  ceiling  = ceiling  . algToDouble
  floor    = floor    . algToDouble

-- | Floating instance (approximate via Double)
instance Floating AlgNum where
  pi      = fromQ (toRational (pi :: Double))
  exp     = algExp
  log     = algLog
  sin     = algSin
  cos     = algCos
  sqrt a  = fromQ (toRational (sqrt (algToDouble a)))
  (**)  a b = fromQ (toRational (algToDouble a ** algToDouble b))
  asin  a = fromQ (toRational (asin  (algToDouble a)))
  acos  a = fromQ (toRational (acos  (algToDouble a)))
  atan  a = fromQ (toRational (atan  (algToDouble a)))
  sinh  a = fromQ (toRational (sinh  (algToDouble a)))
  cosh  a = fromQ (toRational (cosh  (algToDouble a)))
  asinh a = fromQ (toRational (asinh (algToDouble a)))
  acosh a = fromQ (toRational (acosh (algToDouble a)))
  atanh a = fromQ (toRational (atanh (algToDouble a)))

-- | Check if AlgNum is zero
isAlgZero :: AlgNum -> Bool
isAlgZero a = abs (algToDouble a) < 1e-12

instance Eq AlgNum where
  a == b = abs (algToDouble a - algToDouble b) < 1e-12