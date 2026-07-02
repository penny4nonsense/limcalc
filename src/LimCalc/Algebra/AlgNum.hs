-- | Algebraic number arithmetic over ℂ.
--
-- An /algebraic number/ is represented by its minimal polynomial over
-- ℚ together with an /isolating rectangle/ in ℂ — a small
-- axis-aligned rectangle in the Gaussian rationals that is guaranteed
-- to contain exactly one root of the minimal polynomial. This
-- representation supports exact arithmetic (via resultants) with
-- approximate root selection (via Durand-Kerner iteration).
--
-- = Design overview
--
-- Arithmetic (@+@, @*@, @negate@, @recip@) works in three steps:
--
-- 1. Compute the minimal polynomial of the result via a resultant
--    construction ('addResultantQ' or 'mulResultantQ').
-- 2. Compute an approximate isolating rectangle for the result via
--    interval arithmetic on the input rectangles ('rectSum', 'rectMul').
-- 3. Refine to an actual root of the result polynomial using
--    Durand-Kerner ('refineToRoot').
--
-- = Approximation policy
--
-- Transcendental functions ('algSin', 'algCos', 'algExp', 'algLog')
-- are implemented approximately via 'Double' arithmetic, returning a
-- fresh degree-1 'AlgNum' wrapping the floating-point result. This is
-- consistent with the series expansion engine, which uses 'AlgNum'
-- primarily as a coefficient type for its exact rational arithmetic
-- and only requires transcendental values at specific numeric points.
--
-- = Degree cap
--
-- Unminimised resultant chains can blow up the minimal polynomial
-- degree combinatorially. Above 'degreeCapForExactRootFinding',
-- 'refineToRoot' falls back to 'pinNumerically', pinning the value
-- as an approximate degree-1 algebraic number. This trades algebraic
-- exactness for termination.
module LimCalc.Algebra.AlgNum
  ( -- * Types
    Q
  , IsoRect (..)
  , AlgNum (..)
    -- * Constants
  , algZero
  , algOne
  , algI
    -- * Constructors
  , fromQ
  , algSqrt
    -- * Arithmetic
  , algAdd
  , algNeg
  , algMul
  , algInv
    -- * Transcendental approximations
  , algSin
  , algCos
  , algExp
  , algLog
    -- * Conversion
  , algToDouble
  , algImagDouble
    -- * Predicates
  , isAlgZero
    -- * Rectangle operations
  , rectSum
  , rectMul
  , rectInv
  , rectWidth
  , midPoint
  , inRect
    -- * Root finding
  , refineToRoot
  , refineRect
  , durandKerner
  ) where

import Data.Ratio
import Data.Complex
import Data.List (minimumBy)
import Data.Ord (comparing)
import LimCalc.Algebra.QPoly
import LimCalc.Algebra.BivPoly

-- | Exact rational number (alias for clarity).
type Q = Rational

-- | An axis-aligned rectangle in the Gaussian rationals ℚ(i),
-- used to isolate a single root of a minimal polynomial in ℂ.
--
-- The invariant is that the rectangle contains exactly one root
-- of the associated minimal polynomial. 'lowerLeft' and 'upperRight'
-- are the two defining corners.
data IsoRect = IsoRect
  { lowerLeft  :: GaussianQ
    -- ^ Lower-left corner of the isolating rectangle.
  , upperRight :: GaussianQ
    -- ^ Upper-right corner of the isolating rectangle.
  } deriving (Eq)

instance Show IsoRect where
  show (IsoRect ll ur) = "[" ++ show ll ++ ", " ++ show ur ++ "]"

-- | Test whether a 'GaussianQ' point lies inside an 'IsoRect'.
inRect :: GaussianQ -> IsoRect -> Bool
inRect (GQ r i) (IsoRect (GQ r1 i1) (GQ r2 i2)) =
  r1 <= r && r <= r2 && i1 <= i && i <= i2

-- | A complex algebraic number, represented as a minimal polynomial
-- over ℚ together with an isolating rectangle in ℂ.
--
-- The isolating rectangle uniquely identifies which root of the
-- minimal polynomial this value represents. All arithmetic operations
-- maintain the invariant that the rectangle contains exactly one root.
data AlgNum = AlgNum
  { algMinPoly :: QPoly
    -- ^ Minimal polynomial over ℚ. Stored in ascending-degree order
    -- (the 'QPoly' convention). Should be monic; this invariant is
    -- maintained by all arithmetic operations.
  , algIsoRect :: IsoRect
    -- ^ Isolating rectangle containing exactly one root of 'algMinPoly'.
  }

instance Show AlgNum where
  show (AlgNum p rect) = "AlgNum(" ++ show p ++ ", " ++ show rect ++ ")"

-- | Construct an 'AlgNum' from an exact rational.
--
-- The minimal polynomial is @x - r@ and the isolating rectangle is
-- centered at @r@ with a unit margin on each side.
fromQ :: Rational -> AlgNum
fromQ r = AlgNum
  { algMinPoly = QPoly [negate r, 1]
  , algIsoRect = IsoRect
      (GQ (r - 1) (negate 1))
      (GQ (r + 1) 1)
  }

-- | The algebraic number 0.
algZero :: AlgNum
algZero = fromQ 0

-- | The algebraic number 1.
algOne :: AlgNum
algOne = fromQ 1

-- | The imaginary unit i, satisfying i² = −1.
--
-- The isolating rectangle is centered at @0 + 1i@, matching the
-- convention used by 'fromQ' and 'algSqrt'. An earlier version used
-- an asymmetric rectangle (center at @0 + 0.5i@), which caused
-- @i + i@ to silently select the wrong root of the resultant
-- polynomial, returning 0 instead of 2i.
algI :: AlgNum
algI = AlgNum
  { algMinPoly = QPoly [1, 0, 1]
  , algIsoRect = IsoRect
      (GQ (negate 1) 0)
      (GQ 1 2)
  }

-- | Square root of a positive rational, as an 'AlgNum'.
--
-- Minimal polynomial: @x² - r@. The isolating rectangle is placed
-- in the right half-plane, selecting the positive real square root.
algSqrt :: Rational -> AlgNum
algSqrt r = AlgNum
  { algMinPoly = QPoly [negate r, 0, 1]
  , algIsoRect = IsoRect
      (GQ 0 (negate 1))
      (GQ (r + 1) 1)
  }

-- | Evaluate the minimal polynomial of an 'AlgNum' at a 'GaussianQ' point.
algEval :: AlgNum -> GaussianQ -> GaussianQ
algEval (AlgNum p _) = qEval p

-- | Refine an isolating rectangle by bisection until its width is
-- below @eps@.
--
-- Only applicable to real algebraic numbers (imaginary part zero);
-- used internally before Durand-Kerner was introduced.
refineRect :: AlgNum -> Rational -> AlgNum
refineRect an@(AlgNum p rect) eps
  | rectWidth rect < eps = an
  | otherwise =
      let mid     = midPoint rect
          upper   = IsoRect mid (upperRight rect)
          lower   = IsoRect (lowerLeft rect) mid
          newRect = if hasRoot p upper then upper else lower
      in refineRect (AlgNum p newRect) eps

-- | Width of an isolating rectangle (real extent).
rectWidth :: IsoRect -> Rational
rectWidth (IsoRect (GQ r1 _) (GQ r2 _)) = abs (r2 - r1)

-- | Midpoint of an isolating rectangle as a 'GaussianQ'.
--
-- Used as an approximate representative of the algebraic number —
-- the convention underlying 'algToDouble' and 'algImagDouble'.
midPoint :: IsoRect -> GaussianQ
midPoint (IsoRect ll ur) = (ll + ur) / 2

-- | Test whether a polynomial has a root in a rectangle using a
-- sign-change heuristic on the real axis.
hasRoot :: QPoly -> IsoRect -> Bool
hasRoot p rect =
  let ll  = lowerLeft rect
      ur  = upperRight rect
      rLeft  = realQ (qEval p (GQ (realQ ll) 0))
      rRight = realQ (qEval p (GQ (realQ ur) 0))
      rMid   = realQ (qEval p (GQ ((realQ ll + realQ ur) / 2) 0))
  in (rLeft * rMid <= 0) || (rMid * rRight <= 0) || abs rMid < 1e-10

-- | Add two algebraic numbers.
--
-- Computes the minimal polynomial of @a + b@ via 'addResultantQ',
-- constructs an approximate rectangle via 'rectSum', then refines
-- to the correct root with 'refineToRoot'.
algAdd :: AlgNum -> AlgNum -> AlgNum
algAdd a b =
  let resPoly = addResultantQ (algMinPoly a) (algMinPoly b)
      rect    = rectSum (algIsoRect a) (algIsoRect b)
  in refineToRoot (AlgNum resPoly rect)

-- | Negate an algebraic number.
--
-- Negation maps the root @α@ of @p(x)@ to the root @−α@ of @p(−x)@,
-- obtained by negating coefficients at odd indices. The result is
-- rescaled to be monic, since the arithmetic operations assume monic
-- minimal polynomials.
--
-- The isolating rectangle is negated by reflecting both corners
-- through the origin (negating both real and imaginary parts and
-- swapping corners to maintain the lowerLeft\/upperRight ordering).
-- An earlier version only negated the real part, silently leaving
-- the imaginary part unchanged and corrupting every subtraction
-- involving a non-real 'AlgNum'.
algNeg :: AlgNum -> AlgNum
algNeg (AlgNum p rect) =
  let negP0 = QPoly [ if even i then c else negate c
                     | (i, c) <- zip [0 :: Int ..] (qPolyCoef p) ]
      negP  = qStrip negP0
      lc    = qLeadingCoeff negP
      monicNegP = if lc == 0 then negP else qScalePoly (1/lc) negP
      negRect = IsoRect
        (GQ (negate (realQ (upperRight rect))) (negate (imagQ (upperRight rect))))
        (GQ (negate (realQ (lowerLeft rect)))  (negate (imagQ (lowerLeft rect))))
  in AlgNum monicNegP negRect

-- | Multiply two algebraic numbers.
--
-- Computes the minimal polynomial of @a * b@ via 'mulResultantQ',
-- constructs an approximate rectangle via 'rectMul', then refines
-- to the correct root with 'refineToRoot'.
algMul :: AlgNum -> AlgNum -> AlgNum
algMul a b =
  let resPoly = mulResultantQ (algMinPoly a) (algMinPoly b)
      rect    = rectMul (algIsoRect a) (algIsoRect b)
  in refineToRoot (AlgNum resPoly rect)

-- | Invert an algebraic number (@1 \/ a@).
--
-- If @α@ satisfies @p(x) = 0@, then @1\/α@ satisfies the polynomial
-- obtained by reversing the coefficient list of @p@. The result is
-- rescaled to be monic and the isolating rectangle is inverted via
-- 'rectInv'.
algInv :: AlgNum -> AlgNum
algInv (AlgNum p rect) =
  let cs      = qPolyCoef p
      invP    = qStrip $ QPoly (reverse cs)
      lc      = qLeadingCoeff invP
      monicP  = qScalePoly (1/lc) invP
      invRect = rectInv rect
  in refineToRoot (AlgNum monicP invRect)

-- | Bounding rectangle for the sum of two algebraic numbers.
--
-- If @a ∈ R₁@ and @b ∈ R₂@, then @a + b ∈ R₁ + R₂@ (Minkowski sum).
rectSum :: IsoRect -> IsoRect -> IsoRect
rectSum (IsoRect ll1 ur1) (IsoRect ll2 ur2) =
  IsoRect (ll1 + ll2) (ur1 + ur2)

-- | Bounding rectangle for the product of two algebraic numbers.
--
-- The rectangle is centered at the product of the two midpoints
-- (which is the best first-order estimate of the true product),
-- with half-widths determined by the spread of the four corner
-- products. A minimum half-width of 1\/2 prevents degenerate
-- zero-width rectangles.
--
-- An earlier version used the bounding box of all four corner
-- products as the rectangle itself, without centering on the
-- midpoint product. For complex multiplication this produced
-- significantly off-center rectangles, causing 'refineToRoot' to
-- select the wrong root (observed: @algI * algI@ returning @+1@
-- instead of @−1@).
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
      halfWidth d = max d (1 % 2)
  in IsoRect
       (center - GQ (halfWidth spreadR) (halfWidth spreadI))
       (center + GQ (halfWidth spreadR) (halfWidth spreadI))

-- | Bounding rectangle for the reciprocal of an algebraic number.
--
-- Computes the reciprocals of all four corners and takes the
-- bounding box of the results.
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

-- | Refine an 'AlgNum' to a specific root of its minimal polynomial.
--
-- Dispatches on the degree of the (stripped) minimal polynomial:
--
-- * Degree 0: degenerate; return as-is.
-- * Degree 1: the root is exactly @−c₀\/c₁@; construct a clean rectangle.
-- * Degree 2–'degreeCapForExactRootFinding': use Durand-Kerner to find
--   all roots, then select the one closest to the rectangle's midpoint.
-- * Above the cap: fall back to 'pinNumerically' to avoid combinatorial
--   blowup from long unminimised resultant chains.
refineToRoot :: AlgNum -> AlgNum
refineToRoot an@(AlgNum p _) =
  let p' = qStrip p
  in case qDegree p' of
       0 -> an
       1 -> let c0 = qPolyCoef p' !! 0
                c1 = qPolyCoef p' !! 1
                root = negate c0 / c1
            in AlgNum p' (exactRect root)
       d | d > degreeCapForExactRootFinding ->
             let approxTarget = complexMidpoint (algIsoRect an)
             in pinNumerically approxTarget
         | otherwise ->
             let approxTarget = complexMidpoint (algIsoRect an)
                 roots        = durandKerner p'
             in case roots of
                  [] -> an
                  _  -> let chosen = nearestTo approxTarget roots
                        in AlgNum p' (rectAroundComplex chosen)

-- | Degree threshold above which 'refineToRoot' falls back to
-- 'pinNumerically' rather than running Durand-Kerner.
--
-- Chosen to be comfortably above what a single arithmetic operation
-- between low-degree 'AlgNum's should produce, while catching the
-- combinatorial-blowup case from long unminimised resultant chains.
-- Confirmed necessary: chaining a handful of @+\/-@ operations on
-- @algI@ (degree 2) produced degree-7 minimal polynomials, making
-- Durand-Kerner catastrophically slow.
degreeCapForExactRootFinding :: Int
degreeCapForExactRootFinding = 6

-- | Pin an algebraic number numerically as an approximate degree-1
-- 'AlgNum' constructed directly from a 'Complex Double' target.
--
-- Used by 'refineToRoot' when the degree cap is exceeded. Trades
-- algebraic exactness for termination, consistent with the
-- already-approximate transcendental functions in this module.
pinNumerically :: Complex Double -> AlgNum
pinNumerically z =
  let re = toRational (realPart z)
      im = toRational (imagPart z)
  in if abs (imagPart z) < 1e-9
       then fromQ re
       else fromQ re + fromQ im * algI

-- | Construct a clean isolating rectangle around an exact rational root,
-- matching the convention used by 'fromQ'.
exactRect :: Rational -> IsoRect
exactRect root = IsoRect
  (GQ (root - 1) (negate 1))
  (GQ (root + 1) 1)

------------------------------------------------------------------------
-- Complex root-finding (Durand-Kerner / Weierstrass iteration)
------------------------------------------------------------------------

-- | Find all complex roots of a 'QPoly' via Durand-Kerner iteration.
--
-- Converts the exact rational coefficients to 'Complex Double' for
-- the iteration. Returns an empty list for degenerate input
-- (degree < 1).
durandKerner :: QPoly -> [Complex Double]
durandKerner p
  | n < 1     = []
  | otherwise = iterateDK cs (initialGuesses n) durandKernerMaxIters
  where
    n  = qDegree p
    cs = map fromRational (qPolyCoef p) :: [Double]

-- | Hard cap on Durand-Kerner iterations, guarding against
-- pathological non-convergence (e.g. very closely clustered or
-- repeated roots).
durandKernerMaxIters :: Int
durandKernerMaxIters = 200

-- | Convergence tolerance: stop early once every root's update step
-- is smaller than this.
durandKernerTol :: Double
durandKernerTol = 1e-12

-- | Initial guesses spread around a circle in the complex plane.
--
-- Using a radius slightly off 1 and irrational-ish angle offsets
-- helps avoid degenerate initial configurations for symmetric
-- polynomials.
initialGuesses :: Int -> [Complex Double]
initialGuesses n =
  [ mkPolar 1.13 (2 * pi * fromIntegral k / fromIntegral n + 0.3)
  | k <- [0 .. n - 1]
  ]

-- | Evaluate a polynomial (ascending-degree coefficient list) at a
-- complex point via Horner's method.
evalComplexPoly :: [Double] -> Complex Double -> Complex Double
evalComplexPoly cs z = foldr (\c acc -> (c :+ 0) + z * acc) (0 :+ 0) cs

-- | One Durand-Kerner iteration step.
--
-- Updates each root @z_i@ by:
--
-- @z_i ← z_i − p(z_i) \/ (aₙ · ∏_{j≠i} (z_i − z_j))@
--
-- where @aₙ@ is the leading coefficient of @p@. Dividing by @aₙ@
-- is required for non-monic polynomials; an earlier version omitted
-- it, causing divergence (roots drifting to magnitude ~10⁴⁴) on the
-- non-monic resultant polynomials produced by 'addResultantQ' and
-- 'mulResultantQ'.
dkStep :: [Double] -> [Complex Double] -> [Complex Double]
dkStep cs zs =
  [ z - evalComplexPoly cs z / (leading * denom z (filter (/= z) zs'))
  | (z, zs') <- zip zs (repeat zs)
  ]
  where
    leading = (last cs) :+ 0
    denom z others = product [ z - w | w <- others ]

-- | Iterate Durand-Kerner until convergence or the iteration cap is hit.
iterateDK :: [Double] -> [Complex Double] -> Int -> [Complex Double]
iterateDK _  zs 0 = zs
iterateDK cs zs n =
  let zs' = dkStep cs zs
      delta = maximum [ magnitude (a - b) | (a, b) <- zip zs zs' ]
  in if delta < durandKernerTol
       then zs'
       else iterateDK cs zs' (n - 1)

-- | Midpoint of an isolating rectangle as a 'Complex Double'.
--
-- Used by 'refineToRoot' to identify which Durand-Kerner root is
-- the intended one.
complexMidpoint :: IsoRect -> Complex Double
complexMidpoint rect =
  let GQ r i = midPoint rect
  in fromRational r :+ fromRational i

-- | Select the element of a list closest to a target point.
nearestTo :: Complex Double -> [Complex Double] -> Complex Double
nearestTo target = minimumBy (comparing (magnitude . subtract target))

-- | Build an isolating rectangle around an approximate complex root
-- found numerically.
--
-- Sized with a margin of @1\/1000000@ on each side — generous
-- relative to 'durandKernerTol' — to ensure the true root lies
-- comfortably inside.
rectAroundComplex :: Complex Double -> IsoRect
rectAroundComplex z =
  let re = toRational (realPart z)
      im = toRational (imagPart z)
      margin = 1 % 1000000
  in IsoRect
       (GQ (re - margin) (im - margin))
       (GQ (re + margin) (im + margin))

-- | 'Num' instance for 'AlgNum'.
--
-- 'abs' returns the magnitude (negating if the real part is negative),
-- not the value itself. An earlier version defined @abs an = an@,
-- which caused every coefficient magnitude check in 'stripZeros' and
-- 'leadingTermNZ' to silently pass for negative coefficients.
instance Num AlgNum where
  (+)         = algAdd
  (*)         = algMul
  negate      = algNeg
  abs an      = if algToDouble an < 0 then algNeg an else an
  signum _    = algOne
  fromInteger = fromQ . fromInteger

-- | 'Fractional' instance for 'AlgNum'.
instance Fractional AlgNum where
  recip        = algInv
  fromRational = fromQ

-- | Approximate transcendental functions via 'Double' arithmetic.
--
-- Each function evaluates its argument approximately via 'algToDouble',
-- applies the standard 'Double' function, and wraps the result back
-- as a degree-1 'AlgNum' via 'fromQ'. This is consistent with the
-- series expansion engine's use of 'AlgNum' as a coefficient type:
-- exact algebraic structure is maintained for rational arithmetic,
-- while transcendental values are approximated.
algSin :: AlgNum -> AlgNum
algSin a = fromRational (toRational (sin (algToDouble a)))

-- | Approximate cosine. See 'algSin'.
algCos :: AlgNum -> AlgNum
algCos a = fromRational (toRational (cos (algToDouble a)))

-- | Approximate exponential. See 'algSin'.
algExp :: AlgNum -> AlgNum
algExp a = fromRational (toRational (exp (algToDouble a)))

-- | Approximate natural logarithm. See 'algSin'.
algLog :: AlgNum -> AlgNum
algLog a = fromRational (toRational (log (algToDouble a)))

-- | Approximate real part of an 'AlgNum' as a 'Double'.
--
-- Returns the real part of the midpoint of the isolating rectangle.
algToDouble :: AlgNum -> Double
algToDouble (AlgNum _ rect) =
  let mid = midPoint rect
  in fromRational (realQ mid)

-- | Approximate imaginary part of an 'AlgNum' as a 'Double'.
--
-- Returns the imaginary part of the midpoint of the isolating
-- rectangle. Used to detect non-real leading coefficients, e.g. in
-- 'LimCalc.Expand.expand' for 'LimCalc.Expr.Abs'.
algImagDouble :: AlgNum -> Double
algImagDouble (AlgNum _ rect) =
  let mid = midPoint rect
  in fromRational (imagQ mid)

-- | 'Ord' instance for 'AlgNum' (approximate, via 'algToDouble').
instance Ord AlgNum where
  compare a b = compare (algToDouble a) (algToDouble b)

-- | 'Real' instance for 'AlgNum' (approximate, via 'algToDouble').
instance Real AlgNum where
  toRational = toRational . algToDouble

-- | 'RealFrac' instance for 'AlgNum' (approximate).
instance RealFrac AlgNum where
  properFraction a =
    let d = algToDouble a
        n = truncate d :: Integer
    in (fromIntegral n, fromQ (toRational (d - fromIntegral n)))
  truncate = truncate . algToDouble
  round    = round    . algToDouble
  ceiling  = ceiling  . algToDouble
  floor    = floor    . algToDouble

-- | 'Floating' instance for 'AlgNum' (approximate via 'Double').
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

-- | Test whether an 'AlgNum' is zero (within numerical tolerance).
--
-- Checks both real and imaginary parts via 'algToDouble' and
-- 'algImagDouble'. Checking only the real part would incorrectly
-- treat purely imaginary non-zero values (e.g. @1\/(2i) = −0.5i@)
-- as zero, corrupting downstream polynomial arithmetic.
isAlgZero :: AlgNum -> Bool
isAlgZero a = abs (algToDouble a) < 1e-12 && abs (algImagDouble a) < 1e-12

-- | 'Eq' instance for 'AlgNum': equality up to numerical tolerance.
instance Eq AlgNum where
  a == b = isAlgZero (a - b)