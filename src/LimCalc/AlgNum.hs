module LimCalc.AlgNum where

import Data.Ratio
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
algI :: AlgNum
algI = AlgNum
  { algMinPoly = QPoly [1, 0, 1]
  , algIsoRect = IsoRect
      (GQ (negate 1) 0)
      (GQ 1 1)
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
      negRect = IsoRect
        (GQ (negate (realQ (upperRight rect))) (imagQ (lowerLeft rect)))
        (GQ (negate (realQ (lowerLeft rect))) (imagQ (upperRight rect)))
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

-- | Product of two isolating rectangles
rectMul :: IsoRect -> IsoRect -> IsoRect
rectMul r1 r2 =
  let corners = [ lowerLeft r1 * lowerLeft r2
                , lowerLeft r1 * upperRight r2
                , upperRight r1 * lowerLeft r2
                , upperRight r1 * upperRight r2
                ]
      minR = minimum (map realQ corners)
      maxR = maximum (map realQ corners)
      minI = minimum (map imagQ corners)
      maxI = maximum (map imagQ corners)
  in IsoRect (GQ minR minI) (GQ maxR maxI)

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
-- For genuinely higher-degree results (degree >= 2), we still fall
-- back to refineRect/hasRoot, which remains a real-axis-only
-- approximation pending a true complex root isolation procedure
-- (see hasRoot's TODO).
refineToRoot :: AlgNum -> AlgNum
refineToRoot an@(AlgNum p _) =
  let p' = qStrip p
  in case qDegree p' of
       0 -> an  -- degenerate: no roots; leave rectangle as-is
       1 -> let c0 = qPolyCoef p' !! 0
                c1 = qPolyCoef p' !! 1
                root = negate c0 / c1
            in AlgNum p' (exactRect root)
       _ -> refineRect (AlgNum p' (algIsoRect an)) (1 % 1000)

-- | Construct a clean isolating rectangle around an exact rational
-- root, matching the convention used by 'fromQ'.
exactRect :: Rational -> IsoRect
exactRect root = IsoRect
  (GQ (root - 1) (negate 1))
  (GQ (root + 1) 1)

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