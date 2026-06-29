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
  } deriving (Eq)

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
  let corners = [ lowerLeft rect
                , GQ (realQ (upperRight rect)) (imagQ (lowerLeft rect))
                , upperRight rect
                , GQ (realQ (lowerLeft rect)) (imagQ (upperRight rect))
                ]
  in any (\v -> realQ v * realQ (qEval p (midPoint rect)) < 0)
         (map (qEval p) corners)
     || True

-- | Add two algebraic numbers
algAdd :: AlgNum -> AlgNum -> AlgNum
algAdd a b =
  let resPoly = addResultantQ (algMinPoly a) (algMinPoly b)
      rect    = rectSum (algIsoRect a) (algIsoRect b)
  in refineToRoot (AlgNum resPoly rect)

-- | Negate an algebraic number
algNeg :: AlgNum -> AlgNum
algNeg (AlgNum p rect) =
  let negP = QPoly [ if even i then c else negate c
                   | (i, c) <- zip [0 :: Int ..] (qPolyCoef p) ]
      negRect = IsoRect
        (GQ (negate (realQ (upperRight rect))) (imagQ (lowerLeft rect)))
        (GQ (negate (realQ (lowerLeft rect))) (imagQ (upperRight rect)))
  in AlgNum (qStrip negP) negRect

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

-- | Refine until rectangle is tight
refineToRoot :: AlgNum -> AlgNum
refineToRoot an = refineRect an (1 % 1000)

-- | Num instance for AlgNum
instance Num AlgNum where
  (+)         = algAdd
  (*)         = algMul
  negate      = algNeg
  abs an      = an
  signum _    = algOne
  fromInteger = fromQ . fromInteger

-- | Fractional instance for AlgNum
instance Fractional AlgNum where
  recip        = algInv
  fromRational = fromQ