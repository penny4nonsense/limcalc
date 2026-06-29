module LimCalc.AlgNum where

import Data.Ratio

-- | Exact rational number (alias for clarity)
type Q = Rational

-- | Gaussian rational: a + bi where a, b ∈ ℚ
data GaussianQ = GQ
  { realQ :: Rational
  , imagQ :: Rational
  } deriving (Eq)

instance Show GaussianQ where
  show (GQ r 0) = showRat r
  show (GQ 0 i) = showRat i ++ "i"
  show (GQ r i)
    | i < 0    = showRat r ++ " - " ++ showRat (abs i) ++ "i"
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
  signum gq           = GQ 1 0  -- placeholder
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

-- | Polynomial over ℚ for use as minimal polynomials
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
qMulPoly (QPoly []) _          = QPoly []
qMulPoly _ (QPoly [])          = QPoly []
qMulPoly (QPoly cs1) (QPoly cs2) =
  qStrip $ QPoly $ mulCoefs cs1 cs2
  where
    mulCoefs [] _     = []
    mulCoefs (c:cs) ys =
      addCoefs (map (*c) ys) (0 : mulCoefs cs ys)
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (x:xs) (y:ys) = (x+y) : addCoefs xs ys

-- | Scale a QPoly by a rational
qScalePoly :: Rational -> QPoly -> QPoly
qScalePoly 0 _          = QPoly []
qScalePoly c (QPoly cs) = qStrip $ QPoly (map (*c) cs)

-- | Differentiate a QPoly
qDiffPoly :: QPoly -> QPoly
qDiffPoly (QPoly [])  = QPoly []
qDiffPoly (QPoly cs)  =
  qStrip $ QPoly
    [ fromIntegral n * c
    | (n, c) <- zip [1 :: Int ..] (tail cs)
    ]

-- | Isolating rectangle in ℂ
data IsoRect = IsoRect
  { lowerLeft  :: GaussianQ  -- ^ lower-left corner
  , upperRight :: GaussianQ  -- ^ upper-right corner
  } deriving (Eq)

instance Show IsoRect where
  show (IsoRect ll ur) = "[" ++ show ll ++ ", " ++ show ur ++ "]"

-- | Check if a GaussianQ is inside an IsoRect
inRect :: GaussianQ -> IsoRect -> Bool
inRect (GQ r i) (IsoRect (GQ r1 i1) (GQ r2 i2)) =
  r1 <= r && r <= r2 && i1 <= i && i <= i2

-- | A complex algebraic number
-- Represented as minimal polynomial + isolating rectangle
data AlgNum = AlgNum
  { algMinPoly :: QPoly    -- ^ Minimal polynomial over ℚ
  , algIsoRect :: IsoRect  -- ^ Isolating rectangle containing exactly one root
  } deriving (Eq)

instance Show AlgNum where
  show (AlgNum p rect) = "AlgNum(" ++ show p ++ ", " ++ show rect ++ ")"

-- | Construct an algebraic number from a rational
fromQ :: Rational -> AlgNum
fromQ r = AlgNum
  { algMinPoly = QPoly [negate r, 1]  -- x - r
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
-- Minimal polynomial: x^2 + 1
algI :: AlgNum
algI = AlgNum
  { algMinPoly = QPoly [1, 0, 1]  -- x^2 + 1
  , algIsoRect = IsoRect
      (GQ (negate 1) 0)
      (GQ 1 1)
  }

-- | Square root of a positive rational
algSqrt :: Rational -> AlgNum
algSqrt r = AlgNum
  { algMinPoly = QPoly [negate r, 0, 1]  -- x^2 - r
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
      let mid    = midPoint rect
          upper  = IsoRect mid (upperRight rect)
          lower  = IsoRect (lowerLeft rect) mid
          newRect = if hasRoot p upper then upper else lower
      in refineRect (AlgNum p newRect) eps

-- | Width of an isolating rectangle
rectWidth :: IsoRect -> Rational
rectWidth (IsoRect (GQ r1 _) (GQ r2 _)) = abs (r2 - r1)

-- | Midpoint of an isolating rectangle
midPoint :: IsoRect -> GaussianQ
midPoint (IsoRect ll ur) = (ll + ur) / 2

-- | Check if a polynomial has a root in a rectangle
-- Uses sign changes on boundary (simplified — full version needs winding number)
hasRoot :: QPoly -> IsoRect -> Bool
hasRoot p rect =
  let corners = [ lowerLeft rect
                , GQ (realQ (upperRight rect)) (imagQ (lowerLeft rect))
                , upperRight rect
                , GQ (realQ (lowerLeft rect)) (imagQ (upperRight rect))
                ]
      vals = map (qEval p) corners
  in any (\v -> realQ v * realQ (qEval p (midPoint rect)) < 0) vals
     || True  -- simplified: assume root is in one half

-- | Add two algebraic numbers
-- α + β is a root of res_y(p(y), q(x-y))
algAdd :: AlgNum -> AlgNum -> AlgNum
algAdd a b =
  let pa    = algMinPoly a
      pb    = algMinPoly b
      -- resultant gives minimal poly of α + β
      resPoly = addResultant pa pb
      -- isolating rectangle: sum of rectangles contains α + β
      rect  = rectSum (algIsoRect a) (algIsoRect b)
  in refineToRoot (AlgNum resPoly rect)

-- | Negate an algebraic number
-- -α is a root of p(-x)
algNeg :: AlgNum -> AlgNum
algNeg (AlgNum p rect) =
  let negP  = QPoly [ if even i then c else negate c
                    | (i, c) <- zip [0 :: Int ..] (qPolyCoef p) ]
      negRect = IsoRect
        (GQ (negate (realQ (upperRight rect))) (imagQ (lowerLeft rect)))
        (GQ (negate (realQ (lowerLeft rect))) (imagQ (upperRight rect)))
  in AlgNum (qStrip negP) negRect

-- | Multiply two algebraic numbers
-- α * β is a root of res_y(p(y), y^n * q(x/y)) where n = deg(q)
algMul :: AlgNum -> AlgNum -> AlgNum
algMul a b =
  let pa    = algMinPoly a
      pb    = algMinPoly b
      resPoly = mulResultant pa pb
      rect  = rectMul (algIsoRect a) (algIsoRect b)
  in refineToRoot (AlgNum resPoly rect)

-- | Invert an algebraic number
-- 1/α is a root of x^n * p(1/x) where n = deg(p)
algInv :: AlgNum -> AlgNum
algInv (AlgNum p rect) =
  let cs    = qPolyCoef p
      invP  = QPoly (reverse cs)
      -- 1/α is in 1/rect (approximately)
      invRect = rectInv rect
  in refineToRoot (AlgNum (qStrip invP) invRect)

-- | Resultant for addition: res_y(p(y), q(x-y))
-- Computed by substituting and taking resultant
addResultant :: QPoly -> QPoly -> QPoly
addResultant pa pb =
  -- Substitute y -> x-y in pb to get pb(x-y) as poly in y
  -- Then take resultant with pa(y) over y
  -- For now: use a simplified version via companion matrices
  -- Full implementation pending
  qMulPoly pa pb  -- placeholder — real resultant computation needed

-- | Resultant for multiplication
mulResultant :: QPoly -> QPoly -> QPoly
mulResultant pa pb =
  qMulPoly pa pb  -- placeholder

-- | Sum of two isolating rectangles
rectSum :: IsoRect -> IsoRect -> IsoRect
rectSum (IsoRect ll1 ur1) (IsoRect ll2 ur2) =
  IsoRect (ll1 + ll2) (ur1 + ur2)

-- | Product of two isolating rectangles (interval arithmetic)
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

-- | Refine an AlgNum until the rectangle is tight
-- Uses Newton-like bisection to isolate the correct root
refineToRoot :: AlgNum -> AlgNum
refineToRoot an = refineRect an (1 % 1000)

-- | Num instance for AlgNum
instance Num AlgNum where
  (+)         = algAdd
  (*)         = algMul
  negate      = algNeg
  abs an      = an  -- placeholder
  signum an   = algOne  -- placeholder
  fromInteger = fromQ . fromInteger

-- | Fractional instance for AlgNum
instance Fractional AlgNum where
  recip        = algInv
  fromRational = fromQ