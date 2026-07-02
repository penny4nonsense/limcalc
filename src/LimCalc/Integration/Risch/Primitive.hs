-- | Risch algorithm: primitive (logarithmic) case.
--
-- Implements integration of rational functions @p(x)\/q(x)@ over a
-- primitive differential field extension — the core of the Risch
-- algorithm for the case where no exponential extension is present.
--
-- = Algorithm outline
--
-- 1. /Hermite reduction/ ('LimCalc.Algebra.RationalFunction.hermiteReduce'):
--    split @f = g' + h@ where @g@ is a rational function (rational
--    integral part) and @h@ is a proper fraction over a squarefree
--    denominator.
--
-- 2. /Rothstein-Trager/ ('rothsteinTrager'): determine whether @h@
--    has an elementary integral by computing the resultant polynomial
--    @R(z) = res_x(d, a − z·d')@ and checking whether its roots are
--    rational. If so, the integral is @∑ cᵢ · log(gcd(d, a − cᵢ·d'))@.
--    If not, the integral is non-elementary.
--
-- = AlgNum coefficients
--
-- The implementation is generalised from @RatFun Double@ to
-- @RatFun AlgNum@ to support the trig integration path in
-- 'LimCalc.Integration.Risch.Exponential', which introduces genuinely
-- complex coefficients via the Euler substitution @t = e^(ix)@. All
-- existing arithmetic (@hermiteReduce@, @gcdPoly@, etc.) carries
-- through unchanged via 'AlgNum'\'s 'Num' and 'Fractional' instances.
--
-- = Rational vs complex roots
--
-- 'rothsteinTrager' uses 'findRationalRoots' (not 'findComplexRoots')
-- because the classical elementary\/non-elementary criterion for
-- integration over @ℚ(x)@ requires the resultant's roots to be
-- /rational/. Allowing any complex root would always declare the
-- integral elementary (with complex-coefficient logarithms), silently
-- misclassifying genuinely non-elementary-over-the-reals integrals
-- like @1\/(x²+1)@. 'LimCalc.Integration.Risch.Exponential' uses
-- 'findComplexRoots' instead, since that module deliberately extends
-- the field to include @i@.
module LimCalc.Integration.Risch.Primitive
  ( -- * Integration
    integratePrimitive
  , PrimitiveResult (..)
    -- * Rothstein-Trager
  , rothsteinTrager
  , logDerivativeCheck
  , resultantPoly
    -- * Root finding
  , findRationalRoots
  , findComplexRoots
  , snapToRational
    -- * Interpolation
  , interpolate
  , lagrangeBasis
    -- * Partial fractions
  , partialFractions
    -- * Conversion utilities
  , polyToExpr
  , ratFunToExpr
  , algNumToExpr
  , algNumToRational
  , complexToAlgNum
  , addExpr
  ) where

import Data.Complex (Complex, realPart, imagPart)
import Data.List (minimumBy)
import Data.Ratio ((%))
import LimCalc.Algebra.Poly
import LimCalc.Algebra.RationalFunction
import LimCalc.Differentiation.DiffField
import LimCalc.Core.Expr
import LimCalc.Algebra.AlgNum
import LimCalc.Algebra.QPoly (QPoly(..))

-- | Result of the Risch primitive-case integration algorithm.
data PrimitiveResult
  = PrimitiveElementary Expr
    -- ^ The integral is elementary; the 'Expr' is the antiderivative.
  | PrimitiveNonElementary
    -- ^ The integral is non-elementary (Rothstein-Trager found no
    -- rational roots for the resultant polynomial).
  | PrimitiveError String
    -- ^ An internal error occurred.
  deriving (Show, Eq)

-- | Integrate a rational function in the primitive (logarithmic) case.
--
-- Applies Hermite reduction followed by Rothstein-Trager. Returns the
-- antiderivative as an 'Expr' if elementary, or 'PrimitiveNonElementary'
-- if the Rothstein-Trager criterion fails.
integratePrimitive :: RatFun AlgNum -> DiffField -> PrimitiveResult
integratePrimitive rf field =
  let (g, reduced) = hermiteReduce rf
      rtResult     = rothsteinTrager reduced field
  in case rtResult of
       Left "NonElementary" -> PrimitiveNonElementary
       Left err             -> PrimitiveError err
       Right terms          ->
         let logPart = foldr addExpr (Const 0)
               [ Mul (algNumToExpr c) (Log (polyToExpr u))
               | (c, u) <- terms
               ]
             gExpr = ratFunToExpr g
         in PrimitiveElementary (addExpr gExpr logPart)

-- | Rothstein-Trager algorithm for the primitive case.
--
-- Given a proper fraction @a\/d@ over a squarefree denominator @d@,
-- computes the resultant polynomial @R(z) = res_x(d, a − z·d')@
-- and finds its rational roots @c₁, …, cₙ@. The integral is then
-- @∑ cᵢ · log(gcd(d, a − cᵢ·d'))@.
--
-- Detects the log-derivative case @a = c·d'@ directly via
-- 'logDerivativeCheck', bypassing the resultant computation (which
-- degenerates in this case).
--
-- Returns @Left \"NonElementary\"@ if no rational roots are found.
rothsteinTrager :: RatFun AlgNum -> DiffField -> Either String [(AlgNum, Poly AlgNum)]
rothsteinTrager (RatFun a d) field =
  let d' = diffPoly d
  in case logDerivativeCheck a d d' of
       Just c  -> Right [(c, d)]
       Nothing ->
         let rPoly = resultantPoly a d d'
         in case findRationalRoots rPoly of
              Nothing    -> Left "NonElementary"
              Just roots ->
                Right [ (c, gcdPoly d (subPoly a (scalePoly c d')))
                      | c <- roots
                      , not (isAlgZero c)
                      ]

-- | Check whether @a = c · d'@ for some constant @c@.
--
-- This is the log-derivative edge case: @∫ c · (d'\/d) dx = c · log(d)@.
-- The Rothstein-Trager resultant degenerates in this case (the
-- resultant has a root at @z = c@ but the GCD computation gives back
-- all of @d@), so it is detected and handled directly.
logDerivativeCheck :: Poly AlgNum -> Poly AlgNum -> Poly AlgNum -> Maybe AlgNum
logDerivativeCheck a _ d'
  | degree d' < 0 = Nothing
  | otherwise =
      let q = quotPoly a d'
          r = remPoly a d'
      in if degree q == 0 && isZeroPoly r
           then Just (leadingCoeff q)
           else Nothing
  where
    isZeroPoly p = degree p < 0

-- | Compute the Rothstein-Trager resultant polynomial
-- @R(z) = res_x(d, a − z·d')@.
--
-- Evaluates the resultant at @2·deg(d) + 1@ integer points and
-- reconstructs the polynomial by Lagrange interpolation.
resultantPoly :: Poly AlgNum -> Poly AlgNum -> Poly AlgNum -> Poly AlgNum
resultantPoly a d d' =
  let n      = degree d
      points = [ fromInteger (fromIntegral i) | i <- [-n..n] ] :: [AlgNum]
      values = [ resultant d (subPoly a (scalePoly z d')) | z <- points ]
  in interpolate (polyVar d) (zip points values)

-- | Lagrange interpolation through a set of @(x, y)@ pairs.
interpolate :: String -> [(AlgNum, AlgNum)] -> Poly AlgNum
interpolate x points =
  foldr (addPoly . lagrangeBasis x points)
        (zeroPoly x)
        (zip [0..] points)

-- | One Lagrange basis polynomial for the @i@th interpolation point.
lagrangeBasis :: String -> [(AlgNum, AlgNum)] -> (Int, (AlgNum, AlgNum)) -> Poly AlgNum
lagrangeBasis x points (i, (xi, yi)) =
  let others = [ (j, xj) | (j, (xj, _)) <- zip [0..] points, j /= i ]
      basis  = foldr mulBasis (constPoly x 1) others
      scale  = yi / evalPoly basis xi
  in scalePoly scale basis
  where
    mulBasis (_, xj) p =
      mulPoly p (subPoly (Poly x [0, 1]) (constPoly x xj))

-- | Find rational (real, negligible imaginary part) roots of a
-- polynomial with 'AlgNum' coefficients.
--
-- Used by 'rothsteinTrager' for the classical elementary\/non-elementary
-- criterion, which requires rational roots specifically. Built on top
-- of 'findComplexRoots': find all complex roots numerically, then
-- retain only those with negligible imaginary part, snapping each to
-- the nearest rational @p\/q@ with @|q| ≤ 100@.
--
-- Returns 'Nothing' if no rational roots are found and the degree is
-- positive (indicating the integral is non-elementary).
findRationalRoots :: Poly AlgNum -> Maybe [AlgNum]
findRationalRoots p
  | degree p < 0 = Just []
  | otherwise    =
      case findComplexRoots p of
        Nothing -> Nothing
        Just allRoots ->
          let snapped   = [ snapToRational r | r <- allRoots, isReal r ]
              realRoots = [ r | r <- snapped, verifiedRoot r ]
          in if null realRoots && degree p > 0
               then Nothing
               else Just realRoots
  where
    isReal r       = abs (algImagDouble r) < 1e-6
    verifiedRoot r = abs (algToDouble (evalPoly p r)) < 1e-4

-- | Snap an 'AlgNum' to the nearest rational @p\/q@ with @|q| ≤ 100@.
--
-- The Rothstein-Trager resultant has rational coefficients (when the
-- input has rational coefficients), so its roots are exactly rational.
-- Durand-Kerner finds them approximately; this corrects the precision
-- so that downstream GCD computations work correctly.
snapToRational :: AlgNum -> AlgNum
snapToRational r =
  let x    = algToDouble r
      best = minimumBy (\a b -> compare (abs (algToDouble a - x))
                                        (abs (algToDouble b - x)))
               [ fromQ (p' % q')
               | q' <- [1..100]
               , p' <- [ round (x * fromIntegral q') ]
               ]
  in if abs (algToDouble best - x) < 1e-4 then best else r

-- | Find all complex roots of a polynomial with 'AlgNum' coefficients.
--
-- Approximates each 'AlgNum' coefficient as a 'Rational' via its
-- 'Double' midpoint, then applies 'LimCalc.AlgNum.durandKerner'.
-- Returns 'Nothing' only for degenerate input (degree < 1).
--
-- Used directly by 'LimCalc.Risch.Exponential', where the field has
-- been extended to include @i@ and complex roots are the correct
-- criterion. 'rothsteinTrager' uses 'findRationalRoots' instead.
findComplexRoots :: Poly AlgNum -> Maybe [AlgNum]
findComplexRoots p
  | degree p < 0 = Just []
  | otherwise    =
      let approxQPoly = QPoly (map algNumToRational (polyCoef p))
          roots       = durandKerner approxQPoly
      in if null roots && degree p > 0
         then Nothing
         else Just (map complexToAlgNum roots)

-- | Approximate an 'AlgNum' as a 'Rational' via its 'Double' midpoint.
-- Used to construct approximate input for Durand-Kerner.
algNumToRational :: AlgNum -> Rational
algNumToRational = toRational . algToDouble

-- | Lift a 'Complex Double' root (from Durand-Kerner) to an 'AlgNum'.
complexToAlgNum :: Complex Double -> AlgNum
complexToAlgNum z =
  let re = realPart z
      im = imagPart z
  in if abs im < 1e-9
       then fromQ (toRational re)
       else fromQ (toRational re) + fromQ (toRational im) * algI

-- | Convert a 'Poly AlgNum' to an 'Expr'.
polyToExpr :: Poly AlgNum -> Expr
polyToExpr (Poly _ [])  = Const 0
polyToExpr (Poly x cs)  =
  foldr1 Add [ termToExpr x n c
             | (n, c) <- zip [0..] cs
             , not (isAlgZero c) ]
  where
    termToExpr _ 0 c = algNumToExpr c
    termToExpr x 1 c = Mul (algNumToExpr c) (Var x)
    termToExpr x n c = Mul (algNumToExpr c)
                           (Pow (Var x) (Const (fromIntegral (n :: Int))))

-- | Convert a 'RatFun AlgNum' to an 'Expr'.
-- If the denominator is 1, returns just the numerator expression.
ratFunToExpr :: RatFun AlgNum -> Expr
ratFunToExpr (RatFun p q)
  | degree q == 0 && isAlgZero (leadingCoeff q - algOne) = polyToExpr p
  | otherwise = Div (polyToExpr p) (polyToExpr q)

-- | Convert an 'AlgNum' to an 'Expr'.
--
-- Produces @Const re@ for real values, @Mul (Const im) I@ for purely
-- imaginary values, and @Add (Const re) (Mul (Const im) I)@ for
-- complex values. Avoids a spurious @+ 0·I@ in the common real case.
--
-- Note: 'Expr'\'s 'Const' is still @Double@-backed; this is a
-- known limitation documented in 'LimCalc.Expr'.
algNumToExpr :: AlgNum -> Expr
algNumToExpr a =
  let re = algToDouble a
      im = algImagDouble a
  in if abs im < 1e-12
       then Const re
       else if abs re < 1e-12
              then Mul (Const im) I
              else Add (Const re) (Mul (Const im) I)

-- | Add two expressions, simplifying @0 + e = e@.
addExpr :: Expr -> Expr -> Expr
addExpr (Const 0) e = e
addExpr e (Const 0) = e
addExpr e1 e2       = Add e1 e2

-- | Partial fraction decomposition of a proper rational function.
--
-- Uses 'rothsteinTrager' to extract the irreducible factor structure
-- of the denominator. For each linear factor @(x − r)@, computes the
-- residue coefficient via evaluation at @r@. For higher-degree
-- irreducible factors, returns the logarithmic-derivative form
-- @c · d'\/d@.
--
-- Returns @[rf]@ unchanged if the denominator is already irreducible
-- or has no rational-coefficient decomposition.
--
-- Precondition: @deg(p) < deg(q)@. Call
-- 'LimCalc.RationalFunction.ratProperFraction' first for improper inputs.
partialFractions :: RatFun AlgNum -> [RatFun AlgNum]
partialFractions rf@(RatFun p q) =
  let var   = polyVar p
      field = baseField var
  in case rothsteinTrager rf field of
       Left _      -> [rf]
       Right []    -> [rf]
       Right [_]   -> [rf]
       Right pairs ->
         let factors = concatMap (splitFactor var) (dedupFactors (map snd pairs))
         in concatMap (residueTerm var p q) factors
  where
    splitFactor var di
      | degree di <= 1 = [di]
      | otherwise =
          let one      = constPoly var algOne
              subterms = partialFractions (RatFun one di)
          in map (\(RatFun _ d) -> d) subterms

    residueTerm var p q di
      | degree di == 1 =
          let coeffs = polyCoef di
              r      = algMul (algNeg (coeffs !! 0)) (algInv (coeffs !! 1))
              pr     = evalPoly p r
              cofact = quotPoly q di
              cofr   = evalPoly cofact r
              coeff  = algMul pr (algInv cofr)
          in [RatFun (constPoly var coeff) di]
      | otherwise =
          let field = baseField var
          in case rothsteinTrager (RatFun p di) field of
               Right ((c,_):_) ->
                 [RatFun (mulPoly (constPoly var c) (diffPoly di)) di]
               _ -> [RatFun p di]

    dedupFactors :: [Poly AlgNum] -> [Poly AlgNum]
    dedupFactors [] = []
    dedupFactors (x:xs) =
      x : dedupFactors (filter (\y -> not (polyEq x y)) xs)

    polyEq :: Poly AlgNum -> Poly AlgNum -> Bool
    polyEq a b = length (polyCoef a) == length (polyCoef b) &&
      all (\(ca, cb) -> isAlgZero (algAdd ca (algNeg cb)))
          (zip (polyCoef a) (polyCoef b))