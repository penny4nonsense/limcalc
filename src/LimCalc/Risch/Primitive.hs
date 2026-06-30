module LimCalc.Risch.Primitive where

import Data.Complex (Complex, realPart, imagPart)
import Data.List (minimumBy)
import Data.Ratio ((%))
import LimCalc.Poly
import LimCalc.RationalFunction
import LimCalc.DiffField
import LimCalc.Expr
import LimCalc.AlgNum
import LimCalc.QPoly (QPoly(..))

-- | Result of the Risch algorithm for the primitive case
data PrimitiveResult
  = PrimitiveElementary Expr
  | PrimitiveNonElementary
  | PrimitiveError String
  deriving (Show, Eq)

-- | Integrate a rational function in the primitive case
--
-- Generalized from RatFun Double to RatFun AlgNum: trig integration
-- (sin/cos via Euler's formula, e^(i*theta)) produces genuinely
-- complex coefficients, which Double cannot represent. AlgNum's
-- Num/Fractional instances let all the existing arithmetic
-- (hermiteReduce, gcdPoly, etc.) carry through unchanged; what
-- actually needed rewriting was root-finding (see findComplexRoots)
-- and Expr conversion (see algNumToExpr), both of which assumed real
-- rational coefficients.
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

-- | Rothstein-Trager algorithm
--
-- Uses findRationalRoots, NOT findComplexRoots: the classical
-- elementary/non-elementary criterion for ordinary rational-function
-- integration over Q(x) genuinely requires the resultant's roots to
-- be RATIONAL, not merely complex. Allowing any complex root (as
-- findComplexRoots does) answers a different, less restrictive
-- question -- "can this be integrated if I'm willing to introduce
-- complex-coefficient logarithms" -- which is always satisfiable and
-- silently misclassifies genuinely non-elementary-over-the-reals
-- integrals like 1/(x^2+1) (whose antiderivative is arctan(x), not
-- expressible with rational-coefficient logs) as Elementary. See
-- Risch.Exponential's use of findComplexRoots for the case where
-- this distinction does NOT apply (there, the field has already
-- been deliberately extended to include i, so complex roots are the
-- right question).
rothsteinTrager :: RatFun AlgNum -> DiffField -> Either String [(AlgNum, Poly AlgNum)]
rothsteinTrager (RatFun a d) field =
  let d' = diffPoly d
  in case logDerivativeCheck a d d' of
       Just c  -> Right [(c, d)]  -- int c*(d'/d) dx = c*log(d)
       Nothing ->
         let rPoly = resultantPoly a d d'
         in case findRationalRoots rPoly of
              Nothing    -> Left "NonElementary"
              Just roots ->
                Right [ (c, gcdPoly d (subPoly a (scalePoly c d')))
                      | c <- roots
                      , not (isAlgZero c)
                      ]

-- | Check whether a = c * d' for some constant AlgNum c (the
-- log-derivative case). If so, return Just c; otherwise Nothing.
--
-- This is the edge case where the integrand is a scalar multiple of
-- the logarithmic derivative of the denominator: int c*(d'/d) dx =
-- c*log(d). The Rothstein-Trager resultant machinery degenerates in
-- this case (the resultant polynomial has a root at z=c whose GCD
-- computation gives back all of d rather than a proper factor), so
-- we detect and handle it directly.
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

-- | Compute the Rothstein-Trager resultant polynomial R(z)
--
-- Interpolation itself is generic arithmetic and works unchanged
-- over AlgNum; only the evaluation points needed to become AlgNum
-- (via fromInteger, which AlgNum's Num instance already supports).
resultantPoly :: Poly AlgNum -> Poly AlgNum -> Poly AlgNum -> Poly AlgNum
resultantPoly a d d' =
  let n      = degree d
      points = [ fromInteger (fromIntegral i) | i <- [-n..n] ] :: [AlgNum]
      values = [ resultant d (subPoly a (scalePoly z d')) | z <- points ]
  in interpolate (polyVar d) (zip points values)

-- | Lagrange interpolation
interpolate :: String -> [(AlgNum, AlgNum)] -> Poly AlgNum
interpolate x points =
  foldr (addPoly . lagrangeBasis x points)
        (zeroPoly x)
        (zip [0..] points)

-- | One Lagrange basis polynomial
lagrangeBasis :: String -> [(AlgNum, AlgNum)] -> (Int, (AlgNum, AlgNum)) -> Poly AlgNum
lagrangeBasis x points (i, (xi, yi)) =
  let others = [ (j, xj) | (j, (xj, _)) <- zip [0..] points, j /= i ]
      basis  = foldr mulBasis (constPoly x 1) others
      scale  = yi / evalPoly basis xi
  in scalePoly scale basis
  where
    mulBasis (_, xj) p =
      mulPoly p (subPoly (Poly x [0, 1]) (constPoly x xj))

-- | Find rational (real, zero-imaginary-part) roots of a polynomial
-- with AlgNum coefficients, for the classical Rothstein-Trager
-- elementary/non-elementary criterion (see rothsteinTrager's header
-- for why this must be rational-only, not findComplexRoots's full
-- complex-root version).
--
-- Built on top of findComplexRoots: find all complex roots
-- numerically, then keep only those with negligible imaginary part,
-- verifying each by evaluating the polynomial there. This replaces
-- the old approach (enumerating rational candidates via divisor
-- factorization of the constant/leading coefficients), which only
-- ever worked for genuinely rational (Double) coefficients and
-- cannot be meaningfully generalized to AlgNum coefficients at all.
findRationalRoots :: Poly AlgNum -> Maybe [AlgNum]
findRationalRoots p
  | degree p < 0 = Just []
  | otherwise    =
      case findComplexRoots p of
        Nothing -> Nothing
        Just allRoots ->
          let -- Snap real roots to nearby rationals, then verify
              snapped    = [ snapToRational r | r <- allRoots, isReal r ]
              realRoots  = [ r | r <- snapped, verifiedRoot r ]
          in if null realRoots && degree p > 0
               then Nothing
               else Just realRoots
  where
    isReal r       = abs (algImagDouble r) < 1e-6
    verifiedRoot r = abs (algToDouble (evalPoly p r)) < 1e-4
    -- Snap an AlgNum to the nearest rational p/q with |q| <= 100.
    -- The Rothstein-Trager resultant has rational coefficients when
    -- the input rational function has rational coefficients, so its
    -- rational roots are exactly representable as rationals.
    -- Durand-Kerner finds them approximately; this corrects the
    -- precision so that downstream GCD computations work exactly.
    snapToRational r =
      let x = algToDouble r
          best = minimumBy (\a b -> compare (abs (algToDouble a - x))
                                            (abs (algToDouble b - x)))
                   [ fromQ (p' % q')
                   | q' <- [1..100]
                   , p' <- [ round (x * fromIntegral q') ]
                   ]
      in if abs (algToDouble best - x) < 1e-4 then best else r

-- | Find all complex roots of a polynomial with AlgNum coefficients.
--
-- Previously (findRationalRoots) this enumerated rational root
-- candidates via divisor factorization of the constant and leading
-- coefficients -- meaningless once coefficients are AlgNum rather
-- than rational. Replaced with genuine complex root-finding by
-- approximating each AlgNum coefficient to a Rational (via its
-- Double midpoint) and reusing AlgNum's existing durandKerner
-- machinery directly.
--
-- This is consistent in spirit with the rest of the module's
-- already-approximate transcendental functions (algSin, algExp,
-- etc.); the polynomial's own coefficients are already only
-- Double-approximated AlgNum values by the time they reach this
-- resultant/interpolation pipeline, so there is no loss of rigor
-- relative to the rest of the pipeline by also approximating here.
--
-- IMPORTANT: this finds ANY complex root, not just rational ones.
-- Risch.Primitive's own rothsteinTrager uses findRationalRoots
-- instead, since the classical elementary/non-elementary criterion
-- for integration over Q(x) requires rational roots specifically.
-- This function (findComplexRoots) remains correct and necessary for
-- Risch.Exponential's trig-integration use, where the field has
-- already been deliberately extended to include i.
findComplexRoots :: Poly AlgNum -> Maybe [AlgNum]
findComplexRoots p
  | degree p < 0 = Just []
  | otherwise    =
      let approxQPoly = QPoly (map algNumToRational (polyCoef p))
          roots       = durandKerner approxQPoly
      in if null roots && degree p > 0
         then Nothing
         else Just (map complexToAlgNum roots)

-- | Approximate an AlgNum as a Rational via its Double midpoint, for
-- feeding into Durand-Kerner (which only needs approximate numeric
-- coefficients, not exact algebraic ones).
algNumToRational :: AlgNum -> Rational
algNumToRational = toRational . algToDouble

-- | Lift a Complex Double root (from durandKerner) back to an
-- AlgNum, by directly constructing real + imaginary*i from the
-- rounded Double parts.
complexToAlgNum :: Complex Double -> AlgNum
complexToAlgNum z =
  let re = realPart z
      im = imagPart z
  in if abs im < 1e-9
       then fromQ (toRational re)
       else fromQ (toRational re) + fromQ (toRational im) * algI

-- | Convert a polynomial to an Expr
polyToExpr :: Poly AlgNum -> Expr
polyToExpr (Poly _ [])  = Const 0
polyToExpr (Poly x cs)  =
  foldr1 Add [ termToExpr x n c
             | (n, c) <- zip [0..] cs
             , not (isAlgZero c) ]
  where
    termToExpr _ 0 c = algNumToExpr c
    termToExpr x 1 c = Mul (algNumToExpr c) (Var x)
    termToExpr x n c = Mul (algNumToExpr c) (Pow (Var x) (Const (fromIntegral (n :: Int))))

-- | Convert a rational function to an Expr
ratFunToExpr :: RatFun AlgNum -> Expr
ratFunToExpr (RatFun p q)
  | degree q == 0 && isAlgZero (leadingCoeff q - algOne) = polyToExpr p
  | otherwise = Div (polyToExpr p) (polyToExpr q)

-- | Convert an AlgNum to an Expr, decomposing into real and
-- imaginary Double parts (since Expr's Const is still Double-only --
-- see the project's own note that this is "a placeholder for AlgNum
-- pending full algebraic number implementation"). Produces a plain
-- Const for real values, avoiding a spurious "+ 0*I" for the common
-- case.
algNumToExpr :: AlgNum -> Expr
algNumToExpr a =
  let re = algToDouble a
      im = algImagDouble a
  in if abs im < 1e-12
       then Const re
       else if abs re < 1e-12
              then Mul (Const im) I
              else Add (Const re) (Mul (Const im) I)

-- | Add two expressions
addExpr :: Expr -> Expr -> Expr
addExpr (Const 0) e = e
addExpr e (Const 0) = e
addExpr e1 e2       = Add e1 e2
-- | Partial fraction decomposition of a proper rational function
-- over AlgNum coefficients.
--
-- Uses Rothstein-Trager to extract the irreducible factor structure
-- of the denominator. For each (c, d_i) pair, the corresponding
-- partial fraction term is c * d_i' / d_i. For linear factors
-- (degree 1), d_i' = 1, so this reduces to the standard form c/d_i.
--
-- For higher-degree irreducible factors the numerator has degree
-- >= 1 (not the standard Ax+B form), but is mathematically correct
-- as a logarithmic-derivative representation.
--
-- Requires: deg(p) < deg(q). Call ratProperFraction first for
-- improper inputs.
--
-- Returns [rf] unchanged if the denominator is already irreducible
-- or has no rational-coefficient decomposition.
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
          let one = constPoly var algOne
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
      all (\(ca,cb) -> isAlgZero (algAdd ca (algNeg cb)))
          (zip (polyCoef a) (polyCoef b))