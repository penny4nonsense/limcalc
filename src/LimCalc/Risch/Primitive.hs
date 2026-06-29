module LimCalc.Risch.Primitive where

import LimCalc.Poly
import LimCalc.RationalFunction
import LimCalc.DiffField
import LimCalc.Expr

-- | Result of the Risch algorithm for the primitive case
data PrimitiveResult
  = PrimitiveElementary Expr          -- ^ Elementary antiderivative
  | PrimitiveNonElementary            -- ^ Provably non-elementary
  | PrimitiveError String             -- ^ Algorithm error
  deriving (Show, Eq)

-- | Integrate a rational function in the primitive case
-- Field extension: F(θ) where θ = log(f), Dθ = Df/f
-- Input: rational function p/q in F(θ), the differential field
integratePrimitive :: RatFun -> DiffField -> PrimitiveResult
integratePrimitive rf field =
  let (g, reduced) = hermiteReduce rf
      rtResult = rothsteinTrager reduced field
  in case rtResult of
       Left "NonElementary" -> PrimitiveNonElementary
       Left err             -> PrimitiveError err
       Right terms ->
         let logPart = foldr addExpr (Const 0)
               [ Mul (Const c) (Log (polyToExpr u))
               | (c, u) <- terms
               ]
             gExpr = ratFunToExpr g
         in PrimitiveElementary (addExpr gExpr logPart)

-- | Rothstein-Trager algorithm
-- Given a/d with d squarefree, compute the logarithmic part
-- Returns list of (constant, polynomial) pairs for Σ c·log(u)
rothsteinTrager :: RatFun -> DiffField -> Either String [(Double, Poly)]
rothsteinTrager (RatFun a d) field =
  let d'    = diffPoly d
      rPoly = resultantPoly a d d'
  in case findRationalRoots rPoly of
       Nothing    -> Left "NonElementary"  -- signal non-elementary
       Just roots ->
         Right [ (c, gcdPoly d (subPoly a (scalePoly c d')))
               | c <- roots
               , c /= 0
               ]

-- | Compute the Rothstein-Trager resultant polynomial R(z)
-- R(z) = res_x(d(x), a(x) - z·d'(x)) as a polynomial in z
-- This is computed by treating z as a parameter
resultantPoly :: Poly -> Poly -> Poly -> Poly
resultantPoly a d d' =
  -- For degree n denominator, R(z) has degree n
  -- We evaluate R at enough points and interpolate
  let n      = degree d
      points = [ fromIntegral i | i <- [-n..n] ]
      values = [ resultant d (subPoly a (scalePoly z d')) | z <- points ]
  in interpolate (polyVar d) (zip points values)

-- | Lagrange interpolation
-- Given (x_i, y_i) pairs, find the polynomial passing through them
interpolate :: String -> [(Double, Double)] -> Poly
interpolate x points =
  foldr (addPoly . lagrangeBasis x points) (zeroPoly x)
        (zip [0..] points)

-- | One Lagrange basis polynomial
lagrangeBasis :: String -> [(Double, Double)] -> (Int, (Double, Double)) -> Poly
lagrangeBasis x points (i, (xi, yi)) =
  let others = [ (j, xj) | (j, (xj, _)) <- zip [0..] points, j /= i ]
      basis  = foldr mulBasis (constPoly x 1) others
      scale  = yi / evalPoly basis xi
  in scalePoly scale basis
  where
    mulBasis (_, xj) p =
      mulPoly p (subPoly (Poly x [0, 1]) (constPoly x xj))

-- | Find rational roots of a polynomial
-- Uses rational root theorem: tries p/q for factors of constant/leading term
findRationalRoots :: Poly -> Maybe [Double]
findRationalRoots p
  | degree p < 0 = Just []
  | otherwise    =
      let candidates = rationalRootCandidates p
          roots      = filter (\r -> abs (evalPoly p r) < 1e-10) candidates
      in if null roots && degree p > 0
         then Nothing  -- has roots but none rational
         else Just roots

-- | Generate rational root candidates
rationalRootCandidates :: Poly -> [Double]
rationalRootCandidates (Poly _ []) = []
rationalRootCandidates (Poly _ cs) =
  let c0 = abs (head cs)   -- constant term
      cn = abs (last cs)   -- leading coeff
      -- factors of c0 and cn up to reasonable bound
      factors n = [ fromIntegral i | i <- [1..max 10 (round n)]
                  , round n `mod` i == 0 ]
      ps = factors c0
      qs = factors cn
  in nub [ p/q | p <- 0:ps ++ map negate ps
               , q <- qs
               , q /= 0 ]
  where
    nub [] = []
    nub (x:xs) = x : nub (filter (/= x) xs)

-- | Convert a polynomial to an Expr
polyToExpr :: Poly -> Expr
polyToExpr (Poly x [])  = Const 0
polyToExpr (Poly x cs)  =
  foldr1 Add [ termToExpr x n c
             | (n, c) <- zip [0..] cs
             , c /= 0 ]
  where
    termToExpr x 0 c = Const c
    termToExpr x 1 c = Mul (Const c) (Var x)
    termToExpr x n c = Mul (Const c) (Pow (Var x) (Const (fromIntegral n)))

-- | Convert a rational function to an Expr
ratFunToExpr :: RatFun -> Expr
ratFunToExpr (RatFun p q)
  | degree q == 0 && leadingCoeff q == 1 = polyToExpr p
  | otherwise = Div (polyToExpr p) (polyToExpr q)

-- | Add two expressions, simplifying trivially
addExpr :: Expr -> Expr -> Expr
addExpr (Const 0) e = e
addExpr e (Const 0) = e
addExpr e1 e2       = Add e1 e2