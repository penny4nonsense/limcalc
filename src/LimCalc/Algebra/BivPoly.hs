-- | Bivariate polynomial arithmetic and resultant computation over ℚ.
--
-- This module serves one primary purpose: computing the minimal
-- polynomial of a sum or product of two algebraic numbers, via
-- resultant construction. If @α@ satisfies @p(x) = 0@ and @β@
-- satisfies @q(x) = 0@, then @α + β@ is a root of
-- @res_y(p(y), q(x − y))@ and @α · β@ is a root of
-- @res_y(p(y), y^n · q(x\/y))@. These resultants are computed here
-- as 'addResultantQ' and 'mulResultantQ'.
--
-- = Representation
--
-- A 'BivPoly' is a univariate polynomial in @y@ whose coefficients
-- are 'QPoly' values (polynomials in @x@). This represents a
-- bivariate polynomial @p(x, y) = ∑ cᵢ(x) · yⁱ@.
--
-- = Determinant method
--
-- The Sylvester matrix is constructed via 'LimCalc.Algebra.Poly.sylvesterMatrix'
-- (which only uses ring operations), but the determinant is computed
-- via cofactor expansion ('cofactorDet') rather than Gaussian
-- elimination. Gaussian elimination requires division by pivots;
-- for 'QPoly' coefficients this means polynomial quotient, which
-- silently drops remainders and is unsound when the division is not
-- exact. Cofactor expansion only uses ring operations and is always
-- correct, at the cost of @O(n!)@ complexity — acceptable since the
-- matrices involved are small (typically 4×4 or smaller).
--
-- = Squarefree reduction
--
-- The raw resultant typically has repeated factors. 'squarefreeRadical'
-- reduces it to its squarefree part (product of distinct irreducible
-- factors) via Yun's algorithm, since 'LimCalc.Algebra.AlgNum.refineToRoot'
-- only needs the root set, not multiplicities.
module LimCalc.Algebra.BivPoly
  ( -- * Type
    BivPoly (..)
    -- * Constructors
  , zeroBiv
  , constBiv
  , monomialBiv
    -- * Properties
  , bivDegree
  , bivLeadingCoeff
    -- * Normalisation
  , bivStrip
    -- * Arithmetic
  , addBiv
  , negBiv
  , subBiv
  , scaleBiv
  , mulBiv
    -- * Substitution
  , qPolyToBiv
  , xMinusY
  , xMinusYPow
  , substituteXMinusY
    -- * Division
  , pseudoDivBiv
  , divModBiv
  , divBivByQPoly
  , qDivPoly
    -- * Subresultant PRS
  , subresultantPRS
    -- * Resultant
  , addResultantBiv
  , mulResultantBiv
  , addResultantQ
  , mulResultantQ
  , sylvesterResultantRing
  , cofactorDet
  , squarefreeRadical
    -- * Poly QPoly views
  , qPolyAsConstInY
  , asPolyInY
  , substXMinusYAsPolyInY
  , substXOverYAsPolyInY
  ) where

import LimCalc.Algebra.QPoly
import qualified LimCalc.Algebra.Poly as Poly

-- | A bivariate polynomial @p(x, y) = ∑ cᵢ(x) · yⁱ@, represented
-- as a univariate polynomial in @y@ with 'QPoly'-in-@x@ coefficients.
--
-- @bivCoef !! i@ is the coefficient of @yⁱ@, as a 'QPoly' in @x@.
-- The zero polynomial is @BivPoly []@. The last element of a
-- non-empty coefficient list is non-zero (maintained by 'bivStrip').
newtype BivPoly = BivPoly { bivCoef :: [QPoly] }
  deriving (Eq, Show)

-- | The zero bivariate polynomial.
zeroBiv :: BivPoly
zeroBiv = BivPoly []

-- | A constant bivariate polynomial (no @y@ dependence).
constBiv :: QPoly -> BivPoly
constBiv p = BivPoly [p]

-- | A monomial @c(x) · yⁿ@.
monomialBiv :: QPoly -> Int -> BivPoly
monomialBiv c n = BivPoly (replicate n (QPoly []) ++ [c])

-- | Degree in @y@. Returns @−1@ for the zero polynomial.
bivDegree :: BivPoly -> Int
bivDegree (BivPoly []) = -1
bivDegree (BivPoly cs) = length cs - 1

-- | Leading coefficient in @y@ (as a 'QPoly' in @x@).
-- Returns the zero polynomial for the zero bivariate polynomial.
bivLeadingCoeff :: BivPoly -> QPoly
bivLeadingCoeff (BivPoly []) = QPoly []
bivLeadingCoeff (BivPoly cs) = last cs

-- | Remove trailing zero-polynomial coefficients.
bivStrip :: BivPoly -> BivPoly
bivStrip (BivPoly cs) =
  BivPoly (reverse $ dropWhile (== QPoly []) $ reverse cs)

-- | Add two bivariate polynomials.
addBiv :: BivPoly -> BivPoly -> BivPoly
addBiv (BivPoly cs1) (BivPoly cs2) =
  bivStrip $ BivPoly $ addCoefs cs1 cs2
  where
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (x:xs) (y:ys) = qAddPoly x y : addCoefs xs ys

-- | Negate a bivariate polynomial.
negBiv :: BivPoly -> BivPoly
negBiv (BivPoly cs) = BivPoly (map negQPoly cs)
  where negQPoly (QPoly xs) = QPoly (map negate xs)

-- | Subtract two bivariate polynomials.
subBiv :: BivPoly -> BivPoly -> BivPoly
subBiv p q = addBiv p (negBiv q)

-- | Scale a bivariate polynomial by a 'QPoly' factor.
scaleBiv :: QPoly -> BivPoly -> BivPoly
scaleBiv c (BivPoly cs) =
  bivStrip $ BivPoly (map (qMulPoly c) cs)

-- | Multiply two bivariate polynomials (Cauchy product in @y@).
mulBiv :: BivPoly -> BivPoly -> BivPoly
mulBiv (BivPoly []) _ = zeroBiv
mulBiv _ (BivPoly []) = zeroBiv
mulBiv (BivPoly cs1) (BivPoly cs2) =
  bivStrip $ BivPoly $ mulCoefs cs1 cs2
  where
    mulCoefs [] _      = []
    mulCoefs (c:cs) ys =
      addCoefs (map (qMulPoly c) ys) (QPoly [] : mulCoefs cs ys)
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (x:xs) (y:ys) = qAddPoly x y : addCoefs xs ys

-- | Embed a 'QPoly' in @y@ as a 'BivPoly' (treating the 'QPoly'
-- coefficients as rationals, each becoming a constant polynomial in @x@).
qPolyToBiv :: QPoly -> BivPoly
qPolyToBiv (QPoly cs) = BivPoly (map (\c -> QPoly [c]) cs)

-- | The bivariate polynomial @x − y@.
xMinusY :: BivPoly
xMinusY = BivPoly [QPoly [0, 1], QPoly [-1]]

-- | The bivariate polynomial @(x − y)^n@.
xMinusYPow :: Int -> BivPoly
xMinusYPow 0 = BivPoly [QPoly [1]]
xMinusYPow n = mulBiv xMinusY (xMinusYPow (n-1))

-- | Substitute @(x − y)@ for the variable of a 'QPoly', giving
-- @q(x − y)@ as a 'BivPoly'.
substituteXMinusY :: QPoly -> BivPoly
substituteXMinusY (QPoly cs) =
  foldr addBiv zeroBiv
    [ scaleBiv (QPoly [c]) (xMinusYPow n)
    | (n, c) <- zip [0..] cs
    , c /= 0
    ]

-- | Pseudo-division of 'BivPoly' values: multiply @a@ by
-- @lc(b)^delta@ before dividing so that the division is exact
-- over 'QPoly' coefficients. Returns @(delta, quotient, remainder)@.
pseudoDivBiv :: BivPoly -> BivPoly -> (Int, BivPoly, BivPoly)
pseudoDivBiv a b
  | bivDegree b < 0            = error "Division by zero BivPoly"
  | bivDegree a < bivDegree b  = (0, zeroBiv, a)
  | otherwise =
      let delta = bivDegree a - bivDegree b + 1
          lcb   = bivLeadingCoeff b
          a'    = scaleBiv (qPow lcb delta) a
          (q,r) = divModBiv a' b
      in (delta, q, r)

-- | Euclidean division of 'BivPoly' values in @y@.
divModBiv :: BivPoly -> BivPoly -> (BivPoly, BivPoly)
divModBiv a b
  | bivDegree a < bivDegree b = (zeroBiv, a)
  | otherwise = go a zeroBiv
  where
    lc = bivLeadingCoeff b
    db = bivDegree b
    go r acc
      | bivDegree r < db = (acc, r)
      | otherwise =
          let scale  = qQuotPoly (bivLeadingCoeff r) lc
              deg    = bivDegree r - db
              term   = monomialBiv scale deg
              r'     = bivStrip $ subBiv r (mulBiv term b)
          in go r' (addBiv acc term)

-- | Exact polynomial quotient, alias for 'qQuotPoly'.
qDivPoly :: QPoly -> QPoly -> QPoly
qDivPoly = qQuotPoly

-- | Subresultant pseudo-remainder sequence of two 'BivPoly' values.
--
-- Returns the full sequence of subresultants, from which the GCD
-- can be extracted. Used internally by the resultant computation.
subresultantPRS :: BivPoly -> BivPoly -> [BivPoly]
subresultantPRS p q
  | bivDegree p < bivDegree q = subresultantPRS q p
  | otherwise = go p q (QPoly [1]) (QPoly [1])
  where
    go a b _beta psi
      | bivDegree b < 0 = [a]
      | otherwise =
          let delta     = bivDegree a - bivDegree b
              (_, _, r) = pseudoDivBiv a b
              lcb       = bivLeadingCoeff b
              psi'      = if delta == 0
                            then psi
                            else qQuotPoly (negQPow lcb delta)
                                           (qPow psi (max 0 (delta-1)))
              beta'     = qMulPoly (negQPoly lcb) (qPow psi' delta)
              r'        = divBivByQPoly r beta'
          in a : go b r' beta' psi'
    negQPow q n  = qScalePoly ((-1)^n) (qPow q n)
    negQPoly (QPoly cs) = QPoly (map negate cs)

-- | Divide every 'QPoly' coefficient of a 'BivPoly' by a 'QPoly'
-- scalar, asserting that every division is exact.
--
-- The subresultant PRS algorithm constructs @beta'@ so that it
-- exactly divides every coefficient of the remainder at each step.
-- A nonzero remainder here indicates an internal algorithmic error.
-- An earlier stub using "invert and multiply" only worked for
-- degree-0 divisors and silently no-op'd otherwise, corrupting the
-- PRS recursion.
divBivByQPoly :: BivPoly -> QPoly -> BivPoly
divBivByQPoly (BivPoly cs) divisor =
  BivPoly [ checkedQuot c divisor | c <- cs ]
  where
    checkedQuot c d =
      let (quot, rem') = qDivModPoly c d
      in if qStrip rem' == QPoly []
           then quot
           else error
             ("divBivByQPoly: inexact division in subresultant PRS \
              \(nonzero remainder) -- this indicates an algorithmic \
              \error, not expected for valid subresultant input")

-- | Embed a 'QPoly' in @x@ as a constant polynomial in @y@
-- (a @'Poly.Poly' 'QPoly'@ with a single coefficient).
qPolyAsConstInY :: QPoly -> Poly.Poly QPoly
qPolyAsConstInY p = Poly.Poly "y" [p]

-- | View a 'QPoly' @p(t)@ as a @'Poly.Poly' 'QPoly'@ in @y@, with
-- each rational coefficient @cᵢ@ becoming the constant polynomial
-- @QPoly [cᵢ]@.
asPolyInY :: QPoly -> Poly.Poly QPoly
asPolyInY (QPoly cs) = Poly.Poly "y" (map (\c -> QPoly [c]) cs)

-- | Substitute @(x − y)@ for the variable of a 'QPoly', returning
-- @q(x − y)@ as a @'Poly.Poly' 'QPoly'@ in @y@ with 'QPoly'-in-@x@
-- coefficients.
substXMinusYAsPolyInY :: QPoly -> Poly.Poly QPoly
substXMinusYAsPolyInY (QPoly cs) =
  foldr Poly.addPoly (Poly.Poly "y" [])
    [ Poly.scalePoly (QPoly [c]) (xMinusYPolyInY n)
    | (n, c) <- zip [0 :: Int ..] cs
    , c /= 0
    ]
  where
    xMinusYPolyInY :: Int -> Poly.Poly QPoly
    xMinusYPolyInY 0 = Poly.Poly "y" [QPoly [1]]
    xMinusYPolyInY n = Poly.mulPoly xMinusY1 (xMinusYPolyInY (n - 1))
    xMinusY1 = Poly.Poly "y" [QPoly [0, 1], QPoly [-1]]

-- | Substitute @x\/y@ for the variable of @q@, returning
-- @y^deg(q) · q(x\/y)@ as a @'Poly.Poly' 'QPoly'@ in @y@.
--
-- For @q(t) = ∑ cₖ tᵏ@ of degree @d@, this is @∑ cₖ · xᵏ · y^(d−k)@,
-- a genuine polynomial in @y@ (no negative powers) with
-- 'QPoly'-in-@x@ coefficients @cₖ · xᵏ@.
substXOverYAsPolyInY :: QPoly -> Poly.Poly QPoly
substXOverYAsPolyInY q@(QPoly cs) =
  let d = qDegree q
      terms = [ (d - k, QPoly (replicate k 0 ++ [c]))
              | (k, c) <- zip [0 :: Int ..] cs
              , c /= 0
              ]
      maxDeg = maximum (0 : map fst terms)
      coeffsByDeg = [ sum [ c | (k, c) <- terms, k == n ] | n <- [0 .. maxDeg] ]
  in Poly.Poly "y" coeffsByDeg

-- | Compute @res_y(p(y),\ q(x − y))@: the minimal polynomial of
-- @α + β@ given minimal polynomials @p@ of @α@ and @q@ of @β@.
--
-- Uses cofactor expansion ('cofactorDet') rather than Gaussian
-- elimination to avoid the unsoundness of polynomial division over
-- 'QPoly' (see module header). The result is reduced to its
-- squarefree radical via 'squarefreeRadical'.
addResultantBiv :: QPoly -> QPoly -> QPoly
addResultantBiv pa pb =
  let pY  = asPolyInY pa
      qY  = substXMinusYAsPolyInY pb
      res = sylvesterResultantRing pY qY
  in squarefreeRadical res

-- | Compute @res_y(p(y),\ y^n · q(x\/y))@: the minimal polynomial of
-- @α · β@ given minimal polynomials @p@ of @α@ and @q@ of @β@.
--
-- See 'addResultantBiv' for the rationale for using cofactor
-- expansion.
mulResultantBiv :: QPoly -> QPoly -> QPoly
mulResultantBiv pa pb =
  let pY  = asPolyInY pa
      qY  = substXOverYAsPolyInY pb
      res = sylvesterResultantRing pY qY
  in squarefreeRadical res

-- | Resultant via cofactor (Laplace) expansion of the Sylvester
-- matrix, sound over any commutative ring.
--
-- Uses only ring operations (@+@, @−@, @×@), never division, making
-- it correct for 'QPoly' coefficients. Gaussian elimination was
-- previously used but is unsound here because it requires exact
-- polynomial division at each pivot step, which 'qQuotPoly' silently
-- approximates by dropping remainders — producing a visibly wrong
-- resultant (e.g. @x^4@ instead of @x^4 + 4x^2@ for the @i + i@
-- case). The @O(n!)@ complexity is acceptable for the small matrices
-- arising from low-degree minimal polynomials.
sylvesterResultantRing :: Poly.Poly QPoly -> Poly.Poly QPoly -> QPoly
sylvesterResultantRing p q
  | Poly.degree p < 0 || Poly.degree q < 0 = QPoly []
  | otherwise =
      let mat = Poly.sylvesterMatrix p q
      in cofactorDet mat

-- | Determinant via cofactor (Laplace) expansion along the first row.
--
-- Sound over any commutative ring: only uses @+@, @−@, @×@. Intended
-- for small matrices (the Sylvester matrices of low-degree
-- polynomials); @O(n!)@ complexity makes it impractical for large @n@.
cofactorDet :: (Num a) => [[a]] -> a
cofactorDet []          = 1
cofactorDet [[x]]       = x
cofactorDet mat@(row:_) =
  sum [ sign j * (row !! j) * cofactorDet (minor mat 0 j)
      | j <- [0 .. length row - 1]
      ]
  where
    sign j = if even j then 1 else -1
    minor m r c =
      [ [ m !! i !! k | k <- [0 .. length (head m) - 1], k /= c ]
      | i <- [0 .. length m - 1]
      , i /= r
      ]

-- | Reduce a 'QPoly' to its squarefree radical — the product of its
-- distinct irreducible factors, each to the first power.
--
-- The raw resultants from 'addResultantBiv' and 'mulResultantBiv'
-- typically have repeated factors. 'squarefreeRadical' removes
-- multiplicities via Yun's algorithm (through 'LimCalc.Poly.squarefree'),
-- leaving the root set unchanged. Multiplicity is irrelevant for
-- 'LimCalc.AlgNum.refineToRoot', which selects a root by proximity.
squarefreeRadical :: QPoly -> QPoly
squarefreeRadical p =
  let asPolyX = Poly.Poly "x" (qPolyCoef p)
      factors = Poly.squarefree asPolyX
  in case factors of
       [] -> p
       _  -> QPoly (Poly.polyCoef (foldr1 Poly.mulPoly (map fst factors)))

-- | Additive resultant with a fast path for degree-1 inputs.
--
-- For degree-1 polynomials @p(y) = y − a@ and @q(y) = y − b@, the
-- minimal polynomial of @a + b@ is simply @y − (a + b)@, computed
-- directly without matrix construction.
addResultantQ :: QPoly -> QPoly -> QPoly
addResultantQ pa pb
  | qDegree pa == 1 && qDegree pb == 1 =
      let a = negate (head (qPolyCoef pa))
          b = negate (head (qPolyCoef pb))
      in QPoly [negate (a+b), 1]
  | otherwise = addResultantBiv pa pb

-- | Multiplicative resultant with a fast path for degree-1 inputs.
--
-- For degree-1 polynomials @p(y) = y − a@ and @q(y) = y − b@, the
-- minimal polynomial of @a · b@ is simply @y − ab@, computed
-- directly without matrix construction.
mulResultantQ :: QPoly -> QPoly -> QPoly
mulResultantQ pa pb
  | qDegree pa == 1 && qDegree pb == 1 =
      let a = negate (head (qPolyCoef pa))
          b = negate (head (qPolyCoef pb))
      in QPoly [negate (a*b), 1]
  | otherwise = mulResultantBiv pa pb