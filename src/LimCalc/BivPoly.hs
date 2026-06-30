module LimCalc.BivPoly where

import LimCalc.QPoly
import qualified LimCalc.Poly as Poly

-- | Bivariate polynomial in y with QPoly-in-x coefficients
-- BivPoly [c0, c1, c2, ...] represents c0 + c1*y + c2*y^2 + ...
-- where each ci is a QPoly in x
newtype BivPoly = BivPoly { bivCoef :: [QPoly] }
  deriving (Eq, Show)

-- | Zero bivariate polynomial
zeroBiv :: BivPoly
zeroBiv = BivPoly []

-- | Constant bivariate polynomial (no y dependence)
constBiv :: QPoly -> BivPoly
constBiv p = BivPoly [p]

-- | Monomial: c * y^n where c is a QPoly
monomialBiv :: QPoly -> Int -> BivPoly
monomialBiv c n = BivPoly (replicate n (QPoly []) ++ [c])

-- | Degree in y
bivDegree :: BivPoly -> Int
bivDegree (BivPoly []) = -1
bivDegree (BivPoly cs) = length cs - 1

-- | Leading coefficient (highest power of y) as QPoly
bivLeadingCoeff :: BivPoly -> QPoly
bivLeadingCoeff (BivPoly []) = QPoly []
bivLeadingCoeff (BivPoly cs) = last cs

-- | Strip trailing zero polynomials
bivStrip :: BivPoly -> BivPoly
bivStrip (BivPoly cs) =
  BivPoly (reverse $ dropWhile (== QPoly []) $ reverse cs)

-- | Add two bivariate polynomials
addBiv :: BivPoly -> BivPoly -> BivPoly
addBiv (BivPoly cs1) (BivPoly cs2) =
  bivStrip $ BivPoly $ addCoefs cs1 cs2
  where
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (x:xs) (y:ys) = qAddPoly x y : addCoefs xs ys

-- | Negate a bivariate polynomial
negBiv :: BivPoly -> BivPoly
negBiv (BivPoly cs) = BivPoly (map negQPoly cs)
  where negQPoly (QPoly xs) = QPoly (map negate xs)

-- | Subtract two bivariate polynomials
subBiv :: BivPoly -> BivPoly -> BivPoly
subBiv p q = addBiv p (negBiv q)

-- | Scale by a QPoly
scaleBiv :: QPoly -> BivPoly -> BivPoly
scaleBiv c (BivPoly cs) =
  bivStrip $ BivPoly (map (qMulPoly c) cs)

-- | Multiply two bivariate polynomials
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

-- | Convert a QPoly p(y) to a BivPoly
qPolyToBiv :: QPoly -> BivPoly
qPolyToBiv (QPoly cs) = BivPoly (map (\c -> QPoly [c]) cs)

-- | The bivariate polynomial (x - y)
xMinusY :: BivPoly
xMinusY = BivPoly [QPoly [0, 1], QPoly [-1]]

-- | Compute (x - y)^n as a BivPoly
xMinusYPow :: Int -> BivPoly
xMinusYPow 0 = BivPoly [QPoly [1]]
xMinusYPow n = mulBiv xMinusY (xMinusYPow (n-1))

-- | Substitute (x-y) for the variable in q(t) to get q(x-y) as BivPoly
substituteXMinusY :: QPoly -> BivPoly
substituteXMinusY (QPoly cs) =
  foldr addBiv zeroBiv
    [ scaleBiv (QPoly [c]) (xMinusYPow n)
    | (n, c) <- zip [0..] cs
    , c /= 0
    ]

-- | Pseudo-division of BivPolys over QPoly coefficients
pseudoDivBiv :: BivPoly -> BivPoly -> (Int, BivPoly, BivPoly)
pseudoDivBiv a b
  | bivDegree b < 0 = error "Division by zero BivPoly"
  | bivDegree a < bivDegree b = (0, zeroBiv, a)
  | otherwise =
      let delta = bivDegree a - bivDegree b + 1
          lcb   = bivLeadingCoeff b
          a'    = scaleBiv (qPow lcb delta) a
          (q,r) = divModBiv a' b
      in (delta, q, r)

-- | Euclidean division of BivPolys
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

-- | Divide two QPolys exactly (real division, see LimCalc.QPoly.qDivModPoly)
qDivPoly :: QPoly -> QPoly -> QPoly
qDivPoly = qQuotPoly

-- | Subresultant PRS of two BivPolys
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

-- | Divide every coefficient of a BivPoly by a QPoly scalar.
--
-- This division is mathematically guaranteed exact by the
-- subresultant PRS algorithm (beta' is constructed so that it
-- exactly divides every coefficient of r at each step). A nonzero
-- remainder indicates an internal algorithmic error -- previously
-- this used a stubbed "invert and multiply" approach (qInvPoly) that
-- only worked when beta' was degree 0 and silently no-op'd
-- otherwise, corrupting the PRS recursion for any case where beta'
-- was a genuine non-constant polynomial in x (which is the typical
-- case once BivPoly coefficients are themselves QPolys in x, as they
-- are throughout this module).
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

-- | Convert a QPoly (polynomial in x) into a Poly QPoly with a
-- single constant coefficient -- i.e. embed it as a "constant in y"
-- bivariate polynomial, viewed as a univariate polynomial in y over
-- the coefficient ring QPoly (polynomials in x).
qPolyAsConstInY :: QPoly -> Poly.Poly QPoly
qPolyAsConstInY p = Poly.Poly "y" [p]

-- | p(y), the bivariate polynomial viewed as a Poly QPoly in y with
-- QPoly-in-x coefficients pulled directly from p's own coefficients
-- (each coefficient ci of p becomes the constant QPoly [ci]).
asPolyInY :: QPoly -> Poly.Poly QPoly
asPolyInY (QPoly cs) = Poly.Poly "y" (map (\c -> QPoly [c]) cs)

-- | q(x - y) as a Poly QPoly in y: substitute (x - y) for the
-- variable of q, giving a polynomial in y whose coefficients are
-- QPolys in x.
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
    xMinusY1 = Poly.Poly "y" [QPoly [0, 1], QPoly [-1]]  -- x + (-1)*y

-- | y^deg(q) * q(x/y) as a Poly QPoly in y -- the multiplicative
-- substitution. For q(t) = sum c_k t^k (degree d), this is
-- sum c_k * x^k * y^(d-k), a genuine polynomial in y (no negative
-- powers), with QPoly-in-x coefficients c_k * x^k.
substXOverYAsPolyInY :: QPoly -> Poly.Poly QPoly
substXOverYAsPolyInY q@(QPoly cs) =
  let d = qDegree q
      terms = [ (d - k, QPoly (replicate k 0 ++ [c]))  -- c * x^k, as a QPoly in x
              | (k, c) <- zip [0 :: Int ..] cs
              , c /= 0
              ]
      maxDeg = maximum (0 : map fst terms)
      coeffsByDeg = [ sum [ c | (k, c) <- terms, k == n ] | n <- [0 .. maxDeg] ]
  in Poly.Poly "y" coeffsByDeg

-- | Resultant of p(y) and q(x - y) over y -- the additive case.
--
-- Computes the resultant via cofactor expansion of the Sylvester
-- matrix (cofactorDet), NOT Poly.hs's Gaussian-elimination
-- determinant. Gaussian elimination divides by pivots, which for
-- QPoly coefficients means polynomial quotient (qQuotPoly) -- this
-- silently drops any remainder, which is mathematically unsound
-- whenever an elimination step's division isn't exact. This
-- produced a visibly wrong resultant (x^4 instead of the correct
-- x^4 + 4x^2 for the i+i case). Cofactor expansion only needs ring
-- operations (add, multiply, negate), never division, so it's sound
-- over QPoly even though QPoly isn't a true field. The O(n!) cost is
-- acceptable here since the polynomials involved are low degree.
addResultantBiv :: QPoly -> QPoly -> QPoly
addResultantBiv pa pb =
  let pY  = asPolyInY pa
      qY  = substXMinusYAsPolyInY pb
      res = sylvesterResultantRing pY qY
  in squarefreeRadical res

-- | Resultant of p(y) and y^deg(q) * q(x/y) over y -- the
-- multiplicative case. See addResultantBiv for why cofactor
-- expansion is used instead of Poly.hs's Gaussian-elimination
-- resultant.
mulResultantBiv :: QPoly -> QPoly -> QPoly
mulResultantBiv pa pb =
  let pY  = asPolyInY pa
      qY  = substXOverYAsPolyInY pb
      res = sylvesterResultantRing pY qY
  in squarefreeRadical res

-- | Resultant via cofactor expansion of the Sylvester matrix, sound
-- over any commutative ring (only uses +, -, * -- never division).
-- Reuses Poly.hs's sylvesterMatrix construction (that part is purely
-- additive/multiplicative bookkeeping, not division-based, so it's
-- safe to reuse), but computes the determinant via cofactor
-- expansion instead of Poly.hs's Gaussian-elimination determinant.
sylvesterResultantRing :: Poly.Poly QPoly -> Poly.Poly QPoly -> QPoly
sylvesterResultantRing p q
  | Poly.degree p < 0 || Poly.degree q < 0 = QPoly []
  | otherwise =
      let mat = Poly.sylvesterMatrix p q
      in cofactorDet mat

-- | Determinant via cofactor (Laplace) expansion along the first
-- row. Sound over any commutative ring: only uses ring operations,
-- never division. O(n!) -- intended for the small matrices that
-- arise from resultants of low-degree polynomials, not general use.
cofactorDet :: (Num a) => [[a]] -> a
cofactorDet []          = 1  -- determinant of the empty matrix
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

-- | Reduce a QPoly to its squarefree radical (product of distinct
-- irreducible factors, each to the first power), via Poly.hs's
-- already-correct Yun's-algorithm squarefree factorization. The raw
-- resultant from addResultantBiv/mulResultantBiv generally has
-- repeated factors; this removes the multiplicity while leaving the
-- root SET unchanged, which is what the caller's root-finding step
-- needs (it picks one specific root by proximity, so multiplicity is
-- irrelevant once a root is chosen).
squarefreeRadical :: QPoly -> QPoly
squarefreeRadical p =
  let asPolyX = Poly.Poly "x" (qPolyCoef p)
      factors = Poly.squarefree asPolyX
  in case factors of
       [] -> p  -- already squarefree (or zero/degenerate)
       _  -> QPoly (Poly.polyCoef (foldr1 Poly.mulPoly (map fst factors)))

-- | Additive resultant: res_y(p(y), q(x-y))
addResultantQ :: QPoly -> QPoly -> QPoly
addResultantQ pa pb
  | qDegree pa == 1 && qDegree pb == 1 =
      let a = negate (head (qPolyCoef pa))
          b = negate (head (qPolyCoef pb))
      in QPoly [negate (a+b), 1]
  | otherwise = addResultantBiv pa pb

-- | Multiplicative resultant: res_y(p(y), y^n * q(x/y))
mulResultantQ :: QPoly -> QPoly -> QPoly
mulResultantQ pa pb
  | qDegree pa == 1 && qDegree pb == 1 =
      -- p(y) = y - a, q(y) = y - b
      -- minimal poly of a*b is y - a*b
      let a = negate (head (qPolyCoef pa))
          b = negate (head (qPolyCoef pb))
      in QPoly [negate (a*b), 1]
  | otherwise = mulResultantBiv pa pb