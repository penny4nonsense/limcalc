-- | Univariate polynomials over an arbitrary coefficient ring.
--
-- 'Poly' is used throughout limcalc as the coefficient ring for
-- rational functions ('LimCalc.Algebra.RationalFunction'), as the
-- input to resultant computation ('LimCalc.Algebra.BivPoly'), and as
-- the polynomial arithmetic layer underlying the Risch integrator
-- ('LimCalc.Integration.Risch.Primitive').
--
-- = Representation
--
-- Coefficients are stored in /ascending degree order/: index @i@
-- holds the coefficient of @x^i@. The zero polynomial is represented
-- by an empty coefficient list. The invariant is that the last
-- element of a non-empty coefficient list is nonzero (maintained by
-- 'stripTrailing').
--
-- = GCD and resultant
--
-- GCD is computed via the subresultant pseudo-remainder sequence
-- ('subresultantGCD'), which avoids coefficient explosion over
-- general integral domains. The resultant is computed via the
-- Sylvester matrix determinant with Gaussian elimination
-- ('determinant'), which correctly accounts for row-swap sign flips.
module LimCalc.Algebra.Poly
  ( -- * Type
    Poly (..)
    -- * Constructors
  , zeroPoly
  , onePoly
  , constPoly
  , monomialPoly
    -- * Properties
  , degree
  , leadingCoeff
    -- * Normalisation
  , stripTrailing
    -- * Arithmetic
  , addPoly
  , negPoly
  , subPoly
  , mulPoly
  , scalePoly
    -- * Division
  , divModPoly
  , quotPoly
  , remPoly
  , pseudoDivMod
    -- * GCD
  , gcdPoly
  , subresultantGCD
    -- * Calculus
  , diffPoly
  , evalPoly
    -- * Factorisation
  , squarefree
    -- * Resultant
  , resultant
  , sylvesterDet
  , sylvesterMatrix
  , determinant
  ) where

-- | Univariate polynomial over a coefficient ring @a@.
--
-- Coefficients are in ascending degree order: @polyCoef !! i@ is the
-- coefficient of @x^i@. The zero polynomial has @polyCoef = []@.
-- All operations maintain the invariant that the last coefficient is
-- nonzero (via 'stripTrailing').
data Poly a = Poly
  { polyVar  :: String
    -- ^ The name of the indeterminate (used for display only).
  , polyCoef :: [a]
    -- ^ Coefficients in ascending degree order.
  } deriving (Eq)

instance (Show a, Num a, Eq a) => Show (Poly a) where
  show (Poly _ [])  = "0"
  show (Poly x cs)  = concatMap showTerm (reverse $ zip ([0..] :: [Int]) cs)
    where
      showTerm (0, c) = show c
      showTerm (1, c) = show c ++ x
      showTerm (n, c) = show c ++ x ++ "^" ++ show n

-- | The zero polynomial (empty coefficient list).
zeroPoly :: Num a => String -> Poly a
zeroPoly x = Poly x []

-- | The unit polynomial @1@.
onePoly :: Num a => String -> Poly a
onePoly x = Poly x [1]

-- | A constant polynomial @c@ (or the zero polynomial if @c = 0@).
constPoly :: (Num a, Eq a) => String -> a -> Poly a
constPoly x c
  | c == 0    = Poly x []
  | otherwise = Poly x [c]

-- | A monomial @c · x^n@ (or the zero polynomial if @c = 0@).
monomialPoly :: (Num a, Eq a) => String -> a -> Int -> Poly a
monomialPoly x c n
  | c == 0    = zeroPoly x
  | otherwise = Poly x (replicate n 0 ++ [c])

-- | Degree of a polynomial. Returns @−1@ for the zero polynomial.
degree :: Poly a -> Int
degree (Poly _ []) = -1
degree (Poly _ cs) = length cs - 1

-- | Leading coefficient (coefficient of the highest-degree term).
-- Returns @0@ for the zero polynomial.
leadingCoeff :: Num a => Poly a -> a
leadingCoeff (Poly _ []) = 0
leadingCoeff (Poly _ cs) = last cs

-- | Remove trailing zero coefficients to restore the canonical form.
stripTrailing :: (Num a, Eq a) => Poly a -> Poly a
stripTrailing (Poly x cs) = Poly x (reverse $ dropWhile (== 0) $ reverse cs)

-- | Add two polynomials.
addPoly :: (Num a, Eq a) => Poly a -> Poly a -> Poly a
addPoly (Poly x cs1) (Poly _ cs2) =
  stripTrailing $ Poly x $ addCoefs cs1 cs2
  where
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (a:as) (b:bs) = (a+b) : addCoefs as bs

-- | Negate a polynomial (negate every coefficient).
negPoly :: Num a => Poly a -> Poly a
negPoly (Poly x cs) = Poly x (map negate cs)

-- | Subtract two polynomials.
subPoly :: (Num a, Eq a) => Poly a -> Poly a -> Poly a
subPoly p q = addPoly p (negPoly q)

-- | Multiply two polynomials (Cauchy product of coefficient lists).
mulPoly :: (Num a, Eq a) => Poly a -> Poly a -> Poly a
mulPoly (Poly x []) _             = zeroPoly x
mulPoly _ (Poly x [])             = zeroPoly x
mulPoly (Poly x cs1) (Poly _ cs2) =
  stripTrailing $ Poly x $ mulCoefs cs1 cs2
  where
    mulCoefs [] _      = []
    mulCoefs (c:cs) ys =
      addCoefs (map (*c) ys) (0 : mulCoefs cs ys)
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (a:as) (b:bs) = (a+b) : addCoefs as bs

-- | Scale every coefficient by a constant.
scalePoly :: (Num a, Eq a) => a -> Poly a -> Poly a
scalePoly 0 (Poly x _) = zeroPoly x
scalePoly c (Poly x cs) = stripTrailing $ Poly x (map (*c) cs)

-- | Euclidean division: @divModPoly p q = (quot, rem)@ where
-- @p = quot * q + rem@ and @degree rem < degree q@.
--
-- Requires a 'Fractional' coefficient ring so that leading
-- coefficients can be divided exactly.
divModPoly :: (Fractional a, Eq a) => Poly a -> Poly a -> (Poly a, Poly a)
divModPoly p q
  | degree q < 0        = error "Division by zero polynomial"
  | degree p < degree q = (zeroPoly (polyVar p), p)
  | otherwise           = go p (zeroPoly (polyVar p))
  where
    lc = leadingCoeff q
    dq = degree q
    go r acc
      | degree r < dq = (acc, r)
      | otherwise =
          let scale = leadingCoeff r / lc
              deg   = degree r - dq
              term  = monomialPoly (polyVar p) scale deg
              r'    = stripTrailing $ subPoly r (mulPoly term q)
          in go r' (addPoly acc term)

-- | Quotient of Euclidean division.
quotPoly :: (Fractional a, Eq a) => Poly a -> Poly a -> Poly a
quotPoly p q = fst (divModPoly p q)

-- | Remainder of Euclidean division.
remPoly :: (Fractional a, Eq a) => Poly a -> Poly a -> Poly a
remPoly p q = snd (divModPoly p q)

-- | Pseudo-division: multiply @a@ by @lc(b)^delta@ before dividing,
-- where @delta = deg(a) - deg(b) + 1@, to ensure exact integer-ring
-- arithmetic without fractions.
--
-- Returns @(delta, quotient, remainder)@.
pseudoDivMod :: (Fractional a, Eq a) => Poly a -> Poly a -> (Int, Poly a, Poly a)
pseudoDivMod a b
  | degree b < 0        = error "Pseudo-division by zero"
  | degree a < degree b = (0, zeroPoly (polyVar a), a)
  | otherwise =
      let delta = degree a - degree b + 1
          lcb   = leadingCoeff b
          a'    = scalePoly (lcb ^ delta) a
          (q,r) = divModPoly a' b
      in (delta, q, r)

-- | GCD via the subresultant pseudo-remainder sequence.
--
-- Uses 'subresultantGCD' internally to avoid coefficient explosion,
-- then rescales the result to be monic. Returns @1@ if both inputs
-- are coprime.
gcdPoly :: (Fractional a, Eq a) => Poly a -> Poly a -> Poly a
gcdPoly p q
  | degree p < 0        = q
  | degree q < 0        = p
  | degree p < degree q = gcdPoly q p
  | otherwise           = makeMonic $ subresultantGCD p q
  where
    makeMonic r
      | degree r < 0 = onePoly (polyVar p)
      | otherwise    = scalePoly (1 / leadingCoeff r) r

-- | Subresultant pseudo-remainder sequence for GCD computation.
--
-- Avoids the coefficient explosion of the naive Euclidean algorithm
-- by tracking subresultant multipliers @psi@ and @beta@.
subresultantGCD :: (Fractional a, Eq a) => Poly a -> Poly a -> Poly a
subresultantGCD p q = go p q 1 1
  where
    go a b _beta psi
      | degree b < 0 = a
      | otherwise =
          let delta  = degree a - degree b
              (_, _, r) = pseudoDivMod a b
              lcb    = leadingCoeff b
              psi'   = if delta == 0
                         then psi
                         else ((-lcb) ^ delta) /
                              (psi ^ max 0 (delta-1))
              beta'  = (-lcb) * psi' ^ delta
              r'     = scalePoly (1/beta') r
          in go b r' beta' psi'

-- | Formal derivative of a polynomial.
diffPoly :: (Num a, Eq a) => Poly a -> Poly a
diffPoly (Poly x []) = zeroPoly x
diffPoly (Poly x cs) =
  stripTrailing $ Poly x
    [ fromIntegral n * c
    | (n, c) <- zip [1 :: Int ..] (drop 1 cs)
    ]

-- | Evaluate a polynomial at a point using Horner's method.
evalPoly :: Num a => Poly a -> a -> a
evalPoly (Poly _ cs) x = foldr (\c acc -> c + x * acc) 0 cs

-- | Squarefree factorisation via Yun's algorithm.
--
-- Returns a list of @(factor, multiplicity)@ pairs such that
-- @p = product [ factor ^ multiplicity ]@, with each @factor@
-- squarefree and the factors pairwise coprime.
squarefree :: (Fractional a, Eq a) => Poly a -> [(Poly a, Int)]
squarefree p
  | degree p < 0 = []
  | otherwise =
      let p'  = diffPoly p
          g   = gcdPoly p p'
          p1  = quotPoly p g
          p1' = diffPoly p1
      in yun p1 (gcdPoly p1 p1') 1
  where
    yun a b k
      | degree a <= 0 = []
      | otherwise =
          let c  = gcdPoly a b
              a' = quotPoly a c
              b' = subPoly (quotPoly b c) (diffPoly a')
          in if degree a' > 0
             then (a', k) : yun c (gcdPoly c b') (k+1)
             else yun c (gcdPoly c b') (k+1)

-- | Resultant of two polynomials, computed as the determinant of
-- the Sylvester matrix.
--
-- Returns @0@ if either polynomial is zero.
resultant :: (Fractional a, Eq a, Num a) => Poly a -> Poly a -> a
resultant p q
  | degree p < 0 || degree q < 0 = 0
  | otherwise = sylvesterDet p q

-- | Determinant of the Sylvester matrix of two polynomials.
sylvesterDet :: (Fractional a, Eq a) => Poly a -> Poly a -> a
sylvesterDet p q =
  let mat = sylvesterMatrix p q
      sz  = degree p + degree q
  in determinant sz mat

-- | Sylvester matrix of two polynomials.
--
-- An @(m+n) × (m+n)@ matrix where @m = deg(p)@, @n = deg(q)@,
-- whose determinant is the resultant of @p@ and @q@.
sylvesterMatrix :: Num a => Poly a -> Poly a -> [[a]]
sylvesterMatrix p q =
  let m    = degree p
      n    = degree q
      pc   = reverse (polyCoef p)
      qc   = reverse (polyCoef q)
      pRows = [ replicate i 0 ++ pc ++ replicate (n - 1 - i) 0
              | i <- [0..n-1] ]
      qRows = [ replicate i 0 ++ qc ++ replicate (m - 1 - i) 0
              | i <- [0..m-1] ]
  in pRows ++ qRows

-- | Determinant via Gaussian elimination with partial pivoting.
--
-- Tracks the parity of row swaps and applies a sign correction to
-- the result. An earlier version found non-zero pivots by swapping
-- rows but never corrected for the sign flip this introduces —
-- confirmed as a real bug: @determinant 2 [[0,1],[1,0]]@ returned
-- @1@ instead of @−1@ even for plain 'Double' coefficients.
determinant :: (Fractional a, Eq a) => Int -> [[a]] -> a
determinant 0 _   = 1
determinant _ mat = go mat 1 False
  where
    go [] acc swapped = if swapped then negate acc else acc
    go (r:rs) acc swapped =
      case findPivot r rs of
        Nothing     -> 0
        Just (p, rs', didSwap) ->
          let pivot = head p
              rs''  = map (eliminate pivot p) rs'
          in go (map tail rs'') (acc * pivot) (swapped /= didSwap)
    findPivot r rs
      | not (head r == 0) = Just (r, rs, False)
      | otherwise = case break (\r' -> not (head r' == 0)) rs of
          (_, [])             -> Nothing
          (before, (x:after)) -> Just (x, before ++ r : after, True)
    eliminate pivot pivotRow row =
      let factor = head row / pivot
      in zipWith (\a b -> a - factor * b) row pivotRow