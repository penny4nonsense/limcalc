module LimCalc.Poly where

-- | Univariate polynomial over Double coefficients.
-- Coefficients are stored in ascending degree order:
-- [a₀, a₁, a₂, ...] represents a₀ + a₁x + a₂x² + ...
-- The zero polynomial is represented as [].
-- Invariant: the last coefficient is always nonzero (or list is empty).
data Poly = Poly
  { polyVar  :: String    -- ^ The variable name
  , polyCoef :: [Double]  -- ^ Coefficients in ascending degree order
  } deriving (Eq)

instance Show Poly where
  show (Poly x [])  = "0"
  show (Poly x cs)  = showPoly x cs

showPoly :: String -> [Double] -> String
showPoly x cs = concatMap showTerm (reverse $ zip [0..] cs)
  where
    showTerm (0, c) = show c
    showTerm (1, c) = show c ++ x
    showTerm (n, c) = show c ++ x ++ "^" ++ show n

-- | The zero polynomial
zeroPoly :: String -> Poly
zeroPoly x = Poly x []

-- | The unit polynomial (constant 1)
onePoly :: String -> Poly
onePoly x = Poly x [1.0]

-- | Constant polynomial
constPoly :: String -> Double -> Poly
constPoly x 0 = Poly x []
constPoly x c = Poly x [c]

-- | Monomial: c * x^n
monomialPoly :: String -> Double -> Int -> Poly
monomialPoly x c n
  | c == 0    = zeroPoly x
  | otherwise = Poly x (replicate n 0.0 ++ [c])

-- | Degree of a polynomial (-1 for zero polynomial)
degree :: Poly -> Int
degree (Poly _ []) = -1
degree (Poly _ cs) = length cs - 1

-- | Leading coefficient
leadingCoeff :: Poly -> Double
leadingCoeff (Poly _ []) = 0
leadingCoeff (Poly _ cs) = last cs

-- | Strip trailing zeros to maintain invariant
stripTrailing :: Poly -> Poly
stripTrailing (Poly x cs) = Poly x (reverse $ dropWhile (== 0) $ reverse cs)

-- | Add two polynomials
addPoly :: Poly -> Poly -> Poly
addPoly (Poly x cs1) (Poly _ cs2) =
  stripTrailing $ Poly x $ addCoefs cs1 cs2
  where
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (x:xs) (y:ys) = (x+y) : addCoefs xs ys

-- | Negate a polynomial
negPoly :: Poly -> Poly
negPoly (Poly x cs) = Poly x (map negate cs)

-- | Subtract two polynomials
subPoly :: Poly -> Poly -> Poly
subPoly p q = addPoly p (negPoly q)

-- | Multiply two polynomials
mulPoly :: Poly -> Poly -> Poly
mulPoly (Poly x [])  _           = zeroPoly x
mulPoly _            (Poly x []) = zeroPoly x
mulPoly (Poly x cs1) (Poly _ cs2) =
  stripTrailing $ Poly x $ mulCoefs cs1 cs2
  where
    mulCoefs [] _   = []
    mulCoefs (c:cs) ys =
      addCoefs (map (*c) ys) (0.0 : mulCoefs cs ys)
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (x:xs) (y:ys) = (x+y) : addCoefs xs ys

-- | Scale a polynomial by a constant
scalePoly :: Double -> Poly -> Poly
scalePoly 0 (Poly x _)  = zeroPoly x
scalePoly c (Poly x cs) = stripTrailing $ Poly x (map (*c) cs)

-- | Euclidean division: (quotient, remainder)
-- p = q * quot + rem, degree(rem) < degree(q)
divModPoly :: Poly -> Poly -> (Poly, Poly)
divModPoly p q
  | degree q < 0  = error "Division by zero polynomial"
  | degree p < degree q = (zeroPoly (polyVar p), p)
  | otherwise = go p (zeroPoly (polyVar p))
  where
    lc = leadingCoeff q
    dq = degree q
    go r acc
      | degree r < dq = (acc, r)
      | otherwise =
          let scale  = leadingCoeff r / lc
              deg    = degree r - dq
              term   = monomialPoly (polyVar p) scale deg
              r'     = stripTrailing $ subPoly r (mulPoly term q)
          in go r' (addPoly acc term)

-- | Quotient of polynomial division
quotPoly :: Poly -> Poly -> Poly
quotPoly p q = fst (divModPoly p q)

-- | Remainder of polynomial division
remPoly :: Poly -> Poly -> Poly
remPoly p q = snd (divModPoly p q)

-- | Pseudo-division for exact arithmetic
-- Returns (k, q, r) such that lc(b)^k * a = b*q + r
pseudoDivMod :: Poly -> Poly -> (Int, Poly, Poly)
pseudoDivMod a b
  | degree b < 0 = error "Pseudo-division by zero"
  | degree a < degree b = (0, zeroPoly (polyVar a), a)
  | otherwise =
      let delta = degree a - degree b + 1
          lcb   = leadingCoeff b
          a'    = scalePoly (lcb ^ delta) a
          (q,r) = divModPoly a' b
      in (delta, q, r)

-- | GCD via subresultant algorithm
-- Returns monic GCD
gcdPoly :: Poly -> Poly -> Poly
gcdPoly p q
  | degree p < 0 = q
  | degree q < 0 = p
  | degree p < degree q = gcdPoly q p
  | otherwise = makeMonic $ subresultantGCD p q
  where
    makeMonic r
      | degree r < 0 = onePoly (polyVar p)
      | otherwise    = scalePoly (1 / leadingCoeff r) r

-- | Subresultant GCD algorithm
subresultantGCD :: Poly -> Poly -> Poly
subresultantGCD p q = go p q 1 1
  where
    go a b beta psi
      | degree b < 0 = a
      | otherwise =
          let delta  = degree a - degree b
              (_, _, r) = pseudoDivMod a b
              psi'   = if delta == 0
                         then psi
                         else ((-leadingCoeff b) ^ delta) /
                              (psi ^ (delta - 1))
              beta'  = (-leadingCoeff b) * psi' ^ delta
              r'     = scalePoly (1/beta') r
          in go b r' beta' psi'

-- | Differentiate a polynomial
diffPoly :: Poly -> Poly
diffPoly (Poly x [])  = zeroPoly x
diffPoly (Poly x cs)  =
  stripTrailing $ Poly x
    [ fromIntegral n * c
    | (n, c) <- zip [1..] (tail cs)
    ]

-- | Evaluate a polynomial at a point
evalPoly :: Poly -> Double -> Double
evalPoly (Poly _ cs) x = foldr (\c acc -> c + x * acc) 0 cs

-- | Squarefree factorization via Yun's algorithm
-- Returns list of (factor, multiplicity) pairs
squarefree :: Poly -> [(Poly, Int)]
squarefree p
  | degree p < 0 = []
  | otherwise =
      let p'  = diffPoly p
          g   = gcdPoly p p'
          p1  = quotPoly p g
          p1' = diffPoly p1
          h   = subPoly p1' (quotPoly (mulPoly p' g) (mulPoly g g))
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

-- | Resultant of two polynomials via subresultant sequence
resultant :: Poly -> Poly -> Double
resultant p q
  | degree p < 0 || degree q < 0 = 0
  | degree q > degree p = 
      let sign = if odd (degree p * degree q) then -1 else 1
      in sign * resultant q p
  | otherwise = go p q 1
  where
    go a b acc
      | degree b < 0 = acc * leadingCoeff a ^ degree b
      | degree b == 0 = acc * leadingCoeff b ^ degree a
      | otherwise =
          let (_, _, r) = pseudoDivMod a b
              delta     = degree a - degree b
              sign      = if odd (degree a * degree b) then -1 else 1
              lc        = leadingCoeff b ^ (degree a - degree r)
          in go b r (acc * sign * lc)