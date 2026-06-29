module LimCalc.Poly where

-- | Univariate polynomial over a coefficient ring a.
-- Coefficients in ascending degree order.
-- Invariant: last coefficient is nonzero (or list is empty).
data Poly a = Poly
  { polyVar  :: String
  , polyCoef :: [a]
  } deriving (Eq)

instance (Show a, Num a, Eq a) => Show (Poly a) where
  show (Poly _ [])  = "0"
  show (Poly x cs)  = concatMap showTerm (reverse $ zip ([0..] :: [Int]) cs)
    where
      showTerm (0, c) = show c
      showTerm (1, c) = show c ++ x
      showTerm (n, c) = show c ++ x ++ "^" ++ show n

-- | The zero polynomial
zeroPoly :: Num a => String -> Poly a
zeroPoly x = Poly x []

-- | The unit polynomial
onePoly :: Num a => String -> Poly a
onePoly x = Poly x [1]

-- | Constant polynomial
constPoly :: (Num a, Eq a) => String -> a -> Poly a
constPoly x c
  | c == 0    = Poly x []
  | otherwise = Poly x [c]

-- | Monomial: c * x^n
monomialPoly :: (Num a, Eq a) => String -> a -> Int -> Poly a
monomialPoly x c n
  | c == 0    = zeroPoly x
  | otherwise = Poly x (replicate n 0 ++ [c])

-- | Degree (-1 for zero polynomial)
degree :: Poly a -> Int
degree (Poly _ []) = -1
degree (Poly _ cs) = length cs - 1

-- | Leading coefficient
leadingCoeff :: Num a => Poly a -> a
leadingCoeff (Poly _ []) = 0
leadingCoeff (Poly _ cs) = last cs

-- | Strip trailing zeros
stripTrailing :: (Num a, Eq a) => Poly a -> Poly a
stripTrailing (Poly x cs) = Poly x (reverse $ dropWhile (== 0) $ reverse cs)

-- | Add two polynomials
addPoly :: (Num a, Eq a) => Poly a -> Poly a -> Poly a
addPoly (Poly x cs1) (Poly _ cs2) =
  stripTrailing $ Poly x $ addCoefs cs1 cs2
  where
    addCoefs [] ys           = ys
    addCoefs xs []           = xs
    addCoefs (a:as) (b:bs)   = (a+b) : addCoefs as bs

-- | Negate a polynomial
negPoly :: Num a => Poly a -> Poly a
negPoly (Poly x cs) = Poly x (map negate cs)

-- | Subtract two polynomials
subPoly :: (Num a, Eq a) => Poly a -> Poly a -> Poly a
subPoly p q = addPoly p (negPoly q)

-- | Multiply two polynomials
mulPoly :: (Num a, Eq a) => Poly a -> Poly a -> Poly a
mulPoly (Poly x []) _            = zeroPoly x
mulPoly _ (Poly x [])            = zeroPoly x
mulPoly (Poly x cs1) (Poly _ cs2) =
  stripTrailing $ Poly x $ mulCoefs cs1 cs2
  where
    mulCoefs [] _      = []
    mulCoefs (c:cs) ys =
      addCoefs (map (*c) ys) (0 : mulCoefs cs ys)
    addCoefs [] ys         = ys
    addCoefs xs []         = xs
    addCoefs (a:as) (b:bs) = (a+b) : addCoefs as bs

-- | Scale a polynomial by a constant
scalePoly :: (Num a, Eq a) => a -> Poly a -> Poly a
scalePoly 0 (Poly x _) = zeroPoly x
scalePoly c (Poly x cs) = stripTrailing $ Poly x (map (*c) cs)

-- | Euclidean division: (quotient, remainder)
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

-- | Quotient
quotPoly :: (Fractional a, Eq a) => Poly a -> Poly a -> Poly a
quotPoly p q = fst (divModPoly p q)

-- | Remainder
remPoly :: (Fractional a, Eq a) => Poly a -> Poly a -> Poly a
remPoly p q = snd (divModPoly p q)

-- | Pseudo-division for exact arithmetic
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

-- | GCD via subresultant algorithm
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

-- | Subresultant GCD
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

-- | Differentiate a polynomial
diffPoly :: (Num a, Eq a) => Poly a -> Poly a
diffPoly (Poly x []) = zeroPoly x
diffPoly (Poly x cs) =
  stripTrailing $ Poly x
    [ fromIntegral n * c
    | (n, c) <- zip [1 :: Int ..] (drop 1 cs)
    ]

-- | Evaluate a polynomial at a point
evalPoly :: Num a => Poly a -> a -> a
evalPoly (Poly _ cs) x = foldr (\c acc -> c + x * acc) 0 cs

-- | Squarefree factorization via Yun's algorithm
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

-- | Resultant via Sylvester matrix determinant
resultant :: (Fractional a, Eq a, Num a) => Poly a -> Poly a -> a
resultant p q
  | degree p < 0 || degree q < 0 = 0
  | otherwise = sylvesterDet p q

-- | Sylvester matrix determinant
sylvesterDet :: (Fractional a, Eq a) => Poly a -> Poly a -> a
sylvesterDet p q =
  let mat = sylvesterMatrix p q
      sz  = degree p + degree q
  in determinant sz mat

-- | Build Sylvester matrix
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

-- | Gaussian elimination determinant
determinant :: (Fractional a, Eq a) => Int -> [[a]] -> a
determinant 0 _   = 1
determinant _ mat = go mat 1
  where
    go [] acc = acc
    go (r:rs) acc =
      case findPivot r rs of
        Nothing     -> 0
        Just (p, rs') ->
          let pivot = head p
              rs''  = map (eliminate pivot p) rs'
          in go (map tail rs'') (acc * pivot)
    findPivot r rs
      | not (head r == 0) = Just (r, rs)
      | otherwise = case break (\r' -> not (head r' == 0)) rs of
          (_, [])           -> Nothing
          (before, (x:after)) -> Just (x, before ++ r : after)
    eliminate pivot pivotRow row =
      let factor = head row / pivot
      in zipWith (\a b -> a - factor * b) row pivotRow