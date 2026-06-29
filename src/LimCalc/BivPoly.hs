module LimCalc.BivPoly where

import LimCalc.QPoly

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
          let scale  = qDivPoly (bivLeadingCoeff r) lc
              deg    = bivDegree r - db
              term   = monomialBiv scale deg
              r'     = bivStrip $ subBiv r (mulBiv term b)
          in go r' (addBiv acc term)

-- | Divide two QPolys (exact division for constants)
qDivPoly :: QPoly -> QPoly -> QPoly
qDivPoly (QPoly [c]) (QPoly [d]) = QPoly [c / d]
qDivPoly p _                     = p  -- placeholder for general case

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
                            else qDivPoly (negQPow lcb delta)
                                          (qPow psi (max 0 (delta-1)))
              beta'     = qMulPoly (negQPoly lcb) (qPow psi' delta)
              r'        = scaleBiv (qInvPoly beta') r
          in a : go b r' beta' psi'
    negQPow q n  = qScalePoly ((-1)^n) (qPow q n)
    negQPoly (QPoly cs) = QPoly (map negate cs)
    qInvPoly (QPoly [c]) = QPoly [1/c]
    qInvPoly q           = q

-- | Extract constant term from BivPoly
extractConstant :: BivPoly -> QPoly
extractConstant (BivPoly [])    = QPoly []
extractConstant (BivPoly (c:_)) = c

-- | Compute the resultant of p(y) and q(x-y) over y
bivResultant :: QPoly -> QPoly -> QPoly
bivResultant pa pb =
  let paB = qPolyToBiv pa
      pbB = substituteXMinusY pb
      prs = subresultantPRS paB pbB
  in case prs of
       [] -> QPoly []
       _  -> extractConstant (last prs)

-- | Additive resultant: res_y(p(y), q(x-y))
addResultantQ :: QPoly -> QPoly -> QPoly
addResultantQ pa pb
  | qDegree pa == 1 && qDegree pb == 1 =
      let a = negate (head (qPolyCoef pa))
          b = negate (head (qPolyCoef pb))
      in QPoly [negate (a+b), 1]
  | otherwise = bivResultant pa pb

-- | Multiplicative resultant: res_y(p(y), y^n * q(x/y))
mulResultantQ :: QPoly -> QPoly -> QPoly
mulResultantQ pa pb
  | qDegree pa == 1 && qDegree pb == 1 =
      -- p(y) = y - a, q(y) = y - b
      -- minimal poly of a*b is y - a*b
      let a = negate (head (qPolyCoef pa))
          b = negate (head (qPolyCoef pb))
      in QPoly [negate (a*b), 1]
  | otherwise = bivResultant pa pb  -- general case (may loop for now)