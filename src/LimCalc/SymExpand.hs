module LimCalc.SymExpand where

import LimCalc.Expr
import LimCalc.SymPuiseux
import LimCalc.Types
import LimCalc.DiffField (deriveBase)
import LimCalc.Simplify
import Data.Ratio (numerator, denominator)

-- | Symbolically expand f(x + h) as a Puiseux series in h.
-- var: the variable we are expanding in (x)
-- All other variables remain as Var nodes in the coefficients.
symExpand :: Expr -> String -> Either ExpandError SymPuiseuxSeries

-- Constants
symExpand (Const c) _ = Right $ SymPuiseuxSeries [SymPuiseuxTerm 0 (Const c)]
symExpand Pi        _ = Right $ SymPuiseuxSeries [SymPuiseuxTerm 0 Pi]
symExpand E         _ = Right $ SymPuiseuxSeries [SymPuiseuxTerm 0 E]
symExpand I         _ = Right $ SymPuiseuxSeries [SymPuiseuxTerm 0 I]

-- Variable
symExpand (Var name) var
  | name == var =
      -- expanding in this variable: x + h
      Right $ SymPuiseuxSeries
        [ SymPuiseuxTerm 0 (Var name)  -- constant term x
        , SymPuiseuxTerm 1 (Const 1)   -- linear term h
        ]
  | otherwise =
      -- not the expansion variable: stays symbolic
      Right $ SymPuiseuxSeries [SymPuiseuxTerm 0 (Var name)]

-- Addition
symExpand (Add f g) var = do
  sf <- symExpand f var
  sg <- symExpand g var
  return $ addSymSeries sf sg

-- Subtraction
symExpand (Sub f g) var = do
  sf <- symExpand f var
  sg <- symExpand g var
  return $ addSymSeries sf (scaleSymSeries (Const (-1)) sg)

-- Multiplication
symExpand (Mul f g) var = do
  sf <- symExpand f var
  sg <- symExpand g var
  return $ mulSymSeries sf sg

-- Negation
symExpand (Neg f) var = do
  sf <- symExpand f var
  return $ scaleSymSeries (Const (-1)) sf

-- Sin: sin(f(x+h)) via symbolic Taylor series
symExpand (Sin f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = sinSymTaylor c0
  return $ symEvalSeriesAt u s

-- Cos: cos(f(x+h)) via symbolic Taylor series
symExpand (Cos f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = cosSymTaylor c0
  return $ symEvalSeriesAt u s

-- Exp: e^f(x+h) via symbolic Taylor series
symExpand (Exp f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = expSymTaylor c0
  return $ symEvalSeriesAt u s

-- Log
symExpand (Log f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = logSymTaylor c0
  return $ symEvalSeriesAt u s

-- Division
symExpand (Div f g) var = do
  sf <- symExpand f var
  sg <- symExpand g var
  return $ mulSymSeries sf (symInvertSeries sg)

-- Power
symExpand (Pow f (Const r)) var = symExpandPowR f (toRational r) var
symExpand (Pow f (Neg (Const r))) var = symExpandPowR f (toRational (-r)) var
symExpand (Pow _ _) _ = Left $ Unknown "Symbolic exponents not yet supported"

-- Abs
symExpand (Abs _) _ = Left $ Unknown "Abs not yet implemented"

-- Erf: analytic everywhere, Taylor series via iterated derivatives
-- erf(c0 + h) = erf(c0) + erf'(c0)h + erf''(c0)h^2/2 + ...
-- where erf^(n)(x) = (2/sqrt(pi)) * H_{n-1}(x) * e^(-x^2)
-- (Hermite polynomial coefficients), computed via deriveBase.
symExpand (Erf f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = erfSymTaylor c0
  return $ symEvalSeriesAt u s

-- Si: analytic everywhere (Si(0)=0), Taylor series via iterated
-- derivatives. Si^(0)=Si(x), Si^(1)=sin(x)/x, Si^(n) via chain rule.
symExpand (Si f) var = do
  sf <- symExpand f var
  let c0 = symConstantTerm sf
      u  = removeSymTerm 0 sf
      s  = siSymTaylor c0
  return $ symEvalSeriesAt u s

-- Li, Ci, Ei: singular at the natural expansion point (li(x) and
-- Ci(x) have logarithmic singularities; Ei(x) is singular at 0).
-- Full Puiseux expansion would require logarithmic terms similar to
-- how Log is handled, but the Taylor coefficient structure is more
-- complex (involving li/Ci/Ei themselves). Not yet implemented.
symExpand (Li _)  _ = Left $ Unknown
  "Li (logarithmic integral) expansion not yet implemented: \
  \li(x) has a logarithmic singularity at x=0 and x=1"
symExpand (Ci _)  _ = Left $ Unknown
  "Ci (cosine integral) expansion not yet implemented: \
  \Ci(x) has a logarithmic singularity at x=0"
symExpand (Ei _)  _ = Left $ Unknown
  "Ei (exponential integral) expansion not yet implemented: \
  \Ei(x) has a logarithmic singularity at x=0"

-- | How many terms
depth :: Int
depth = 6

-- | Get the constant term of a symbolic series
symConstantTerm :: SymPuiseuxSeries -> Expr
symConstantTerm (SymPuiseuxSeries []) = Const 0
symConstantTerm (SymPuiseuxSeries (t:_))
  | symExp t == 0 = symCoeff t
  | otherwise     = Const 0

-- | Evaluate S = Σ aₙ·t^n by substituting t = u (symbolic)
symEvalSeriesAt :: SymPuiseuxSeries -> SymPuiseuxSeries -> SymPuiseuxSeries
symEvalSeriesAt u (SymPuiseuxSeries ts) =
  foldr addSymSeries zeroSym
    [ scaleSymSeries (symCoeff t) (symPowSeries u (symExp t))
    | t <- ts
    ]

-- | Raise a symbolic series to an integer power
symPowSeries :: SymPuiseuxSeries -> Rational -> SymPuiseuxSeries
symPowSeries _ 0 = SymPuiseuxSeries [SymPuiseuxTerm 0 (Const 1)]
symPowSeries u n
  | n > 0 && denominator n == 1 =
      let k = fromIntegral (numerator n) :: Int
      in foldl mulSymSeries (SymPuiseuxSeries [SymPuiseuxTerm 0 (Const 1)]) (replicate k u)
  | otherwise = SymPuiseuxSeries []

-- | Symbolic Taylor series for sin around symbolic x
-- sin(x + h) = sin(x) + cos(x)h - sin(x)h²/2 - cos(x)h³/6 + ...
sinSymTaylor :: Expr -> SymPuiseuxSeries
sinSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (sinSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    sinSymCoeff n =
      let bases  = cycle [Sin x, Cos x, Neg (Sin x), Neg (Cos x)]
          facts  = scanl (*) 1 [1..] :: [Int]
          base   = bases !! n
          factor = fromIntegral (facts !! n)
      in if factor == 1
         then base
         else Div base (Const factor)

-- | Symbolic Taylor series for cos around symbolic x
cosSymTaylor :: Expr -> SymPuiseuxSeries
cosSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (cosSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    cosSymCoeff n =
      let bases  = cycle [Cos x, Neg (Sin x), Neg (Cos x), Sin x]
          facts  = scanl (*) 1 [1..] :: [Int]
          base   = bases !! n
          factor = fromIntegral (facts !! n)
      in if factor == 1
         then base
         else Div base (Const factor)

-- | Symbolic Taylor series for exp around symbolic x
-- exp(x + h) = exp(x) + exp(x)h + exp(x)h²/2 + ...
expSymTaylor :: Expr -> SymPuiseuxSeries
expSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (expSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    expSymCoeff n =
      let facts  = scanl (*) 1 [1..] :: [Int]
          factor = fromIntegral (facts !! n)
      in if factor == 1
         then Exp x
         else Div (Exp x) (Const factor)

-- | Symbolic Taylor series for log around symbolic x
-- log(x + h) = log(x) + h/x - h²/(2x²) + ...
logSymTaylor :: Expr -> SymPuiseuxSeries
logSymTaylor x = SymPuiseuxSeries $ take depth $
  SymPuiseuxTerm 0 (Log x) :
  [ SymPuiseuxTerm (fromIntegral n) (logSymCoeff n)
  | n <- [1..] :: [Int] ]
  where
    logSymCoeff n =
      let sign   = if even n then Const (-1) else Const 1
          factor = Const (fromIntegral n)
          xpow   = Pow x (Const (fromIntegral n))
      in Div sign (Mul factor xpow)

-- | Symbolic Taylor series for erf around symbolic x.
-- erf(x + h) = erf(x) + erf'(x)h + erf''(x)h^2/2! + ...
-- The n-th derivative is computed by iterating deriveBase.
-- Coefficients grow in size for large n due to the Hermite polynomial
-- structure, but remain correct and work well for the primary use
-- case of diff (Erf f) var (which only needs the n=1 coefficient).
erfSymTaylor :: Expr -> SymPuiseuxSeries
erfSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (erfSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    -- Memoized list of nth derivatives: [erf(x), erf'(x), erf''(x), ...]
    derivs = iterate (\e -> simplify (deriveBase e)) (Erf x)
    erfSymCoeff 0 = Erf x
    erfSymCoeff n =
      let d      = derivs !! n
          facts  = scanl (*) 1 [1..] :: [Int]
          factor = fromIntegral (facts !! n)
      in if factor == 1 then d else Div d (Const factor)

-- | Symbolic Taylor series for Si around symbolic x.
-- Si(x + h) = Si(x) + Si'(x)h + Si''(x)h^2/2! + ...
-- Si'(x) = sin(x)/x, subsequent derivatives via iterating deriveBase.
siSymTaylor :: Expr -> SymPuiseuxSeries
siSymTaylor x = SymPuiseuxSeries $ take depth
  [ SymPuiseuxTerm (fromIntegral n) (siSymCoeff n)
  | n <- [0..] :: [Int] ]
  where
    derivs = iterate (\e -> simplify (deriveBase e)) (Si x)
    siSymCoeff 0 = Si x
    siSymCoeff n =
      let d      = derivs !! n
          facts  = scanl (*) 1 [1..] :: [Int]
          factor = fromIntegral (facts !! n)
      in if factor == 1 then d else Div d (Const factor)
symInvertSeries :: SymPuiseuxSeries -> SymPuiseuxSeries
symInvertSeries s =
  case symLeadingTerm s of
    Nothing -> SymPuiseuxSeries []
    Just lt ->
      let alpha = symExp lt
          a     = symCoeff lt
          w     = symNormalizeW s lt alpha a
          negw  = scaleSymSeries (Const (-1)) w
          geo   = symGeometricSeries negw
          shift = negate alpha
      in shiftSymExponents shift
           (scaleSymSeries (Div (Const 1) a) geo)

-- | Geometric series 1/(1+w) = Σ (-w)^n
symGeometricSeries :: SymPuiseuxSeries -> SymPuiseuxSeries
symGeometricSeries u =
  let upows = take (depth+1) $ iterate (truncateSymSeries depth . mulSymSeries u)
                                       (SymPuiseuxSeries [SymPuiseuxTerm 0 (Const 1)])
  in truncateSymSeries depth $ foldr addSymSeries zeroSym upows

-- | Normalize w = s/(a*h^alpha) - 1
symNormalizeW :: SymPuiseuxSeries -> SymPuiseuxTerm -> Rational -> Expr -> SymPuiseuxSeries
symNormalizeW (SymPuiseuxSeries ts) _lt alpha a =
  let shifted = [ SymPuiseuxTerm (symExp t - alpha) (Div (symCoeff t) a) | t <- ts ]
  in removeSymTerm 0 (SymPuiseuxSeries shifted)

-- | Symbolic power expansion f^r
symExpandPowR :: Expr -> Rational -> String -> Either ExpandError SymPuiseuxSeries
symExpandPowR f r var = do
  s <- symExpand f var
  let s' = truncateSymSeries depth s
  case symLeadingTerm s' of
    Nothing -> Right zeroSym
    Just lt ->
      let alpha = symExp lt
          a     = symCoeff lt
          w     = truncateSymSeries depth (symNormalizeW s' lt alpha a)
          binom = symBinomialSeries r w
          shift = alpha * r
      in Right $ truncateSymSeries depth
           (shiftSymExponents shift
             (scaleSymSeries (Pow a (Const (fromRational r))) binom))

-- | Symbolic binomial series (1+w)^r
symBinomialSeries :: Rational -> SymPuiseuxSeries -> SymPuiseuxSeries
symBinomialSeries r w =
  let bcs   = symBinomCoeffs r depth
      wpows = take (depth+1) $ iterate (truncateSymSeries depth . mulSymSeries w)
                                       (SymPuiseuxSeries [SymPuiseuxTerm 0 (Const 1)])
  in truncateSymSeries depth $ foldr addSymSeries zeroSym
       [ scaleSymSeries c wp
       | (c, wp) <- zip bcs wpows
       ]

-- | Symbolic binomial coefficients as Expr
symBinomCoeffs :: Rational -> Int -> [Expr]
symBinomCoeffs r n = take (n+1) $ scanl step (Const 1) [0..n-1]
  where
    step acc k =
      let rk = fromRational r - fromIntegral k
          kk = fromIntegral (k+1)
      in Div (Mul acc (Const rk)) (Const kk)