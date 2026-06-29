module LimCalc.Expand where

import LimCalc.Expr
import LimCalc.Puiseux

-- | Expand f(x₀ + h) as a Puiseux series in h.
-- This is the core engine — all calculus operations derive from this.
expand :: Expr -> Double -> PuiseuxSeries

-- Constants are trivially constant series
expand (Const c) _ = PuiseuxSeries [PuiseuxTerm 0 c]

-- Pi and E as numeric constants for now
expand Pi _ = PuiseuxSeries [PuiseuxTerm 0 pi]
expand E  _ = PuiseuxSeries [PuiseuxTerm 0 (exp 1)]

-- Imaginary unit — placeholder, will need AlgNum
expand I  _ = PuiseuxSeries [PuiseuxTerm 0 0]  -- TODO: complex

-- Var x at x₀: f(x₀ + h) = x₀ + h
-- So the series is x₀*h^0 + 1*h^1
expand (Var _) x0 = PuiseuxSeries
  [ PuiseuxTerm 0 x0   -- constant term x₀
  , PuiseuxTerm 1 1.0  -- linear term h
  ]

-- Addition: expand both and add series
expand (Add f g) x0 = addSeries (expand f x0) (expand g x0)

-- Subtraction
expand (Sub f g) x0 = addSeries (expand f x0) (scaleSeries (-1) (expand g x0))

-- Multiplication: expand both and multiply series
expand (Mul f g) x0 = mulSeries (expand f x0) (expand g x0)

-- Negation
expand (Neg f) x0 = scaleSeries (-1) (expand f x0)

-- Division: f/g = f * g^(-1)
-- For now: expand g, invert series, multiply
expand (Div f g) x0 = mulSeries (expand f x0) (invertSeries (expand g x0))

-- Power: f^g -- the hard case, handled separately
expand (Pow f g) x0 = expandPow f g x0

-- Exp: e^f(x₀+h) -- compose exp series with expansion of f
expand (Exp f) x0 = composeSeries expSeries (expand f x0) x0

-- Log: ln(f(x₀+h))
expand (Log f) x0 = composeSeries logSeries (expand f x0) x0

-- Sin: sin(f(x₀+h))
expand (Sin f) x0 = composeSeries sinSeries (expand f x0) x0

-- Cos: cos(f(x₀+h))
expand (Cos f) x0 = composeSeries cosSeries (expand f x0) x0

-- Abs: non-analytic, handle separately
expand (Abs f) x0 = expandAbs f x0

-- | How many terms to compute in series expansions
depth :: Int
depth = 8

-- | Series for sin(u) around u=0: u - u³/6 + u⁵/120 - ...
sinSeries :: PuiseuxSeries
sinSeries = PuiseuxSeries
  [ PuiseuxTerm 1 1.0
  , PuiseuxTerm 3 (-1/6)
  , PuiseuxTerm 5 (1/120)
  , PuiseuxTerm 7 (-1/5040)
  ]

-- | Series for cos(u) around u=0: 1 - u²/2 + u⁴/24 - ...
cosSeries :: PuiseuxSeries
cosSeries = PuiseuxSeries
  [ PuiseuxTerm 0 1.0
  , PuiseuxTerm 2 (-1/2)
  , PuiseuxTerm 4 (1/24)
  , PuiseuxTerm 6 (-1/720)
  ]

-- | Series for exp(u) around u=0: 1 + u + u²/2 + u³/6 + ...
expSeries :: PuiseuxSeries
expSeries = PuiseuxSeries
  [ PuiseuxTerm 0 1.0
  , PuiseuxTerm 1 1.0
  , PuiseuxTerm 2 (1/2)
  , PuiseuxTerm 3 (1/6)
  , PuiseuxTerm 4 (1/24)
  ]

-- | Series for log(1+u) around u=0: u - u²/2 + u³/3 - ...
logSeries :: PuiseuxSeries
logSeries = PuiseuxSeries
  [ PuiseuxTerm 1 1.0
  , PuiseuxTerm 2 (-1/2)
  , PuiseuxTerm 3 (1/3)
  , PuiseuxTerm 4 (-1/4)
  , PuiseuxTerm 5 (1/5)
  ]

-- | Compose a known series S with an expansion E.
-- Computes S(E(h)) by substituting E into S term by term.
-- E must have zero constant term (i.e. E(0) = 0) for this to work directly.
-- For the general case we shift by the constant term first.
composeSeries :: PuiseuxSeries -> PuiseuxSeries -> Double -> PuiseuxSeries
composeSeries s e x0 =
  let c0    = constantTerm e          -- constant term of e
      e'    = shiftSeries (-c0) e     -- e shifted so e'(0) = 0
      s'    = shiftInput c0 s         -- s evaluated at c0 + e'
  in truncateSeries depth s'
  where
    shiftInput c0 (PuiseuxSeries ts) =
      -- substitute (c0 + e') into s by evaluating s at c0 first
      -- then adding correction terms -- stub for now
      PuiseuxSeries ts  -- TODO: full composition

-- | Get the constant term (h^0 coefficient) of a series
constantTerm :: PuiseuxSeries -> Double
constantTerm (PuiseuxSeries []) = 0
constantTerm (PuiseuxSeries (t:ts))
  | pExp t == 0 = coeff t
  | otherwise   = 0

-- | Shift a series by adding a constant to its constant term
shiftSeries :: Double -> PuiseuxSeries -> PuiseuxSeries
shiftSeries c s = addSeries s (PuiseuxSeries [PuiseuxTerm 0 c])

-- | Truncate a series to n terms
truncateSeries :: Int -> PuiseuxSeries -> PuiseuxSeries
truncateSeries n (PuiseuxSeries ts) = PuiseuxSeries (take n ts)

-- | Invert a series: compute 1/s
-- Stub for now
invertSeries :: PuiseuxSeries -> PuiseuxSeries
invertSeries s = s  -- TODO

-- | Handle Pow case: f^g
-- Stub for now
expandPow :: Expr -> Expr -> Double -> PuiseuxSeries
expandPow f g x0 = PuiseuxSeries []  -- TODO

-- | Handle Abs case
expandAbs :: Expr -> Double -> PuiseuxSeries
expandAbs f x0 = PuiseuxSeries []  -- TODO