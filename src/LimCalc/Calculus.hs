module LimCalc.Calculus where

import LimCalc.Expr
import LimCalc.Puiseux
import LimCalc.Types
import LimCalc.Expand
import LimCalc.SymExpand
import LimCalc.SymPuiseux

-- | Symbolic derivative of f with respect to var.
-- Returns an Expr representing f'(x) — the derivative as a function.
diff :: Expr -> String -> Either ExpandError Expr
diff f var = do
  series <- symExpand f var
  return $ symCoeffAt 1 series

-- | The derivative of f with respect to var at point.
-- Defined as the coefficient of h^1 in the Puiseux expansion of f(point + h*var).
-- This is lim_{h→0} (f(x+h) - f(x))/h emerging analytically from the series.
derivative :: Expr -> Point -> String -> Either ExpandError Double
derivative f point var = do
  series <- expand f point var
  return $ coeffAt 1 series

-- | Extract the coefficient of h^n from a series
coeffAt :: Rational -> PuiseuxSeries -> Double
coeffAt n (PuiseuxSeries ts) =
  case filter (\t -> pExp t == n) ts of
    []    -> 0
    (t:_) -> coeff t

-- | The nth derivative of f at point with respect to var.
-- The nth derivative is n! times the coefficient of h^n in the expansion.
nthDerivative :: Int -> Expr -> Point -> String -> Either ExpandError Double
nthDerivative n f point var = do
  series <- expand f point var
  return $ factorial n * coeffAt (fromIntegral n) series

-- | Factorial
factorial :: Int -> Double
factorial 0 = 1
factorial n = fromIntegral n * factorial (n-1)

-- | Partial derivative of f with respect to var at point.
-- Same as derivative — the Point handles the multivariate case.
partial :: Expr -> Point -> String -> Either ExpandError Double
partial = derivative