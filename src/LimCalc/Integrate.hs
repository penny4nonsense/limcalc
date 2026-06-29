module LimCalc.Integrate where

import LimCalc.Expr
import LimCalc.Types
import LimCalc.Expand
import LimCalc.Puiseux
import LimCalc.SymExpand
import LimCalc.SymPuiseux
import LimCalc.Simplify
import qualified Data.Map.Strict as Map

-- | Symbolically integrate f with respect to var.
-- Integration is term-by-term on the Puiseux series:
-- h^n → h^(n+1)/(n+1)  for n ≠ -1
-- h^(-1) → log(h)
-- Returns an Expr representing the antiderivative (no constant of integration).
integrate :: Expr -> String -> Either ExpandError Expr
integrate f var = do
  series <- symExpand f var
  return $ simplify $ sumTerms (integrateSymSeries series)

-- | Integrate a symbolic Puiseux series term by term
integrateSymSeries :: SymPuiseuxSeries -> SymPuiseuxSeries
integrateSymSeries (SymPuiseuxSeries ts) =
  SymPuiseuxSeries $ map integrateTerm ts

-- | Integrate a single term: c * h^n → c * h^(n+1)/(n+1)
-- Special case: c * h^(-1) → c * log(h)
integrateTerm :: SymPuiseuxTerm -> SymPuiseuxTerm
integrateTerm (SymPuiseuxTerm n c)
  | n == -1   = SymPuiseuxTerm 0 (Mul c (Log (Var "__h__")))  -- log term
  | otherwise = SymPuiseuxTerm (n + 1) (Div c (Const (fromRational (n + 1))))

-- | Sum all terms of a symbolic series into a single Expr
sumTerms :: SymPuiseuxSeries -> Expr
sumTerms (SymPuiseuxSeries [])     = Const 0
sumTerms (SymPuiseuxSeries (t:ts)) =
  foldl (\acc t' -> Add acc (termToExpr t')) (termToExpr t) ts

-- | Convert a symbolic term back to an Expr
-- The series variable is the expansion variable
termToExpr :: SymPuiseuxTerm -> Expr
termToExpr (SymPuiseuxTerm n c)
  | n == 0    = c
  | n == 1    = Mul c (Var "__h__")
  | otherwise = Mul c (Pow (Var "__h__") (Const (fromRational n)))

-- | Definite integral of f from a to b with respect to var
-- Uses the antiderivative evaluated at b minus a
definiteIntegral :: Expr -> String -> Double -> Double -> Either ExpandError Double
definiteIntegral f var a b = do
  series <- expand f (Map.fromList [(var, a)]) var
  return $ integrateNumeric series (b - a)

-- | Numerically integrate a Puiseux series from 0 to h
-- by summing term by term
integrateNumeric :: PuiseuxSeries -> Double -> Double
integrateNumeric (PuiseuxSeries ts) h =
  sum [ integrateNumericTerm t h | t <- ts ]

-- | Integrate a single numeric term from 0 to h
integrateNumericTerm :: PuiseuxTerm -> Double -> Double
integrateNumericTerm (PuiseuxTerm n c) h
  | n == -1   = c * log h
  | otherwise = c * h ** (fromRational n + 1) / (fromRational n + 1)