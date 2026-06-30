module LimCalc.Puiseux where

import Data.List (sortBy)
import Data.Ord (comparing)

-- | A term in a Puiseux series: coefficient * h^pExp
data PuiseuxTerm a = PuiseuxTerm
  { pExp  :: Rational
  , coeff :: a
  } deriving (Show, Eq)

-- | A Puiseux series
newtype PuiseuxSeries a = PuiseuxSeries
  { terms :: [PuiseuxTerm a]
  } deriving (Show, Eq)

-- | Sort terms by exponent
normalize :: PuiseuxSeries a -> PuiseuxSeries a
normalize (PuiseuxSeries ts) =
  PuiseuxSeries $ sortBy (comparing pExp) ts

-- | Leading term
leadingTerm :: PuiseuxSeries a -> Maybe (PuiseuxTerm a)
leadingTerm (PuiseuxSeries [])    = Nothing
leadingTerm (PuiseuxSeries (t:_)) = Just t

-- | Leading term with nonzero coefficient
leadingTermNZ :: (Floating a, Ord a) => PuiseuxSeries a -> Maybe (PuiseuxTerm a)
leadingTermNZ (PuiseuxSeries ts) =
  case filter (\t -> abs (coeff t) > 1e-12) ts of
    []    -> Nothing
    (t:_) -> Just t

-- | Remove terms with near-zero coefficients
stripZeros :: (Floating a, Ord a) => PuiseuxSeries a -> PuiseuxSeries a
stripZeros (PuiseuxSeries ts) =
  PuiseuxSeries $ filter (\t -> abs (coeff t) > 1e-12) ts

-- | Combine terms with equal exponents
combineLike :: (Floating a, Ord a) => PuiseuxSeries a -> PuiseuxSeries a -> PuiseuxSeries a
combineLike (PuiseuxSeries s1) (PuiseuxSeries s2) =
  stripZeros $ normalize $ PuiseuxSeries $ mergePlus s1 s2
  where
    mergePlus [] ys = ys
    mergePlus xs [] = xs
    mergePlus (x:xs) (y:ys)
      | pExp x == pExp y =
          PuiseuxTerm (pExp x) (coeff x + coeff y)
            : mergePlus xs ys
      | pExp x < pExp y  = x : mergePlus xs (y:ys)
      | otherwise        = y : mergePlus (x:xs) ys

-- | Add two Puiseux series
addSeries :: (Floating a, Ord a) => PuiseuxSeries a -> PuiseuxSeries a -> PuiseuxSeries a
addSeries s1 s2 = combineLike s1 s2

-- | Scale a series by a constant
scaleSeries :: (Floating a, Ord a) => a -> PuiseuxSeries a -> PuiseuxSeries a
scaleSeries c (PuiseuxSeries ts) =
  stripZeros $ PuiseuxSeries $ map (\t -> t { coeff = c * coeff t }) ts

-- | Multiply two Puiseux series (Cauchy product)
--
-- combineLike's mergePlus assumes both of its inputs already have at
-- most one term per exponent (that invariant is what lets addSeries
-- just call it directly on two series). The raw Cauchy cross-product
-- below does NOT have that property -- multiple (t1, t2) pairs can
-- land on the same exponent -- so dumping it straight into a single
-- combineLike call (merged against an empty series) was a no-op:
-- mergePlus xs [] = xs, meaning no consolidation ever happened, and
-- duplicate-exponent terms survived side by side after sorting.
--
-- Folding each cross term through combineLike one at a time restores
-- the intended invariant at every step, so duplicate exponents are
-- actually summed.
mulSeries :: (Floating a, Ord a) => PuiseuxSeries a -> PuiseuxSeries a -> PuiseuxSeries a
mulSeries (PuiseuxSeries s1) (PuiseuxSeries s2) =
  stripZeros $ normalize $
    foldr combineLike (PuiseuxSeries [])
      [ PuiseuxSeries [PuiseuxTerm (pExp t1 + pExp t2) (coeff t1 * coeff t2)]
      | t1 <- s1, t2 <- s2
      ]

-- | Get the constant term (h^0 coefficient)
constantTerm :: (Floating a, Ord a) => PuiseuxSeries a -> a
constantTerm (PuiseuxSeries []) = 0
constantTerm (PuiseuxSeries (t:_))
  | pExp t == 0 = coeff t
  | otherwise   = 0