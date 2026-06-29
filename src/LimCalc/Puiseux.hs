module LimCalc.Puiseux where

import Data.List (sortBy)
import Data.Ord (comparing)

-- | A term in a Puiseux series: coefficient * h^pExp
-- pExp is rational (allowing fractional powers like h^(1/2))
-- coeff is Double for now (placeholder for complex AlgNum)
data PuiseuxTerm = PuiseuxTerm
  { pExp  :: Rational  -- ^ The exponent (p/q for some integers p, q)
  , coeff :: Double    -- ^ The coefficient (placeholder for AlgNum)
  } deriving (Show, Eq)

-- | A Puiseux series is a list of terms sorted by ascending exponent.
newtype PuiseuxSeries = PuiseuxSeries
  { terms :: [PuiseuxTerm]
  } deriving (Show, Eq)

-- | Sort terms by exponent
normalize :: PuiseuxSeries -> PuiseuxSeries
normalize (PuiseuxSeries ts) =
  PuiseuxSeries $ sortBy (comparing pExp) ts

-- | The leading term of a series (lowest exponent)
leadingTerm :: PuiseuxSeries -> Maybe PuiseuxTerm
leadingTerm (PuiseuxSeries [])    = Nothing
leadingTerm (PuiseuxSeries (t:_)) = Just t

-- | Add two Puiseux series
addSeries :: PuiseuxSeries -> PuiseuxSeries -> PuiseuxSeries
addSeries (PuiseuxSeries s1) (PuiseuxSeries s2) =
  normalize $ PuiseuxSeries $ mergePlus s1 s2
  where
    mergePlus [] ys = ys
    mergePlus xs [] = xs
    mergePlus (x:xs) (y:ys)
      | pExp x == pExp y =
          let c = coeff x + coeff y
          in if c == 0
             then mergePlus xs ys
             else PuiseuxTerm (pExp x) c : mergePlus xs ys
      | pExp x < pExp y  = x : mergePlus xs (y:ys)
      | otherwise        = y : mergePlus (x:xs) ys

-- | Scale a series by a constant
scaleSeries :: Double -> PuiseuxSeries -> PuiseuxSeries
scaleSeries c (PuiseuxSeries ts) =
  PuiseuxSeries $ map (\t -> t { coeff = c * coeff t }) ts

-- | Multiply two Puiseux series (Cauchy product)
mulSeries :: PuiseuxSeries -> PuiseuxSeries -> PuiseuxSeries
mulSeries (PuiseuxSeries s1) (PuiseuxSeries s2) =
  normalize $ PuiseuxSeries
    [ PuiseuxTerm (pExp t1 + pExp t2) (coeff t1 * coeff t2)
    | t1 <- s1, t2 <- s2
    ]