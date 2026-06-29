module LimCalc.SymPuiseux where

import Data.List (sortBy)
import Data.Ord (comparing)
import LimCalc.Expr

-- | A term in a symbolic Puiseux series: symCoeff * h^symExp
-- symCoeff is an Expr — coefficients are symbolic expressions
data SymPuiseuxTerm = SymPuiseuxTerm
  { symExp   :: Rational  -- ^ The exponent
  , symCoeff :: Expr      -- ^ The coefficient as an expression
  } deriving (Show, Eq)

-- | A symbolic Puiseux series
newtype SymPuiseuxSeries = SymPuiseuxSeries
  { symTerms :: [SymPuiseuxTerm]
  } deriving (Show, Eq)

-- | Zero symbolic series
zeroSym :: SymPuiseuxSeries
zeroSym = SymPuiseuxSeries []

-- | Sort terms by exponent
normalizeSym :: SymPuiseuxSeries -> SymPuiseuxSeries
normalizeSym (SymPuiseuxSeries ts) =
  SymPuiseuxSeries $ sortBy (comparing symExp) ts

-- | Combine terms with equal exponents by adding coefficients
combineLikeSym :: SymPuiseuxSeries -> SymPuiseuxSeries -> SymPuiseuxSeries
combineLikeSym (SymPuiseuxSeries s1) (SymPuiseuxSeries s2) =
  normalizeSym $ SymPuiseuxSeries $ mergePlus s1 s2
  where
    mergePlus [] ys = ys
    mergePlus xs [] = xs
    mergePlus (x:xs) (y:ys)
      | symExp x == symExp y =
          SymPuiseuxTerm (symExp x) (Add (symCoeff x) (symCoeff y))
            : mergePlus xs ys
      | symExp x < symExp y  = x : mergePlus xs (y:ys)
      | otherwise            = y : mergePlus (x:xs) ys

-- | Add two symbolic series
addSymSeries :: SymPuiseuxSeries -> SymPuiseuxSeries -> SymPuiseuxSeries
addSymSeries s1 s2 = combineLikeSym s1 s2

-- | Scale a symbolic series by an Expr
scaleSymSeries :: Expr -> SymPuiseuxSeries -> SymPuiseuxSeries
scaleSymSeries c (SymPuiseuxSeries ts) =
  SymPuiseuxSeries $ map (\t -> t { symCoeff = Mul c (symCoeff t) }) ts

-- | Multiply two symbolic series (Cauchy product)
mulSymSeries :: SymPuiseuxSeries -> SymPuiseuxSeries -> SymPuiseuxSeries
mulSymSeries (SymPuiseuxSeries s1) (SymPuiseuxSeries s2) =
  let raw = [ SymPuiseuxTerm (symExp t1 + symExp t2) (Mul (symCoeff t1) (symCoeff t2))
            | t1 <- s1, t2 <- s2 ]
      sorted = normalizeSym (SymPuiseuxSeries raw)
  in foldTerms sorted
  where
    foldTerms (SymPuiseuxSeries []) = SymPuiseuxSeries []
    foldTerms (SymPuiseuxSeries (t:ts)) =
      let (same, rest) = span (\t' -> symExp t' == symExp t) ts
          combined = foldl (\acc t' -> Add acc (symCoeff t')) (symCoeff t) same
          SymPuiseuxSeries tail = foldTerms (SymPuiseuxSeries rest)
      in SymPuiseuxSeries (SymPuiseuxTerm (symExp t) combined : tail)

-- | Get the coefficient of h^n — sum all terms with that exponent
symCoeffAt :: Rational -> SymPuiseuxSeries -> Expr
symCoeffAt n (SymPuiseuxSeries ts) =
  case filter (\t -> symExp t == n) ts of
    []     -> Const 0
    [t]    -> symCoeff t
    (t:ts) -> foldl (\acc t' -> Add acc (symCoeff t')) (symCoeff t) ts

-- | Leading term
symLeadingTerm :: SymPuiseuxSeries -> Maybe SymPuiseuxTerm
symLeadingTerm (SymPuiseuxSeries [])    = Nothing
symLeadingTerm (SymPuiseuxSeries (t:_)) = Just t

-- | Truncate to n terms
truncateSymSeries :: Int -> SymPuiseuxSeries -> SymPuiseuxSeries
truncateSymSeries n (SymPuiseuxSeries ts) =
  SymPuiseuxSeries (take n ts)

-- | Shift all exponents by delta
shiftSymExponents :: Rational -> SymPuiseuxSeries -> SymPuiseuxSeries
shiftSymExponents delta (SymPuiseuxSeries ts) =
  SymPuiseuxSeries [ SymPuiseuxTerm (symExp t + delta) (symCoeff t) | t <- ts ]

-- | Remove term with given exponent
removeSymTerm :: Rational -> SymPuiseuxSeries -> SymPuiseuxSeries
removeSymTerm e (SymPuiseuxSeries ts) =
  SymPuiseuxSeries $ filter (\t -> symExp t /= e) ts