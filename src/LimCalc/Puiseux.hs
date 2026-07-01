module LimCalc.Puiseux where

import Data.List (sortBy)
import Data.Ord  (comparing)

-- | A term in a log-Puiseux series: coeff * h^lpExp * log(h)^lpLog
data LogPuiseuxTerm a = LogPuiseuxTerm
  { lpCoeff :: a
  , lpExp   :: Rational  -- ^ exponent of h
  , lpLog   :: Int       -- ^ exponent of log(h); 0 = pure power term
  } deriving (Show, Eq)

-- | A log-Puiseux series: finite list of terms in canonical order
newtype LogPuiseuxSeries a = LogPuiseuxSeries
  { lterms :: [LogPuiseuxTerm a]
  } deriving (Show, Eq)

-- | Smart constructor: pure power term (no log factor)
pureTerm :: Rational -> a -> LogPuiseuxTerm a
pureTerm e c = LogPuiseuxTerm c e 0

-- | Smart constructor: term with log factor
logTerm :: Rational -> Int -> a -> LogPuiseuxTerm a
logTerm e k c = LogPuiseuxTerm c e k

-- | Canonical ordering: ascending (lpExp, lpLog)
termOrd :: LogPuiseuxTerm a -> LogPuiseuxTerm a -> Ordering
termOrd = comparing (\t -> (lpExp t, lpLog t))

-- | Sort into canonical order
normalize :: LogPuiseuxSeries a -> LogPuiseuxSeries a
normalize (LogPuiseuxSeries ts) = LogPuiseuxSeries (sortBy termOrd ts)

-- | Leading term (lowest lpExp, then lowest lpLog within that)
leadingTerm :: LogPuiseuxSeries a -> Maybe (LogPuiseuxTerm a)
leadingTerm (LogPuiseuxSeries [])    = Nothing
leadingTerm (LogPuiseuxSeries (t:_)) = Just t

-- | Leading term with nonzero coefficient
leadingTermNZ :: (Floating a, Ord a) => LogPuiseuxSeries a -> Maybe (LogPuiseuxTerm a)
leadingTermNZ (LogPuiseuxSeries ts) =
  case filter (\t -> abs (lpCoeff t) > 1e-12) ts of
    []    -> Nothing
    (t:_) -> Just t

-- | Remove terms with near-zero coefficients
stripZeros :: (Floating a, Ord a) => LogPuiseuxSeries a -> LogPuiseuxSeries a
stripZeros (LogPuiseuxSeries ts) =
  LogPuiseuxSeries $ filter (\t -> abs (lpCoeff t) > 1e-12) ts

-- | Combine like terms: same (lpExp, lpLog) pair
combineLike :: (Floating a, Ord a)
            => LogPuiseuxSeries a -> LogPuiseuxSeries a -> LogPuiseuxSeries a
combineLike (LogPuiseuxSeries s1) (LogPuiseuxSeries s2) =
  stripZeros $ normalize $ LogPuiseuxSeries $ mergePlus s1 s2
  where
    key t = (lpExp t, lpLog t)
    mergePlus [] ys = ys
    mergePlus xs [] = xs
    mergePlus (x:xs) (y:ys) = case compare (key x) (key y) of
      EQ -> LogPuiseuxTerm (lpCoeff x + lpCoeff y) (lpExp x) (lpLog x)
              : mergePlus xs ys
      LT -> x : mergePlus xs (y:ys)
      GT -> y : mergePlus (x:xs) ys

-- | Add two log-Puiseux series
addSeries :: (Floating a, Ord a)
          => LogPuiseuxSeries a -> LogPuiseuxSeries a -> LogPuiseuxSeries a
addSeries = combineLike

-- | Scale by a constant
scaleSeries :: (Floating a, Ord a)
            => a -> LogPuiseuxSeries a -> LogPuiseuxSeries a
scaleSeries c (LogPuiseuxSeries ts) =
  stripZeros $ LogPuiseuxSeries
    [ LogPuiseuxTerm (c * lpCoeff t) (lpExp t) (lpLog t) | t <- ts ]

-- | Multiply two log-Puiseux series.
-- Product rule: (c1 h^p1 log^k1) * (c2 h^p2 log^k2)
--             = c1*c2 * h^(p1+p2) * log^(k1+k2)
-- Type is closed under multiplication.
mulSeries :: (Floating a, Ord a)
          => LogPuiseuxSeries a -> LogPuiseuxSeries a -> LogPuiseuxSeries a
mulSeries (LogPuiseuxSeries s1) (LogPuiseuxSeries s2) =
  stripZeros $ normalize $
    foldr combineLike (LogPuiseuxSeries [])
      [ LogPuiseuxSeries
          [ LogPuiseuxTerm (lpCoeff t1 * lpCoeff t2)
                           (lpExp t1 + lpExp t2)
                           (lpLog t1 + lpLog t2) ]
      | t1 <- s1, t2 <- s2
      ]

-- | Constant term: coefficient of h^0 * log(h)^0
constantTerm :: (Floating a, Ord a) => LogPuiseuxSeries a -> a
constantTerm (LogPuiseuxSeries []) = 0
constantTerm (LogPuiseuxSeries (t:_))
  | lpExp t == 0 && lpLog t == 0 = lpCoeff t
  | otherwise                    = 0

-- | Extract the pure h^0 coefficient regardless of log power
-- (used for log-expansion where the constant may carry a log term)
constantCoeff :: (Floating a, Ord a) => Rational -> LogPuiseuxSeries a -> [(Int, a)]
constantCoeff e (LogPuiseuxSeries ts) =
  [ (lpLog t, lpCoeff t) | t <- ts, lpExp t == e ]

-- | Remove all terms at a given (lpExp, lpLog) pair
removeTerm :: Rational -> Int -> LogPuiseuxSeries a -> LogPuiseuxSeries a
removeTerm e k (LogPuiseuxSeries ts) =
  LogPuiseuxSeries $ filter (\t -> not (lpExp t == e && lpLog t == k)) ts

-- | Remove all terms at a given lpExp (all log powers)
removeExp :: Rational -> LogPuiseuxSeries a -> LogPuiseuxSeries a
removeExp e (LogPuiseuxSeries ts) =
  LogPuiseuxSeries $ filter (\t -> lpExp t /= e) ts

-- | Truncate to terms at or below a given order (lpExp)
truncateToOrder :: Rational -> LogPuiseuxSeries a -> LogPuiseuxSeries a
truncateToOrder maxOrder (LogPuiseuxSeries ts) =
  LogPuiseuxSeries $ filter (\t -> lpExp t <= maxOrder) ts

-- | Shift all lpExp values by delta
shiftExponents :: Rational -> LogPuiseuxSeries a -> LogPuiseuxSeries a
shiftExponents delta (LogPuiseuxSeries ts) =
  LogPuiseuxSeries [ LogPuiseuxTerm (lpCoeff t) (lpExp t + delta) (lpLog t)
                   | t <- ts ]

-- | True if the series has any log terms (lpLog > 0)
hasLogTerms :: LogPuiseuxSeries a -> Bool
hasLogTerms (LogPuiseuxSeries ts) = any (\t -> lpLog t > 0) ts