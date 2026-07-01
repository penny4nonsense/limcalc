-- | Log-Puiseux series arithmetic.
--
-- A /log-Puiseux series/ is a finite formal sum of terms of the form
--
-- @c · h^p · log(h)^k@
--
-- where @c@ is a coefficient, @p ∈ ℚ@ is the /power exponent/, and
-- @k ∈ ℕ@ is the /log exponent/. Pure Puiseux series (no log factors)
-- are the special case @k = 0@ throughout.
--
-- This type arose as the natural representation for expansions of
-- functions with logarithmic singularities: @log(h)@, @Ci(h)@, and
-- @Ei(h)@ all have genuine @log(h)@ terms near @h = 0@ that cannot
-- be represented as pure power series. The log-Puiseux type is closed
-- under addition and multiplication (since @log^k · log^j = log^(k+j)@),
-- making it a suitable foundation for the expansion engine.
--
-- /Canonical form/: terms are stored in ascending lexicographic order
-- on @(lpExp, lpLog)@, with near-zero coefficients (below @1e-12@
-- in absolute value) stripped. All arithmetic operations maintain
-- this invariant.
--
-- /Li and iterated logs/: the logarithmic integral @li(x)@ has a
-- @log(log(h))@ singularity at @x = 0@, which is outside this type.
-- That case is documented as a known gap; the type is otherwise
-- sufficient for all special functions currently in 'Expr'.
module LimCalc.Puiseux
  ( -- * Term type
    LogPuiseuxTerm (..)
    -- * Series type
  , LogPuiseuxSeries (..)
    -- * Smart constructors
  , pureTerm
  , logTerm
    -- * Normalisation
  , normalize
  , stripZeros
    -- * Accessors
  , leadingTerm
  , leadingTermNZ
  , constantTerm
  , constantCoeff
  , hasLogTerms
    -- * Arithmetic
  , addSeries
  , scaleSeries
  , mulSeries
    -- * Structural operations
  , removeTerm
  , removeExp
  , truncateToOrder
  , shiftExponents
  ) where

import Data.List (sortBy)
import Data.Ord  (comparing)

-- | A single term in a log-Puiseux series: @c · h^p · log(h)^k@.
--
-- The fields are:
--
-- * 'lpCoeff' — the coefficient @c@
-- * 'lpExp'   — the rational power @p@ of @h@
-- * 'lpLog'   — the non-negative integer power @k@ of @log(h)@;
--               @k = 0@ gives a pure power term with no log factor
data LogPuiseuxTerm a = LogPuiseuxTerm
  { lpCoeff :: a
    -- ^ Coefficient of the term.
  , lpExp   :: Rational
    -- ^ Exponent of @h@. May be negative (pole), zero (constant),
    -- positive integer (regular power), or positive rational
    -- (branch point / fractional power).
  , lpLog   :: Int
    -- ^ Exponent of @log(h)@. Zero for pure power terms.
    -- Positive values arise from expansions of @log@, @Ci@, @Ei@,
    -- and similar functions near singular points.
  } deriving (Show, Eq)

-- | A log-Puiseux series: a finite list of terms in canonical order.
--
-- Canonical order is ascending lexicographic on @('lpExp', 'lpLog')@,
-- with no two terms sharing the same @(lpExp, lpLog)@ pair and no
-- term with a near-zero coefficient. All arithmetic operations
-- maintain this invariant via 'combineLike' and 'stripZeros'.
newtype LogPuiseuxSeries a = LogPuiseuxSeries
  { lterms :: [LogPuiseuxTerm a]
    -- ^ The list of terms in canonical order.
  } deriving (Show, Eq)

-- | Construct a pure power term @c · h^p@ (no log factor).
pureTerm :: Rational -> a -> LogPuiseuxTerm a
pureTerm e c = LogPuiseuxTerm c e 0

-- | Construct a log term @c · h^p · log(h)^k@.
logTerm :: Rational
        -> Int       -- ^ Log exponent @k@ (must be non-negative).
        -> a
        -> LogPuiseuxTerm a
logTerm e k c = LogPuiseuxTerm c e k

-- | Canonical term ordering: ascending lexicographic on @(lpExp, lpLog)@.
termOrd :: LogPuiseuxTerm a -> LogPuiseuxTerm a -> Ordering
termOrd = comparing (\t -> (lpExp t, lpLog t))

-- | Sort a series into canonical order.
normalize :: LogPuiseuxSeries a -> LogPuiseuxSeries a
normalize (LogPuiseuxSeries ts) = LogPuiseuxSeries (sortBy termOrd ts)

-- | Return the leading term (smallest 'lpExp', then smallest 'lpLog'),
-- or 'Nothing' if the series is empty.
leadingTerm :: LogPuiseuxSeries a -> Maybe (LogPuiseuxTerm a)
leadingTerm (LogPuiseuxSeries [])    = Nothing
leadingTerm (LogPuiseuxSeries (t:_)) = Just t

-- | Return the leading term whose coefficient is non-negligible
-- (@|c| > 1e-12@), or 'Nothing' if all coefficients are near zero.
--
-- Used to determine the order of a pole or zero at the expansion point.
leadingTermNZ :: (Floating a, Ord a) => LogPuiseuxSeries a -> Maybe (LogPuiseuxTerm a)
leadingTermNZ (LogPuiseuxSeries ts) =
  case filter (\t -> abs (lpCoeff t) > 1e-12) ts of
    []    -> Nothing
    (t:_) -> Just t

-- | Remove all terms whose coefficient satisfies @|c| ≤ 1e-12@.
--
-- Applied after every arithmetic operation to keep the series in
-- canonical form and prevent accumulation of floating-point noise.
stripZeros :: (Floating a, Ord a) => LogPuiseuxSeries a -> LogPuiseuxSeries a
stripZeros (LogPuiseuxSeries ts) =
  LogPuiseuxSeries $ filter (\t -> abs (lpCoeff t) > 1e-12) ts

-- | Merge two series in canonical order, summing coefficients of
-- terms with equal @(lpExp, lpLog)@ keys.
--
-- Both input series must already be in canonical order (the invariant
-- maintained by all arithmetic operations). The output is normalised
-- and zero-stripped.
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

-- | Add two log-Puiseux series.
addSeries :: (Floating a, Ord a)
          => LogPuiseuxSeries a -> LogPuiseuxSeries a -> LogPuiseuxSeries a
addSeries = combineLike

-- | Scale every term by a constant factor.
scaleSeries :: (Floating a, Ord a)
            => a -> LogPuiseuxSeries a -> LogPuiseuxSeries a
scaleSeries c (LogPuiseuxSeries ts) =
  stripZeros $ LogPuiseuxSeries
    [ LogPuiseuxTerm (c * lpCoeff t) (lpExp t) (lpLog t) | t <- ts ]

-- | Multiply two log-Puiseux series (Cauchy product).
--
-- The product rule for individual terms is:
--
-- @(c₁ · h^p₁ · log^k₁) · (c₂ · h^p₂ · log^k₂)
--     = c₁c₂ · h^(p₁+p₂) · log^(k₁+k₂)@
--
-- The type is closed under multiplication: the log exponent of a
-- product is the sum of the factors' log exponents.
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

-- | Extract the coefficient of the @h^0 · log(h)^0@ term (the true
-- constant term), or @0@ if no such term is present.
constantTerm :: (Floating a, Ord a) => LogPuiseuxSeries a -> a
constantTerm (LogPuiseuxSeries []) = 0
constantTerm (LogPuiseuxSeries (t:_))
  | lpExp t == 0 && lpLog t == 0 = lpCoeff t
  | otherwise                    = 0

-- | Extract all @(logExponent, coefficient)@ pairs at a given power
-- exponent, regardless of log power.
--
-- Used in the @Log@ expansion to recover the full structure at
-- exponent 0 when the constant term may carry a log factor.
constantCoeff :: (Floating a, Ord a) => Rational -> LogPuiseuxSeries a -> [(Int, a)]
constantCoeff e (LogPuiseuxSeries ts) =
  [ (lpLog t, lpCoeff t) | t <- ts, lpExp t == e ]

-- | Remove the term at a specific @(lpExp, lpLog)@ pair.
removeTerm :: Rational -> Int -> LogPuiseuxSeries a -> LogPuiseuxSeries a
removeTerm e k (LogPuiseuxSeries ts) =
  LogPuiseuxSeries $ filter (\t -> not (lpExp t == e && lpLog t == k)) ts

-- | Remove all terms at a given power exponent, across all log powers.
--
-- Used in 'composeSeries' to strip the constant part of an expansion
-- before using it as the perturbation argument.
removeExp :: Rational -> LogPuiseuxSeries a -> LogPuiseuxSeries a
removeExp e (LogPuiseuxSeries ts) =
  LogPuiseuxSeries $ filter (\t -> lpExp t /= e) ts

-- | Retain only terms whose power exponent does not exceed @maxOrder@.
--
-- Used throughout the expansion engine to bound series to a finite
-- working precision.
truncateToOrder :: Rational -> LogPuiseuxSeries a -> LogPuiseuxSeries a
truncateToOrder maxOrder (LogPuiseuxSeries ts) =
  LogPuiseuxSeries $ filter (\t -> lpExp t <= maxOrder) ts

-- | Add @delta@ to every power exponent, leaving log exponents unchanged.
--
-- Used to implement the shift @h^alpha · (...)@ → @(...)@ that
-- factors out the leading power of a series in 'normalizeW',
-- 'invertSeries', and 'expandPowR'.
shiftExponents :: Rational -> LogPuiseuxSeries a -> LogPuiseuxSeries a
shiftExponents delta (LogPuiseuxSeries ts) =
  LogPuiseuxSeries [ LogPuiseuxTerm (lpCoeff t) (lpExp t + delta) (lpLog t)
                   | t <- ts ]

-- | Return 'True' if the series contains any term with @lpLog > 0@.
--
-- A series for which this returns 'False' is a pure Puiseux series
-- (no log factors), and can be handled by the simpler pure-power
-- code paths where applicable.
hasLogTerms :: LogPuiseuxSeries a -> Bool
hasLogTerms (LogPuiseuxSeries ts) = any (\t -> lpLog t > 0) ts