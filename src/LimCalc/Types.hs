-- | Shared types used across the limcalc expansion and calculus pipeline.
--
-- This module exists to break import cycles: 'ExpandResult' is needed
-- by both 'LimCalc.Expand' and 'LimCalc.Calculus', and 'Point' is
-- needed by 'LimCalc.Expand' and 'LimCalc.MultivariateLimit'. Placing
-- them here lets all consumers import a single lightweight module
-- without pulling in the full expansion engine.
module LimCalc.Types
  ( -- * Points
    Point
    -- * Expansion errors
  , ExpandError (..)
    -- * Expansion results
  , ExpandResult
  ) where

import Data.Map.Strict (Map)
import LimCalc.Puiseux
import LimCalc.AlgNum

-- | A point in a (possibly multivariate) domain: a map from variable
-- names to their 'AlgNum' values.
--
-- Used as the base point @x₀@ in 'LimCalc.Expand.expand': the
-- expansion computes the log-Puiseux series of @f(x₀ + h)@ in @h@.
type Point = Map String AlgNum

-- | Errors that can arise during log-Puiseux series expansion.
--
-- The constructors are ordered roughly by informativeness:
--
-- * 'Singularity' and 'Undefined' signal that the function is
--   genuinely not defined or not expandable at the given point.
-- * 'NonAnalytic' is a /certified/ negative result: the engine has
--   proved that no local Puiseux series exists (e.g. @|x|@ at @x=0@).
-- * 'DomainError' signals that the input is outside the function's
--   domain in a way that prevents expansion.
-- * 'Unknown' signals that the engine could not determine the
--   expansion — a gap in coverage, not a mathematical impossibility.
data ExpandError
  = Singularity String
    -- ^ The function has a pole or essential singularity at the
    -- expansion point (e.g. @1\/x@ at @x=0@). The series is still
    -- computed — with negative leading exponent — but this error is
    -- returned when the singularity makes further processing
    -- (e.g. division by zero series) impossible.
  | Undefined String
    -- ^ The function is not defined at the expansion point
    -- (e.g. @log(x)@ at @x ≤ 0@, @log@ of a series with a pole).
  | DomainError String
    -- ^ The input is outside the function's natural domain in a way
    -- that is distinct from a pole or logarithmic singularity
    -- (e.g. @|f|@ when @f@ has a non-real leading coefficient).
  | Unknown String
    -- ^ The expansion engine does not have an implementation for
    -- this case. This is a coverage gap, not a mathematical result:
    -- the function may well be expandable, but the engine cannot
    -- currently do it (e.g. @Li(x)@ near @x=0@, which requires
    -- @log(log(h))@).
  | NonAnalytic String
    -- ^ The function provably has no local log-Puiseux series
    -- representation at this point. This is a certified negative
    -- result: the engine has determined that no derivative exists
    -- here. Distinct from 'Unknown' (which means the engine gave up)
    -- and 'Singularity' (which means a pole). Example: @|x|@ at
    -- @x=0@, where the function has a genuine kink.
  deriving (Show, Eq)

-- | The result of a log-Puiseux series expansion: either an
-- 'ExpandError' explaining why expansion failed, or a
-- 'LogPuiseuxSeries' over 'AlgNum' coefficients.
type ExpandResult = Either ExpandError (LogPuiseuxSeries AlgNum)