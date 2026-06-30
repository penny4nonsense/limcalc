module LimCalc.Types where

import Data.Map.Strict (Map)
import LimCalc.Puiseux
import LimCalc.AlgNum

-- | A point in n-dimensional space
type Point = Map String AlgNum

-- | Errors during expansion
data ExpandError
  = Singularity String
  | Undefined String
  | DomainError String
  | Unknown String
  -- | A point where the function provably has no local Puiseux
  -- series representation -- a genuine analytic kink (e.g. |x| at
  -- x=0), as opposed to Unknown (couldn't determine) or Singularity
  -- (pole/blow-up). This is a definite, certified answer: the
  -- engine knows no derivative exists here, not that it failed to
  -- find one.
  | NonAnalytic String
  deriving (Show, Eq)

-- | Result of expansion
type ExpandResult = Either ExpandError (PuiseuxSeries AlgNum)