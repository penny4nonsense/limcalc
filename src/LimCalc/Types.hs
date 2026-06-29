module LimCalc.Types where

import Data.Map.Strict (Map)
import LimCalc.Puiseux

-- | A point in n-dimensional space, mapping variable names to values
type Point = Map String Double

-- | Errors that can occur during series expansion
data ExpandError
  = Singularity String    -- ^ Pole or essential singularity
  | Undefined String      -- ^ Log of negative, sqrt of negative, etc.
  | DomainError String    -- ^ Would require complex result
  | Unknown String        -- ^ Couldn't determine
  deriving (Show, Eq)

-- | Result of an expansion attempt
type ExpandResult = Either ExpandError PuiseuxSeries