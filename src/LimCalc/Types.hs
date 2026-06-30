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
  deriving (Show, Eq)

-- | Result of expansion
type ExpandResult = Either ExpandError (PuiseuxSeries AlgNum)