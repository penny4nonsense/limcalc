module LimCalc.Limit where

import LimCalc.Expr
import LimCalc.Types
import LimCalc.Expand
import LimCalc.Puiseux
import LimCalc.AlgNum
import qualified Data.Map.Strict as Map

-- | Compute lim_{var → x0} f
limit :: Expr -> String -> Double -> LimitResult Double
limit f var x0 =
  let point = Map.fromList [(var, fromQ (toRational x0))]
  in case expand f point var of
       Left err -> LimitError err
       Right series ->
         case leadingTermNZ (stripZeros series) of
           Nothing -> Exists 0
           Just lt ->
             if pExp lt < 0
               then Pole (pExp lt)
               else Exists (algToDouble (constantTerm series))

-- | Result of a limit computation
data LimitResult a
  = Exists a
  | Pole Rational
  | DoesNotExist String
  | LimitError ExpandError
  deriving (Show, Eq)

-- | Right-hand limit
limitRight :: Expr -> String -> Double -> LimitResult Double
limitRight = limit

-- | Left-hand limit
limitLeft :: Expr -> String -> Double -> LimitResult Double
limitLeft = limit