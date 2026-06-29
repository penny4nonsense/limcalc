module LimCalc.Limit where

import LimCalc.Expr
import LimCalc.Types
import LimCalc.Expand
import LimCalc.Puiseux
import qualified Data.Map.Strict as Map

-- | Compute lim_{var → x0} f
-- The limit is the constant term of the Puiseux expansion.
-- If the leading exponent is negative, the limit is infinite (pole).
-- If the series is empty, the limit could not be determined.
limit :: Expr -> String -> Double -> LimitResult Double
limit f var x0 =
  let point = Map.fromList [(var, x0)]
  in case expand f point var of
       Left err -> LimitError err
       Right series ->
         case leadingTermNZ (stripZeros series) of
           Nothing -> Exists 0  -- zero function
           Just lt ->
             if pExp lt < 0
               then Pole (pExp lt)
               else Exists (constantTerm series)

-- | Result of a limit computation
data LimitResult a
  = Exists a                    -- ^ Limit exists and equals a
  | Pole Rational               -- ^ Pole of given order
  | DoesNotExist String         -- ^ Limit does not exist (with reason)
  | LimitError ExpandError      -- ^ Expansion failed
  deriving (Show, Eq)

-- | Compute lim_{var → x0+} f (right-hand limit)
-- The leading term exponent tells us the approach direction
limitRight :: Expr -> String -> Double -> LimitResult Double
limitRight = limit  -- for now same as limit; directional limits TODO

-- | Compute lim_{var → x0-} f (left-hand limit)
limitLeft :: Expr -> String -> Double -> LimitResult Double
limitLeft = limit  -- TODO: directional limits