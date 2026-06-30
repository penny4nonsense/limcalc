module LimCalc.Calculus where

import LimCalc.Expr
import LimCalc.Puiseux
import LimCalc.Types
import LimCalc.Expand
import LimCalc.AlgNum
import LimCalc.SymExpand
import LimCalc.SymPuiseux
import qualified Data.Map.Strict as Map

-- | Numeric derivative of f with respect to var at point
derivative :: Expr -> Map.Map String Double -> String -> Either ExpandError Double
derivative f pointD var =
  let point = Map.map (fromQ . toRational) pointD
  in case expand f point var of
       Left err     -> Left err
       Right series -> Right $ algToDouble (coeffAt 1 series)

-- | Extract the coefficient of h^n from a series
coeffAt :: Rational -> PuiseuxSeries AlgNum -> AlgNum
coeffAt n (PuiseuxSeries ts) =
  case filter (\t -> pExp t == n) ts of
    []    -> algZero
    (t:_) -> coeff t

-- | nth derivative
nthDerivative :: Int -> Expr -> Map.Map String Double -> String -> Either ExpandError Double
nthDerivative n f pointD var =
  let point = Map.map (fromQ . toRational) pointD
  in case expand f point var of
       Left err     -> Left err
       Right series -> Right $ algToDouble $
         fromIntegral (factorial n) * coeffAt (fromIntegral n) series

-- | Factorial
factorial :: Int -> Int
factorial 0 = 1
factorial n = n * factorial (n-1)

-- | Partial derivative
partial :: Expr -> Map.Map String Double -> String -> Either ExpandError Double
partial = derivative

-- | Symbolic derivative
diff :: Expr -> String -> Either ExpandError Expr
diff f var = do
  series <- symExpand f var
  return $ symCoeffAt 1 series