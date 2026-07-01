module LimCalc.Calculus where
import LimCalc.Expr
import LimCalc.Puiseux
import LimCalc.Types
import LimCalc.Expand hiding (factorial)
import LimCalc.AlgNum
import LimCalc.SymExpand
import LimCalc.SymPuiseux
import LimCalc.DiffField
import LimCalc.Simplify
import qualified Data.Map.Strict as Map

-- | Numeric derivative of f with respect to var at point
derivative :: Expr -> Map.Map String Double -> String -> Either ExpandError Double
derivative f pointD var =
  let point = Map.map (fromQ . toRational) pointD
  in case expand f point var of
       Left err     -> Left err
       Right series -> Right $ algToDouble (coeffAt 1 series)

-- | Extract the coefficient of h^n * log(h)^0 from a log-Puiseux series
coeffAt :: Rational -> LogPuiseuxSeries AlgNum -> AlgNum
coeffAt n (LogPuiseuxSeries ts) =
  case filter (\t -> lpExp t == n && lpLog t == 0) ts of
    []    -> algZero
    (t:_) -> lpCoeff t

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

-- | Symbolic partial derivative of f with respect to var, treating
-- all other variables as constants.
partialDiff :: Expr -> String -> Either ExpandError Expr
partialDiff f var =
  Right $ simplify $ deriveExpr (baseField var) f

-- | Gradient of a scalar-valued expression with respect to a list of variables.
gradient :: Expr -> [String] -> Either ExpandError [Expr]
gradient f vars = mapM (partialDiff f) vars

-- | Jacobian matrix of a vector-valued function with respect to a list of variables.
jacobian :: [Expr] -> [String] -> Either ExpandError [[Expr]]
jacobian fs vars = mapM (\f -> gradient f vars) fs

-- | Hessian matrix of a scalar-valued expression.
hessian :: Expr -> [String] -> Either ExpandError [[Expr]]
hessian f vars = do
  grad <- gradient f vars
  mapM (\dfdxi -> gradient dfdxi vars) grad