module LimCalc.Calculus where

import LimCalc.Expr
import LimCalc.Puiseux
import LimCalc.Types
import LimCalc.Expand
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

-- | Symbolic partial derivative of f with respect to var, treating
-- all other variables as constants.
--
-- Uses DiffField.deriveExpr directly (the algebraic chain-rule
-- approach) rather than symExpand. This correctly handles iterated
-- partial derivatives (e.g. d^2f/dx^2 = d/dx(df/dx)) because
-- deriveExpr never introduces spurious division or expansion
-- artifacts -- it treats every Var other than var as a constant and
-- applies the usual algebraic rules. The series-expansion approach
-- (symExpand) was incorrect for iterated partials because expanding
-- the already-differentiated result around x again introduced x0
-- terms that polluted the coefficient extraction.
partialDiff :: Expr -> String -> Either ExpandError Expr
partialDiff f var =
  Right $ simplify $ deriveExpr (baseField var) f

-- | Gradient of a scalar-valued expression with respect to a list of
-- variables. Returns [df/dx1, df/dx2, ..., df/dxn].
gradient :: Expr -> [String] -> Either ExpandError [Expr]
gradient f vars = mapM (partialDiff f) vars

-- | Jacobian matrix of a vector-valued function with respect to a
-- list of variables. Entry (i,j) = dfi/dxj.
-- Returns a list of rows, each row being the gradient of one fi.
jacobian :: [Expr] -> [String] -> Either ExpandError [[Expr]]
jacobian fs vars = mapM (\f -> gradient f vars) fs

-- | Hessian matrix of a scalar-valued expression.
-- Entry (i,j) = d^2f/(dxi dxj). Symmetric by Clairaut's theorem
-- (for smooth f), so H[i][j] = H[j][i].
hessian :: Expr -> [String] -> Either ExpandError [[Expr]]
hessian f vars = do
  grad <- gradient f vars
  mapM (\dfdxi -> gradient dfdxi vars) grad