-- | Calculus operations: differentiation (numeric and symbolic),
-- and multivariate calculus (gradient, Jacobian, Hessian).
--
-- = Two differentiation paths
--
-- limcalc provides two symbolic differentiation functions:
--
-- * 'diff': the /series extraction path/. Expands @f(x + h)@ as a
--   'LimCalc.Series.SymPuiseux.SymPuiseuxSeries' via
--   'LimCalc.Series.SymExpand.symExpand' and reads off the @h^1@
--   coefficient. Produces a symbolic 'Expr' for the derivative.
--   Used for univariate symbolic differentiation.
--
-- * 'partialDiff': the /algebraic chain-rule path/. Differentiates
--   by structural recursion via
--   'LimCalc.Differentiation.DiffField.deriveExpr'. Correct for
--   iterated partial derivatives: differentiating the result of a
--   previous 'partialDiff' works because 'deriveExpr' treats all
--   variables other than the specified one as constants, with no
--   expansion artifacts.
--
-- The series path ('diff') generates Taylor coefficients by iterating
-- 'LimCalc.Differentiation.DiffField.deriveBase', making it a
-- consequence of the algebraic path rather than an independent
-- implementation.
--
-- = Numeric differentiation
--
-- 'derivative' and 'nthDerivative' evaluate derivatives numerically
-- at a specific point by expanding @f(x₀ + h)@ as a
-- 'LimCalc.Series.Puiseux.LogPuiseuxSeries' via
-- 'LimCalc.Series.Expand.expand' and reading off the appropriate
-- coefficient.
module LimCalc.Differentiation.Calculus
  ( -- * Numeric differentiation
    derivative
  , nthDerivative
  , partial
    -- * Symbolic differentiation
  , diff
  , partialDiff
    -- * Multivariate calculus
  , gradient
  , jacobian
  , hessian
    -- * Helpers
  , coeffAt
  , factorial
  ) where

import LimCalc.Core.Expr
import LimCalc.Series.Puiseux
import LimCalc.Core.Types
import LimCalc.Series.Expand hiding (factorial)
import LimCalc.Algebra.AlgNum
import LimCalc.Series.SymExpand
import LimCalc.Series.SymPuiseux
import LimCalc.Differentiation.DiffField
import LimCalc.Core.Simplify
import qualified Data.Map.Strict as Map

-- | Numerically evaluate the derivative of @f@ with respect to @var@
-- at a specific point.
--
-- Expands @f(x₀ + h)@ as a log-Puiseux series and returns the
-- coefficient of @h^1 · log(h)^0@ as a 'Double'.
derivative :: Expr -> Map.Map String Double -> String -> Either ExpandError Double
derivative f pointD var =
  let point = Map.map (fromQ . toRational) pointD
  in case expand f point var of
       Left err     -> Left err
       Right series -> Right $ algToDouble (coeffAt 1 series)

-- | Extract the coefficient of @h^n · log(h)^0@ from a log-Puiseux series.
--
-- Returns 'algZero' if no such term exists. Used by 'derivative' and
-- 'nthDerivative' to read off Taylor coefficients from the expansion.
coeffAt :: Rational -> LogPuiseuxSeries AlgNum -> AlgNum
coeffAt n (LogPuiseuxSeries ts) =
  case filter (\t -> lpExp t == n && lpLog t == 0) ts of
    []    -> algZero
    (t:_) -> lpCoeff t

-- | Numerically evaluate the @n@th derivative of @f@ at a specific point.
--
-- Uses the relation @f^(n)(x₀) = n! · [h^n] f(x₀ + h)@, where
-- @[h^n]@ denotes the coefficient of @h^n · log(h)^0@ in the
-- log-Puiseux expansion.
nthDerivative :: Int -> Expr -> Map.Map String Double -> String -> Either ExpandError Double
nthDerivative n f pointD var =
  let point = Map.map (fromQ . toRational) pointD
  in case expand f point var of
       Left err     -> Left err
       Right series -> Right $ algToDouble $
         fromIntegral (factorial n) * coeffAt (fromIntegral n) series

-- | Factorial function used by 'nthDerivative'.
factorial :: Int -> Int
factorial 0 = 1
factorial n = n * factorial (n-1)

-- | Numeric partial derivative. Alias for 'derivative'.
partial :: Expr -> Map.Map String Double -> String -> Either ExpandError Double
partial = derivative

-- | Symbolic derivative of @f@ with respect to @var@.
--
-- Uses the series extraction path: expands @f(x + h)@ symbolically
-- via 'LimCalc.SymExpand.symExpand' and returns the @h^1@ coefficient
-- as an 'Expr'. The result is not automatically simplified; apply
-- 'LimCalc.Simplify.simplify' as needed.
diff :: Expr -> String -> Either ExpandError Expr
diff f var = do
  series <- symExpand f var
  return $ symCoeffAt 1 series

-- | Symbolic partial derivative of @f@ with respect to @var@,
-- treating all other variables as constants.
--
-- Uses the algebraic chain-rule path via
-- 'LimCalc.DiffField.deriveExpr'. This correctly handles iterated
-- partial derivatives: the result of a previous 'partialDiff' can
-- be differentiated again without introducing expansion artifacts.
-- The result is simplified via 'LimCalc.Simplify.simplify'.
partialDiff :: Expr -> String -> Either ExpandError Expr
partialDiff f var =
  Right $ simplify $ deriveExpr (baseField var) f

-- | Gradient of a scalar-valued expression: @[∂f\/∂x₁, …, ∂f\/∂xₙ]@.
gradient :: Expr -> [String] -> Either ExpandError [Expr]
gradient f vars = mapM (partialDiff f) vars

-- | Jacobian matrix of a vector-valued function.
--
-- Entry @(i, j)@ is @∂fᵢ\/∂xⱼ@. Returns a list of rows, each row
-- being the gradient of one component @fᵢ@.
jacobian :: [Expr] -> [String] -> Either ExpandError [[Expr]]
jacobian fs vars = mapM (\f -> gradient f vars) fs

-- | Hessian matrix of a scalar-valued expression.
--
-- Entry @(i, j)@ is @∂²f\/(∂xᵢ ∂xⱼ)@. Symmetric by Clairaut's
-- theorem for smooth @f@.
hessian :: Expr -> [String] -> Either ExpandError [[Expr]]
hessian f vars = do
  grad <- gradient f vars
  mapM (\dfdxi -> gradient dfdxi vars) grad