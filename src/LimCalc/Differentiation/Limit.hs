-- | Univariate limit computation via log-Puiseux series expansion.
--
-- Limits are computed by expanding @f(x₀ + h)@ as a log-Puiseux
-- series in @h@ via 'LimCalc.Series.Expand.expand', then reading off
-- the behaviour as @h → 0⁺@:
--
-- * If the series is identically zero (all terms vanish), the limit is 0.
-- * If the leading term has negative exponent, the function has a pole.
-- * Otherwise, the limit is the constant term (the @h^0 · log(h)^0@
--   coefficient).
--
-- Note that terms with @lpLog > 0@ (e.g. @log(h)@ terms from @Ci@
-- or @Ei@ near 0) are handled correctly: @log(h) → −∞@ as @h → 0⁺@,
-- so such a term dominates over the constant. However, the current
-- implementation reads off the constant term of the series regardless
-- of log terms; for functions like @Ci(x)@ near @x = 0@, the
-- @gamma + log(h)@ structure means the "limit" diverges, and the
-- expansion will reflect this via the log term's coefficient.
--
-- = One-sided limits
--
-- 'limitRight' and 'limitLeft' are currently identical to 'limit'
-- (both expand around @x₀@ from the right). True left-hand limits
-- would require expanding around @x₀ - h@ instead; this is a known
-- gap.
module LimCalc.Differentiation.Limit
  ( -- * Limit computation
    limit
  , limitRight
  , limitLeft
    -- * Result type
  , LimitResult (..)
  ) where

import LimCalc.Core.Expr
import LimCalc.Core.Types
import LimCalc.Series.Expand
import LimCalc.Series.Puiseux
import LimCalc.Algebra.AlgNum
import qualified Data.Map.Strict as Map

-- | Compute @lim_{var → x₀} f@ via log-Puiseux series expansion.
--
-- Returns:
--
-- * @'Exists' v@ — the limit exists and equals @v@.
-- * @'Pole' p@ — the function has a pole of order @p@ at @x₀@
--   (leading exponent @p < 0@).
-- * @'LimitError' err@ — the expansion failed with 'ExpandError' @err@.
--   This includes 'NonAnalytic' for functions like @|x|@ at @x = 0@,
--   and 'Unknown' for unimplemented cases.
-- * @'DoesNotExist'@ — not currently produced by the univariate engine
--   (used by 'LimCalc.MultivariateLimit').
limit :: Expr -> String -> Double -> LimitResult Double
limit f var x0 =
  let point = Map.fromList [(var, fromQ (toRational x0))]
  in case expand f point var of
       Left err    -> LimitError err
       Right series ->
         case leadingTermNZ (stripZeros series) of
           Nothing -> Exists 0
           Just lt ->
             if lpExp lt < 0
               then Pole (lpExp lt)
               else Exists (algToDouble (constantTerm series))

-- | The result of a limit computation.
data LimitResult a
  = Exists a
    -- ^ The limit exists and equals the given value.
  | Pole Rational
    -- ^ The function has a pole at the limit point; the 'Rational'
    -- is the leading exponent of the Laurent\/Puiseux expansion
    -- (negative for a pole).
  | DoesNotExist String
    -- ^ The limit does not exist. The 'String' gives a reason.
    -- Currently only produced by 'LimCalc.MultivariateLimit'.
  | LimitError ExpandError
    -- ^ The expansion engine returned an error. This includes
    -- 'NonAnalytic' for genuine kinks and 'Unknown' for unimplemented
    -- cases.
  deriving (Show, Eq)

-- | Right-hand limit @lim_{var → x₀⁺} f@.
--
-- Currently identical to 'limit': the expansion is always from the
-- right (h → 0⁺). True directional limits are a known gap.
limitRight :: Expr -> String -> Double -> LimitResult Double
limitRight = limit

-- | Left-hand limit @lim_{var → x₀⁻} f@.
--
-- Currently identical to 'limit'. A genuine left-hand limit would
-- require expanding around @x₀ − h@ instead of @x₀ + h@.
limitLeft :: Expr -> String -> Double -> LimitResult Double
limitLeft = limit