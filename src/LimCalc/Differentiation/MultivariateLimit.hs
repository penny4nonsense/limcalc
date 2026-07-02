-- | Multivariate limit computation via path-based probing.
--
-- Multivariate limits are computed by substituting a set of
-- parametric paths @(x(t), y(t))@ into the function and evaluating
-- the resulting univariate limit as @t → 0@ using
-- 'LimCalc.Limit.limit'. If all paths agree on a finite value, the
-- limit is reported as existing; if any two paths disagree, the limit
-- is reported as not existing.
--
-- = Completeness caveat
--
-- Agreement on all probed paths does /not/ constitute a proof that
-- the limit exists — it is evidence. A true proof would require
-- showing the limit is the same along /all/ paths, which is generally
-- undecidable. For most elementary functions the standard path set is
-- sufficient in practice. The classic pathological case
-- @f(x,y) = x²y\/(x⁴+y²)@ at the origin — which is 0 on every line
-- but 1\/2 on the parabola @y = x²@ — /is/ caught correctly, because
-- the parabolic path is included in 'standardPaths'.
--
-- = Current limitations
--
-- Only univariate and bivariate limits are implemented. The bivariate
-- case uses 7 standard paths: the two coordinate axes, three lines
-- of slopes 1, −1, 2, the parabola @y − y₀ = (x − x₀)²@, and the
-- cubic @y − y₀ = (x − x₀)³@.
module LimCalc.Differentiation.MultivariateLimit
  ( -- * Result type
    MultivariateLimitResult (..)
    -- * Limit computation
  , multivariateLimit
  , limitAlongPath
    -- * Path construction
  , standardPaths
  , checkPaths
    -- * Expression substitution
  , substExpr
    -- * Utilities
  , withinTolerance
  ) where

import LimCalc.Core.Expr
import LimCalc.Differentiation.Limit
import LimCalc.Core.Simplify
import qualified Data.Map.Strict as Map

-- | Result of a multivariate limit computation.
data MultivariateLimitResult
  = MVExists Double
    -- ^ The limit exists and equals the given value (up to
    -- 'withinTolerance' agreement across all probed paths).
  | MVDoesNotExist String
    -- ^ The limit does not exist; the 'String' gives a reason
    -- (e.g. differing values along different paths, or a pole
    -- along some paths but not others).
  | MVError String
    -- ^ The computation failed (e.g. unsupported number of variables,
    -- or all paths hit singularities).
  deriving (Show, Eq)

-- | Substitute variables in an 'Expr' according to an environment map.
--
-- @substExpr env f@ replaces every @'Var' x@ in @f@ with @env[x]@,
-- recursing into all subexpressions. Variables not in @env@ are
-- left unchanged.
substExpr :: Map.Map String Expr -> Expr -> Expr
substExpr env (Var x)    = Map.findWithDefault (Var x) x env
substExpr _   (Const c)  = Const c
substExpr _   Pi         = Pi
substExpr _   E          = E
substExpr _   I          = I
substExpr env (Add f g)  = Add    (substExpr env f) (substExpr env g)
substExpr env (Sub f g)  = Sub    (substExpr env f) (substExpr env g)
substExpr env (Mul f g)  = Mul    (substExpr env f) (substExpr env g)
substExpr env (Div f g)  = Div    (substExpr env f) (substExpr env g)
substExpr env (Pow f g)  = Pow    (substExpr env f) (substExpr env g)
substExpr env (Neg f)    = Neg    (substExpr env f)
substExpr env (Abs f)    = Abs    (substExpr env f)
substExpr env (Exp f)    = Exp    (substExpr env f)
substExpr env (Log f)    = Log    (substExpr env f)
substExpr env (Sin f)    = Sin    (substExpr env f)
substExpr env (Cos f)    = Cos    (substExpr env f)
substExpr env (Erf f)    = Erf    (substExpr env f)
substExpr env (Li f)     = Li     (substExpr env f)
substExpr env (Si f)     = Si     (substExpr env f)
substExpr env (Ci f)     = Ci     (substExpr env f)
substExpr env (Ei f)     = Ei     (substExpr env f)
substExpr env (Arcsin f) = Arcsin (substExpr env f)
substExpr env (Arccos f) = Arccos (substExpr env f)
substExpr env (Arctan f) = Arctan (substExpr env f)

-- | Compute the limit of @f@ along a parametric path toward a point.
--
-- @path@ maps each variable to an 'Expr' in @\"t\"@ such that each
-- expression evaluates to the target value as @t → 0@. For example,
-- the path @(x, y) = (x₀ + t, y₀ + mt)@ (a line of slope @m@) is:
--
-- @Map.fromList [(\"x\", Add (Const x0) (Var \"t\")), (\"y\", Add (Const y0) (Mul (Const m) (Var \"t\")))]@
--
-- The path is substituted into @f@, simplified, then the univariate
-- limit as @t → 0@ is evaluated.
limitAlongPath :: Expr -> Map.Map String Expr -> LimitResult Double
limitAlongPath f path =
  let f' = simplify (substExpr path f)
  in limit f' "t" 0

-- | Compute @lim_{vars → x0s} f@ via path-based probing.
--
-- Dispatches on the number of variables:
--
-- * 1 variable: delegates directly to 'LimCalc.Differentiation.Limit.limit'.
-- * 2 variables: probes 'standardPaths' and checks consistency via
--   'checkPaths'.
-- * Other: returns 'MVError' (not yet implemented).
multivariateLimit :: Expr -> [String] -> [Double] -> MultivariateLimitResult
multivariateLimit f vars x0s
  | length vars /= length x0s =
      MVError "multivariateLimit: number of variables must match number of target values"
  | length vars == 1 =
      case limit f (head vars) (head x0s) of
        Exists v       -> MVExists v
        Pole _         -> MVDoesNotExist "Pole at limit point"
        DoesNotExist r -> MVDoesNotExist r
        LimitError e   -> MVError (show e)
  | length vars == 2 =
      let [vx, vy] = vars
          [x0, y0] = x0s
      in checkPaths f vx vy x0 y0 (standardPaths vx vy x0 y0)
  | otherwise =
      MVError "multivariateLimit: only 1 and 2 variable cases currently implemented"

-- | The standard set of paths toward @(x₀, y₀)@ used by
-- 'multivariateLimit'.
--
-- Each path is a @Map@ from variable names to 'Expr' in @\"t\"@,
-- representing a curve through @(x₀, y₀)@ as @t → 0@:
--
-- * x-axis: @x = x₀ + t, y = y₀@
-- * y-axis: @x = x₀, y = y₀ + t@
-- * Line slope 1: @x = x₀ + t, y = y₀ + t@
-- * Line slope −1: @x = x₀ + t, y = y₀ − t@
-- * Line slope 2: @x = x₀ + t, y = y₀ + 2t@
-- * Parabola: @x = x₀ + t, y = y₀ + t²@
-- * Cubic: @x = x₀ + t, y = y₀ + t³@
standardPaths :: String -> String -> Double -> Double -> [Map.Map String Expr]
standardPaths vx vy x0 y0 =
  let cx = Const x0
      cy = Const y0
      t  = Var "t"
  in [ Map.fromList [(vx, Add cx t), (vy, cy)]
     , Map.fromList [(vx, cx), (vy, Add cy t)]
     , Map.fromList [(vx, Add cx t), (vy, Add cy t)]
     , Map.fromList [(vx, Add cx t), (vy, Sub cy t)]
     , Map.fromList [(vx, Add cx t), (vy, Add cy (Mul (Const 2) t))]
     , Map.fromList [(vx, Add cx t), (vy, Add cy (Pow t (Const 2)))]
     , Map.fromList [(vx, Add cx t), (vy, Add cy (Pow t (Const 3)))]
     ]

-- | Evaluate limits along each path and check for consistency.
--
-- * If all non-error paths give the same finite value → 'MVExists'.
-- * If any two finite-valued paths disagree → 'MVDoesNotExist'.
-- * If some paths give finite values and others give poles →
--   'MVDoesNotExist'.
-- * If all paths hit singularities or errors → 'MVError'.
--
-- 'LimitError' results along individual paths are skipped: a path
-- may pass through a point where the function is undefined
-- (e.g. @sin(xy)\/(xy)@ along @y = 0@) even though the overall
-- limit exists. Only the finite and pole results determine the outcome.
checkPaths :: Expr -> String -> String -> Double -> Double
           -> [Map.Map String Expr] -> MultivariateLimitResult
checkPaths f _vx _vy _x0 _y0 paths =
  let results = map (limitAlongPath f) paths
      finite  = [ v | Exists v <- results ]
      poles   = [ r | Pole r   <- results ]
      errors  = [ e | LimitError e <- results ]
  in case (finite, poles) of
       ([], []) -> MVError ("No valid path results; all paths hit " ++
                            "singularities or errors. First error: " ++
                            show (head errors))
       ([], _)  -> MVDoesNotExist "Pole along at least one path"
       _ ->
         case poles of
           (_:_) -> MVDoesNotExist
                      "Limit does not exist: finite along some paths, \
                      \pole along others"
           [] ->
             let v  = head finite
                 vs = tail finite
             in case filter (not . withinTolerance v) vs of
                  []               -> MVExists v
                  (disagreeing:_)  -> MVDoesNotExist $
                    "Limit does not exist: different values along \
                    \different paths (e.g. " ++ show v ++
                    " vs " ++ show disagreeing ++ ")"

-- | Tolerance for comparing limit values from different paths.
-- Two values are considered equal if @|a − b| < 1e-6@.
withinTolerance :: Double -> Double -> Bool
withinTolerance a b = abs (a - b) < 1e-6