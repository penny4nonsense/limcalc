module LimCalc.MultivariateLimit where

import LimCalc.Expr
import LimCalc.Limit
import LimCalc.Simplify
import qualified Data.Map.Strict as Map

-- | Result of a multivariate limit computation
data MultivariateLimitResult
  = MVExists Double          -- ^ Limit exists and equals this value
  | MVDoesNotExist String    -- ^ Limit does not exist (reason given)
  | MVError String           -- ^ Could not compute (technical failure)
  deriving (Show, Eq)

-- | Substitute variables in an Expr with sub-expressions.
-- substExpr env f replaces every Var x in f with env[x], recursing
-- into all sub-expressions. Variables not in env are left unchanged.
substExpr :: Map.Map String Expr -> Expr -> Expr
substExpr env (Var x)     = Map.findWithDefault (Var x) x env
substExpr _   (Const c)   = Const c
substExpr _   Pi          = Pi
substExpr _   E           = E
substExpr _   I           = I
substExpr env (Add f g)   = Add (substExpr env f) (substExpr env g)
substExpr env (Sub f g)   = Sub (substExpr env f) (substExpr env g)
substExpr env (Mul f g)   = Mul (substExpr env f) (substExpr env g)
substExpr env (Div f g)   = Div (substExpr env f) (substExpr env g)
substExpr env (Pow f g)   = Pow (substExpr env f) (substExpr env g)
substExpr env (Neg f)     = Neg (substExpr env f)
substExpr env (Abs f)     = Abs (substExpr env f)
substExpr env (Exp f)     = Exp (substExpr env f)
substExpr env (Log f)     = Log (substExpr env f)
substExpr env (Sin f)     = Sin (substExpr env f)
substExpr env (Cos f)     = Cos (substExpr env f)
substExpr env (Erf f)     = Erf (substExpr env f)
substExpr env (Li f)      = Li  (substExpr env f)
substExpr env (Si f)      = Si  (substExpr env f)
substExpr env (Ci f)      = Ci  (substExpr env f)
substExpr env (Ei f)      = Ei  (substExpr env f)
substExpr env (Arcsin f)  = Arcsin (substExpr env f)
substExpr env (Arccos f)  = Arccos (substExpr env f)
substExpr env (Arctan f)  = Arctan (substExpr env f)

-- | Evaluate the limit of f along a parametric path toward a point.
--
-- path: a map from each variable to an Expr in "t" such that each
-- expression evaluates to x0[var] as t → 0. E.g. for the path
-- (x,y)=(x0+t, y0+m*t) (a line with slope m through (x0,y0)):
--   path = {"x": Const x0 + Var "t", "y": Const y0 + Mul (Const m) (Var "t")}
--
-- Substitutes the path into f, simplifies, then takes t → 0 using
-- the existing univariate limit machinery.
limitAlongPath :: Expr -> Map.Map String Expr -> LimitResult Double
limitAlongPath f path =
  let f' = simplify (substExpr path f)
  in limit f' "t" 0

-- | Compute the multivariate limit of f as (vars) → (x0s).
--
-- Strategy: probe a standard set of paths through the target point.
-- If all paths give the same finite value → MVExists that value.
-- If any two paths give different values → MVDoesNotExist.
-- If a path gives a pole or error → reported accordingly.
--
-- The path set (for 2 variables x,y → (x0,y0)) covers:
--   - Along each coordinate axis
--   - Along lines y-y0 = m*(x-x0) for m in {1,-1,2}
--   - Along the curve y-y0 = (x-x0)^2
-- These catch most common cases of non-existence (e.g. xy/(x^2+y^2)
-- at the origin gives different values along y=0 and y=x).
--
-- NOTE: Agreement on all probed paths does NOT constitute a proof
-- that the limit exists -- it is evidence. A true proof would
-- require showing the limit along ALL paths, which is generally
-- undecidable. For most elementary functions the probed paths are
-- sufficient, but pathological cases exist (e.g. the classic
-- f(x,y) = x^2*y/(x^4+y^2) at the origin is 0 on every line but
-- 1/2 on the parabola y=x^2 -- our path set does include y=x^2,
-- so this specific case IS caught correctly).
multivariateLimit :: Expr -> [String] -> [Double] -> MultivariateLimitResult
multivariateLimit f vars x0s
  | length vars /= length x0s =
      MVError "multivariateLimit: number of variables must match number of target values"
  | length vars == 1 =
      -- Delegate to the univariate limit directly
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

-- | The standard set of paths toward (x0, y0).
-- Each path is a Map from variable names to Expr in "t",
-- representing a curve through (x0,y0) as t→0.
standardPaths :: String -> String -> Double -> Double -> [Map.Map String Expr]
standardPaths vx vy x0 y0 =
  let cx = Const x0
      cy = Const y0
      t  = Var "t"
  in [ -- Along x-axis: x = x0+t, y = y0
       Map.fromList [(vx, Add cx t), (vy, cy)]
       -- Along y-axis: x = x0, y = y0+t
     , Map.fromList [(vx, cx), (vy, Add cy t)]
       -- Along line y-y0 = x-x0 (slope 1)
     , Map.fromList [(vx, Add cx t), (vy, Add cy t)]
       -- Along line y-y0 = -(x-x0) (slope -1)
     , Map.fromList [(vx, Add cx t), (vy, Sub cy t)]
       -- Along line y-y0 = 2*(x-x0) (slope 2)
     , Map.fromList [(vx, Add cx t), (vy, Add cy (Mul (Const 2) t))]
       -- Along parabola y-y0 = (x-x0)^2
     , Map.fromList [(vx, Add cx t), (vy, Add cy (Pow t (Const 2)))]
       -- Along cubic y-y0 = (x-x0)^3
     , Map.fromList [(vx, Add cx t), (vy, Add cy (Pow t (Const 3)))]
     ]

-- | Evaluate limits along each path and check for consistency.
-- Returns MVExists v if all non-error paths give v,
-- MVDoesNotExist if any two disagree.
checkPaths :: Expr -> String -> String -> Double -> Double
           -> [Map.Map String Expr] -> MultivariateLimitResult
checkPaths f vx vy x0 y0 paths =
  let results = map (limitAlongPath f) paths
      finite  = [ v | Exists v <- results ]
      poles   = [ r | Pole r   <- results ]
      -- LimitErrors (singularities) along individual paths are
      -- skipped: a path may pass through a point where the function
      -- is undefined (e.g. sin(xy)/(xy) along y=0 gives sin(0)/0
      -- which is a "singularity" for the expansion engine even
      -- though the limit exists). Only error out if NO valid paths.
      errors  = [ e | LimitError e <- results ]
  in case (finite, poles) of
       ([], []) -> MVError ("No valid path results; all paths hit " ++
                            "singularities or errors. First error: " ++
                            show (head errors))
       ([], _)  -> MVDoesNotExist "Pole along at least one path"
       _        ->
         case poles of
           (_:_) -> MVDoesNotExist
                      "Limit does not exist: finite along some paths, \
                      \pole along others"
           []    ->
             let v = head finite
                 vs = tail finite
             in case filter (not . withinTolerance v) vs of
                  [] -> MVExists v
                  (disagreeing:_) -> MVDoesNotExist $
                    "Limit does not exist: different values along \
                    \different paths (e.g. " ++ show v ++
                    " vs " ++ show disagreeing ++ ")"

-- | Tolerance for comparing limit values from different paths.
-- Two values are considered equal if they agree to 6 decimal places.
withinTolerance :: Double -> Double -> Bool
withinTolerance a b = abs (a - b) < 1e-6
