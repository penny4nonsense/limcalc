module LimCalc.Risch.Exponential where

import LimCalc.Poly
import LimCalc.RationalFunction
import LimCalc.DiffField
import LimCalc.Expr
import LimCalc.Risch.Primitive
import LimCalc.Simplify

-- | Result of the Risch algorithm for the exponential case
data ExponentialResult
  = ExponentialElementary Expr      -- ^ Elementary antiderivative
  | ExponentialNonElementary        -- ^ Provably non-elementary
  | ExponentialError String         -- ^ Algorithm error
  deriving (Show, Eq)

-- | Integrate in the exponential case
-- Field extension: F(θ) where θ = exp(f), Dθ = θ·Df
-- The key difference from primitive: Dθ/θ = Df is in F, not F(θ)
integrateExponential :: RatFun -> DiffField -> ExponentialResult
integrateExponential rf field =
  -- Step 1: Hermite reduction
  let (g, reduced) = hermiteReduce rf
  -- Step 2: Handle the reduced part
      result = integrateExponentialReduced reduced field
  in case result of
       Left err    -> ExponentialError err
       Right terms ->
         let logPart = foldr addExpr (Const 0)
               [ Mul (Const c) (Log (polyToExpr u))
               | (c, u) <- terms ]
             gExpr = ratFunToExpr g
         in ExponentialElementary (addExpr gExpr logPart)

-- | Handle the reduced part in the exponential case
-- For the exponential case, we need to handle polynomials in θ specially
integrateExponentialReduced :: RatFun -> DiffField -> Either String [(Double, Poly)]
integrateExponentialReduced (RatFun a d) field =
  -- In the exponential case, if d has no factors involving θ,
  -- we can use the same Rothstein-Trager approach
  -- Otherwise we need the LODE (linear ODE) approach
  let d' = diffPoly d
      g  = gcdPoly d d'
  in if degree g == 0
     -- d is squarefree, use Rothstein-Trager
     then rothsteinTragerExp a d field
     -- d has repeated factors, need further reduction
     else Left "Repeated factors in exponential case — not yet implemented"

-- | Rothstein-Trager for the exponential case
rothsteinTragerExp :: Poly -> Poly -> DiffField -> Either String [(Double, Poly)]
rothsteinTragerExp a d field =
  let d'    = diffPoly d
      rPoly = resultantPoly a d d'
  in case findRationalRoots rPoly of
       Nothing    -> Left "Non-rational logarithmic constants — non-elementary"
       Just roots ->
         Right [ (c, gcdPoly d (subPoly a (scalePoly c d')))
               | c <- roots
               , c /= 0
               ]

-- | Integrate exp(f) directly
-- int exp(f) = exp(f)/f' if f' is constant, else non-elementary in general
integrateExp :: Expr -> String -> Either String Expr
integrateExp f var =
  let f'      = simplify (deriveBase f)
      isConst = not (containsVar var f')
  in if isConst
     then case f' of
            Const 0 -> Left "Derivative of exponent is zero"
            Const c -> Right $ Div (Exp f) (Const c)
            _       -> Left "Non-constant derivative — non-elementary"
     else Left "Non-elementary exponential integral"

-- | Check if an expression contains a variable
containsVar :: String -> Expr -> Bool
containsVar v (Var x)     = v == x
containsVar v (Const _)   = False
containsVar v Pi          = False
containsVar v E           = False
containsVar v I           = False
containsVar v (Add f g)   = containsVar v f || containsVar v g
containsVar v (Sub f g)   = containsVar v f || containsVar v g
containsVar v (Mul f g)   = containsVar v f || containsVar v g
containsVar v (Div f g)   = containsVar v f || containsVar v g
containsVar v (Pow f g)   = containsVar v f || containsVar v g
containsVar v (Neg f)     = containsVar v f
containsVar v (Exp f)     = containsVar v f
containsVar v (Log f)     = containsVar v f
containsVar v (Sin f)     = containsVar v f
containsVar v (Cos f)     = containsVar v f
containsVar v (Abs f)     = containsVar v f