module LimCalc.Risch.Exponential where

import LimCalc.Poly
import LimCalc.RationalFunction
import LimCalc.DiffField
import LimCalc.Expr
import LimCalc.Simplify
import LimCalc.Risch.Primitive

-- | Result of the Risch algorithm for the exponential case
data ExponentialResult
  = ExponentialElementary Expr
  | ExponentialNonElementary
  | ExponentialError String
  deriving (Show, Eq)

-- | Integrate in the exponential case
integrateExponential :: RatFun Double -> DiffField -> ExponentialResult
integrateExponential rf field =
  let (g, reduced) = hermiteReduce rf
      result       = integrateExponentialReduced reduced field
  in case result of
       Left err    -> ExponentialError err
       Right terms ->
         let logPart = foldr addExpr (Const 0)
               [ Mul (Const c) (Log (polyToExpr u))
               | (c, u) <- terms ]
             gExpr = ratFunToExpr g
         in ExponentialElementary (addExpr gExpr logPart)

-- | Handle the reduced part in the exponential case
integrateExponentialReduced :: RatFun Double -> DiffField -> Either String [(Double, Poly Double)]
integrateExponentialReduced (RatFun a d) field =
  let d' = diffPoly d
      g  = gcdPoly d d'
  in if degree g == 0
     then rothsteinTragerExp a d field
     else Left "Repeated factors in exponential case — not yet implemented"

-- | Rothstein-Trager for the exponential case
rothsteinTragerExp :: Poly Double -> Poly Double -> DiffField -> Either String [(Double, Poly Double)]
rothsteinTragerExp a d field =
  let d'    = diffPoly d
      rPoly = resultantPoly a d d'
  in case findRationalRoots rPoly of
       Nothing    -> Left "NonElementary"
       Just roots ->
         Right [ (c, gcdPoly d (subPoly a (scalePoly c d')))
               | c <- roots
               , c /= 0
               ]

-- | Integrate exp(f) directly
integrateExp :: Expr -> String -> Either String Expr
integrateExp f var =
  let f'     = simplify (deriveBase f)
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
containsVar _ (Const _)   = False
containsVar _ Pi          = False
containsVar _ E           = False
containsVar _ I           = False
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