module LimCalc.Risch where

import LimCalc.Expr
import LimCalc.Poly
import LimCalc.RationalFunction
import LimCalc.DiffField
import LimCalc.Simplify
import LimCalc.Risch.Primitive
import LimCalc.Risch.Exponential

-- | Result of the Risch integration algorithm
data RischResult
  = Elementary Expr         -- ^ Elementary antiderivative found
  | NonElementary           -- ^ Provably no elementary antiderivative
  | NotImplemented String   -- ^ Case not yet implemented
  | RischError String       -- ^ Algorithm error
  deriving (Show, Eq)

-- | Top-level Risch integration
-- Analyzes the expression, builds the differential field tower,
-- and dispatches to the appropriate case
rischIntegrate :: Expr -> String -> RischResult
rischIntegrate f var =
  case classifyIntegrand f var of
    RationalCase rf ->
      let field = baseField var
      in case integratePrimitive rf field of
           PrimitiveElementary e  -> Elementary (simplify e)
           PrimitiveNonElementary -> NonElementary
           PrimitiveError s       -> RischError s

    ExponentialCase expArg ->
      case integrateExp expArg var of
        Right e                                   -> Elementary (simplify e)
        Left "Non-elementary exponential integral" -> NonElementary
        Left msg                                   -> NotImplemented msg

    LogCase logArg ->
      -- int log(f) dx = x·log(f) - int x·f'/f dx
      -- Integration by parts
      let inner    = Mul (Var var) (Div (deriveBase logArg) logArg)
          innerResult = rischIntegrate inner var
      in case innerResult of
           Elementary e ->
             Elementary (simplify (Sub (Mul (Var var) (Log logArg)) e))
           other -> other

    TrigCase ->
      -- Convert sin/cos to exp via Euler's formula and integrate
      NotImplemented "Trig integration via exponential extension"

    PolynomialCase poly ->
      -- Integrate polynomial term by term
      Elementary (simplify (integratePolynomial poly var))

    GeneralCase ->
      NotImplemented "General case — full tower not yet implemented"

-- | Classification of integrand
data IntegrandClass
  = RationalCase RatFun       -- ^ Pure rational function
  | ExponentialCase Expr      -- ^ exp(f) times rational
  | LogCase Expr              -- ^ log(f) times rational
  | TrigCase                  -- ^ trig functions
  | PolynomialCase Expr       -- ^ pure polynomial
  | GeneralCase               -- ^ general case
  deriving (Show)

-- | Classify an integrand
classifyIntegrand :: Expr -> String -> IntegrandClass
classifyIntegrand f var
  | isPoly f var    = PolynomialCase f
  | isRational f var = RationalCase (exprToRatFun f var)
  | isExpForm f var  = case extractExpArg f of
                         Just arg -> ExponentialCase arg
                         Nothing  -> GeneralCase
  | isLogForm f var  = case extractLogArg f of
                         Just arg -> LogCase arg
                         Nothing  -> GeneralCase
  | hasTrig f        = TrigCase
  | otherwise        = GeneralCase

-- | Check if expression is a polynomial in var
isPoly :: Expr -> String -> Bool
isPoly (Const _)   _   = True
isPoly (Var x)     var = x == var
isPoly (Add f g)   var = isPoly f var && isPoly g var
isPoly (Sub f g)   var = isPoly f var && isPoly g var
isPoly (Mul f g)   var = isPoly f var && isPoly g var
isPoly (Neg f)     var = isPoly f var
isPoly (Pow (Var x) (Const n)) var = x == var && n >= 0 && n == fromIntegral (round n)
isPoly _           _   = False

-- | Check if expression is rational in var
isRational :: Expr -> String -> Bool
isRational (Div f g) var = isPoly f var && isPoly g var
isRational f         var = isPoly f var

-- | Check if expression is of the form exp(f)
isExpForm :: Expr -> String -> Bool
isExpForm (Exp _)       _ = True
isExpForm (Mul f (Exp _)) _ = True
isExpForm (Mul (Exp _) f) _ = True
isExpForm _             _ = False

-- | Check if expression involves log
isLogForm :: Expr -> String -> Bool
isLogForm (Log _)         _ = True
isLogForm (Mul _ (Log _)) _ = True
isLogForm (Mul (Log _) _) _ = True
isLogForm _               _ = False

-- | Check if expression involves trig
hasTrig :: Expr -> Bool
hasTrig (Sin _)   = True
hasTrig (Cos _)   = True
hasTrig (Add f g) = hasTrig f || hasTrig g
hasTrig (Mul f g) = hasTrig f || hasTrig g
hasTrig (Div f g) = hasTrig f || hasTrig g
hasTrig (Neg f)   = hasTrig f
hasTrig _         = False

-- | Extract argument of exp
extractExpArg :: Expr -> Maybe Expr
extractExpArg (Exp f)         = Just f
extractExpArg (Mul _ (Exp f)) = Just f
extractExpArg (Mul (Exp f) _) = Just f
extractExpArg _               = Nothing

-- | Extract argument of log
extractLogArg :: Expr -> Maybe Expr
extractLogArg (Log f)         = Just f
extractLogArg (Mul _ (Log f)) = Just f
extractLogArg (Mul (Log f) _) = Just f
extractLogArg _               = Nothing

-- | Convert expression to rational function
exprToRatFun :: Expr -> String -> RatFun
exprToRatFun (Div f g) var = ratFun (exprToPoly f var) (exprToPoly g var)
exprToRatFun f         var = ratFun (exprToPoly f var) (onePoly var)

-- | Convert polynomial expression to Poly
exprToPoly :: Expr -> String -> Poly
exprToPoly (Const c)   var = constPoly var c
exprToPoly (Var x)     var
  | x == var  = Poly var [0, 1]
  | otherwise = constPoly var 0
exprToPoly (Add f g)   var = addPoly (exprToPoly f var) (exprToPoly g var)
exprToPoly (Sub f g)   var = subPoly (exprToPoly f var) (exprToPoly g var)
exprToPoly (Mul f g)   var = mulPoly (exprToPoly f var) (exprToPoly g var)
exprToPoly (Neg f)     var = negPoly (exprToPoly f var)
exprToPoly (Pow (Var x) (Const n)) var
  | x == var  = monomialPoly var 1 (round n)
  | otherwise = constPoly var 0
exprToPoly _           var = constPoly var 0

-- | Integrate a polynomial expression term by term
integratePolynomial :: Expr -> String -> Expr
integratePolynomial f var =
  let p   = exprToPoly f var
      p'  = integratePoly p
  in polyToExpr p'

-- | Integrate a polynomial term by term
-- int Σ aₙxⁿ dx = Σ aₙxⁿ⁺¹/(n+1)
integratePoly :: Poly -> Poly
integratePoly (Poly x []) = zeroPoly x
integratePoly (Poly x cs) =
  Poly x $ 0 : [ c / fromIntegral (n+1)
                | (n, c) <- zip [0..] cs ]