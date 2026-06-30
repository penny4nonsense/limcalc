module LimCalc.Risch where

import LimCalc.Expr
import LimCalc.Poly
import LimCalc.RationalFunction
import LimCalc.DiffField
import LimCalc.Simplify
import LimCalc.Risch.Primitive
import LimCalc.Risch.Exponential
import LimCalc.AlgNum

-- | Result of the Risch integration algorithm
data RischResult
  = Elementary Expr
  | NonElementary
  | NotImplemented String
  | RischError String
  deriving (Show, Eq)

-- | Top-level Risch integration
rischIntegrate :: Expr -> String -> RischResult
rischIntegrate f var =
  case recognizeSpecialIntegral f var of
    Just result -> result
    Nothing -> rischIntegrateClassified f var

-- | Recognize the five classical integrals whose antiderivatives are
-- the standard non-elementary special functions, returning the
-- closed form directly rather than falling through to NonElementary.
--
-- These don't fit cleanly into the existing classification scheme:
--   - 1/log(x) doesn't match isLogForm (which only recognizes
--     Log f, g*log(f), or log(f)*g -- not log(f) as a denominator).
--   - sin(x)/x, cos(x)/x, e^x/x are classified as TrigCase or
--     ExponentialCase via hasTrig/isExpForm's Div clauses, but
--     exprToTrigRatFun and the exponential-case machinery only
--     handle sin/cos/exp combined arithmetically with EACH OTHER,
--     not divided by the bare integration variable itself.
--
-- Scoped specifically to the argument being the bare integration
-- variable (e.g. sin(x)/x, not sin(2x)/x) -- the classical special
-- functions are defined for f(x)/x with the SAME variable in the
-- numerator's argument and the denominator, and generalizing beyond
-- that (e.g. via substitution) is a separate piece of work.
recognizeSpecialIntegral :: Expr -> String -> Maybe RischResult
recognizeSpecialIntegral f var =
  case f of
    -- int e^(-x^2) dx = (sqrt(pi)/2) * erf(x)
    Exp (Neg (Pow (Var v) (Const 2))) | v == var ->
      Just $ Elementary $ simplify $
        Mul (Div (Pow Pi (Const 0.5)) (Const 2)) (Erf (Var var))

    -- int 1/log(x) dx = li(x)
    Div (Const 1) (Log (Var v)) | v == var ->
      Just $ Elementary (Li (Var var))

    -- int sin(x)/x dx = Si(x)
    Div (Sin (Var v)) (Var v') | v == var && v' == var ->
      Just $ Elementary (Si (Var var))

    -- int cos(x)/x dx = Ci(x)
    Div (Cos (Var v)) (Var v') | v == var && v' == var ->
      Just $ Elementary (Ci (Var var))

    -- int e^x/x dx = Ei(x)
    Div (Exp (Var v)) (Var v') | v == var && v' == var ->
      Just $ Elementary (Ei (Var var))

    _ -> Nothing

rischIntegrateClassified :: Expr -> String -> RischResult
rischIntegrateClassified f var =
  case classifyIntegrand f var of
    RationalCase rf ->
      let field = baseField var
      in case integratePrimitive (ratFunDoubleToAlgNum rf) field of
           PrimitiveElementary e  -> Elementary (simplify e)
           PrimitiveNonElementary -> NonElementary
           PrimitiveError s       -> RischError s

    ExponentialCase expArg ->
      case integrateExp expArg var of
        Right e                                    -> Elementary (simplify e)
        Left "Non-elementary exponential integral" -> NonElementary
        Left msg                                   -> NotImplemented msg

    LogCase g h ->
      -- Integration by parts: int g*log(h) dx = G*log(h) - int G*(h'/h) dx
      -- where G = int g dx. Plain log(h) is the g=1 case.
      --
      -- Previously this discarded g entirely (always computed
      -- int log(h) dx regardless of what g was) and never simplified
      -- the recursive 'inner' expression before reclassifying it --
      -- which broke even the plain log(x) case, since x*(1/x)
      -- doesn't structurally look like a polynomial to
      -- classifyIntegrand without simplification collapsing it to 1
      -- first.
      case rischIntegrate (simplify g) var of
        Elementary gAntideriv ->
          let inner = simplify (Mul gAntideriv (Div (deriveBase h) h))
          in case rischIntegrate inner var of
               Elementary e ->
                 Elementary (simplify (Sub (Mul gAntideriv (Log h)) e))
               other -> other
        other -> other

    TrigCase ->
      case exprToTrigRatFun f var of
        Nothing -> NotImplemented
          "Trig integration only supports rational expressions in \
          \sin(var) and cos(var) directly (not e.g. sin(2*x) or \
          \sin(x)^2 written via Pow, or trig of a different variable)"
        Just rf ->
          let theta = Mul I (Var var)
              field = addExtension (baseField var) (Exponential theta)
          in case integrateExponential rf field of
               ExponentialElementary e   -> Elementary (foldEuler (simplify e))
               ExponentialNonElementary  -> NonElementary
               ExponentialError msg      -> NotImplemented msg

    PolynomialCase poly ->
      Elementary (simplify (integratePolynomial poly var))

    GeneralCase ->
      NotImplemented "General case — full tower not yet implemented"

-- | Convert a Sin/Cos/Add/Mul/Div/Neg/Const expression (exactly the
-- grammar hasTrig recognizes) into a rational function in
-- t = exp(i*var), via Euler's formula:
--   sin(var) = (t^2 - 1) / (2*i*t)
--   cos(var) = (t^2 + 1) / (2*t)
exprToTrigRatFun :: Expr -> String -> Maybe (RatFun AlgNum)
exprToTrigRatFun expr var = go expr
  where
    tVar = "t"
    tSquaredMinusOne = Poly tVar [negate algOne, algZero, algOne]
    tSquaredPlusOne  = Poly tVar [algOne, algZero, algOne]
    twoT  = Poly tVar [algZero, fromQ 2]
    twoIT = Poly tVar [algZero, fromQ 2 * algI]
    onePolyT = Poly tVar [algOne]

    go (Sin (Var v)) | v == var = Just (ratFun tSquaredMinusOne twoIT)
    go (Cos (Var v)) | v == var = Just (ratFun tSquaredPlusOne twoT)
    go (Const c) = Just (ratFun (Poly tVar [fromRational (toRational c)]) onePolyT)
    go (Add f g) = addRat <$> go f <*> go g
    go (Sub f g) = subRat <$> go f <*> go g
    go (Mul f g) = mulRat <$> go f <*> go g
    go (Div f g) = divRat <$> go f <*> go g
    go (Neg f)   = negRat <$> go f
    go _         = Nothing

-- | Classification of integrand
data IntegrandClass
  = RationalCase (RatFun Double)
  | ExponentialCase Expr
  | LogCase Expr Expr  -- ^ multiplier g, log argument h: integrand is g * log(h)
  | TrigCase
  | PolynomialCase Expr
  | GeneralCase
  deriving (Show)

-- | Classify an integrand
classifyIntegrand :: Expr -> String -> IntegrandClass
classifyIntegrand f var
  | isPoly f var     = PolynomialCase f
  | isRational f var = RationalCase (exprToRatFun f var)
  | isExpForm f var  = case extractExpArg f of
                         Just arg -> ExponentialCase arg
                         Nothing  -> GeneralCase
  | isLogForm f var  = case extractLogParts f of
                         Just (g, h) -> LogCase g h
                         Nothing     -> GeneralCase
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
isPoly (Pow (Var x) (Const n)) var =
  x == var && n >= 0 && n == fromIntegral (round n :: Int)
isPoly _           _   = False

-- | Check if expression is rational in var
isRational :: Expr -> String -> Bool
isRational (Div f g) var = isPoly f var && isPoly g var
isRational f         var = isPoly f var

-- | Check if expression is of the form exp(f)
isExpForm :: Expr -> String -> Bool
isExpForm (Exp _)         _ = True
isExpForm (Mul _ (Exp _)) _ = True
isExpForm (Mul (Exp _) _) _ = True
isExpForm _               _ = False

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

-- | Extract (multiplier, log argument) from a log-form expression.
-- Plain Log f is treated as multiplier 1.
extractLogParts :: Expr -> Maybe (Expr, Expr)
extractLogParts (Log f)         = Just (Const 1, f)
extractLogParts (Mul g (Log f)) = Just (g, f)
extractLogParts (Mul (Log f) g) = Just (g, f)
extractLogParts _               = Nothing

-- | Convert expression to rational function
exprToRatFun :: Expr -> String -> RatFun Double
exprToRatFun (Div f g) var = ratFun (exprToPoly f var) (exprToPoly g var)
exprToRatFun f         var = ratFun (exprToPoly f var) (onePoly var)

-- | Convert polynomial expression to Poly Double
exprToPoly :: Expr -> String -> Poly Double
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
  polyToExpr (polyDoubleToAlgNum (integratePoly (exprToPoly f var)))

-- | Integrate a polynomial term by term
integratePoly :: Poly Double -> Poly Double
integratePoly (Poly x []) = zeroPoly x
integratePoly (Poly x cs) =
  Poly x $ 0 : [ c / fromIntegral (n+1 :: Int)
                | (n, c) <- zip [0..] cs ]

-- | Convert a Poly Double to a Poly AlgNum, at the boundary between
-- this module's real-coefficient classification/polynomial-
-- integration logic and Risch.Primitive's AlgNum-generalized
-- machinery (needed for trig integration's genuinely complex
-- coefficients; ordinary polynomial/rational-function integration
-- in x never needs complex numbers, so this module's own logic
-- correctly stays in Double throughout and only converts at the
-- call boundary).
polyDoubleToAlgNum :: Poly Double -> Poly AlgNum
polyDoubleToAlgNum (Poly x cs) = Poly x (map (fromRational . toRational) cs)

-- | Convert a RatFun Double to a RatFun AlgNum, same rationale as
-- polyDoubleToAlgNum.
ratFunDoubleToAlgNum :: RatFun Double -> RatFun AlgNum
ratFunDoubleToAlgNum (RatFun p q) = RatFun (polyDoubleToAlgNum p) (polyDoubleToAlgNum q)