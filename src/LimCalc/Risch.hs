-- | Top-level Risch integration algorithm.
--
-- This module is the entry point for symbolic integration. It
-- classifies the integrand, dispatches to the appropriate sub-algorithm,
-- and applies post-processing (Euler folding, log simplification) to
-- produce human-readable output.
--
-- = Classification
--
-- Integrands are classified by 'classifyIntegrand' into:
--
-- * 'PolynomialCase': integrate term by term.
-- * 'RationalCase': Hermite reduction + Rothstein-Trager
--   ('LimCalc.Risch.Primitive.integratePrimitive').
-- * 'ExponentialCase': direct integration of @exp(f)@ when @f' = const@
--   ('LimCalc.Risch.Exponential.integrateExp').
-- * 'LogCase': integration by parts — @∫ g·log(h) dx = G·log(h) − ∫ G·(h'\/h) dx@
--   where @G = ∫ g dx@.
-- * 'TrigCase': Euler substitution @t = e^(ix)@, converting @sin\/cos@
--   to a rational function in @t@, integrated via the exponential-case
--   machinery.
-- * 'GeneralCase': not yet implemented.
--
-- = Special function recognition
--
-- Before classification, 'recognizeSpecialIntegral' pattern-matches
-- against a fixed set of integrands whose antiderivatives are standard
-- non-elementary special functions (@erf@, @li@, @Si@, @Ci@, @Ei@)
-- or algebraic functions (@arcsin@, @arctan@, @arcsinh@, @arccosh@,
-- and related forms). These don't fit cleanly into the classification
-- scheme and are handled by direct pattern matching.
--
-- = Trig integration
--
-- @sin(x)@ and @cos(x)@ are integrated via the Euler substitution
-- @t = e^(ix)@, which converts the integrand to a rational function
-- in @t@. The exponential-case Risch machinery ('LimCalc.Risch.Exponential')
-- then integrates this rational function. The result is converted back
-- to @sin\/cos@ form via 'LimCalc.Simplify.foldEuler'.
module LimCalc.Risch
  ( -- * Top-level integration
    rischIntegrate
  , RischResult (..)
    -- * Classification
  , classifyIntegrand
  , IntegrandClass (..)
    -- * Special function recognition
  , recognizeSpecialIntegral
    -- * Trig conversion
  , exprToTrigRatFun
    -- * Integrand predicates
  , isPoly
  , isRational
  , isExpForm
  , isLogForm
  , hasTrig
    -- * Argument extraction
  , extractExpArg
  , extractLogParts
    -- * Expression conversion
  , exprToRatFun
  , exprToPoly
  , integratePolynomial
  , integratePoly
  , polyDoubleToAlgNum
  , ratFunDoubleToAlgNum
  ) where

import LimCalc.Expr
import LimCalc.Poly
import LimCalc.RationalFunction
import LimCalc.DiffField
import LimCalc.Simplify
import LimCalc.Risch.Primitive
import LimCalc.Risch.Exponential
import LimCalc.AlgNum

-- | The result of Risch integration.
data RischResult
  = Elementary Expr
    -- ^ The integral is elementary; the 'Expr' is the antiderivative
    -- (without a constant of integration).
  | NonElementary
    -- ^ The integral provably has no elementary antiderivative.
  | NotImplemented String
    -- ^ The integrand falls into a case not yet implemented.
    -- The 'String' gives a reason.
  | RischError String
    -- ^ An internal error occurred during integration.
  deriving (Show, Eq)

-- | Top-level symbolic integration of @f@ with respect to @var@.
--
-- First checks 'recognizeSpecialIntegral' for direct pattern matches,
-- then dispatches via 'classifyIntegrand'.
rischIntegrate :: Expr -> String -> RischResult
rischIntegrate f var =
  case recognizeSpecialIntegral f var of
    Just result -> result
    Nothing     -> rischIntegrateClassified f var

-- | Pattern-match against integrands whose antiderivatives are
-- standard special functions or algebraic expressions.
--
-- Covers:
--
-- * @∫ e^(−x²) dx = (√π\/2) · erf(x)@
-- * @∫ 1\/log(x) dx = li(x)@
-- * @∫ sin(x)\/x dx = Si(x)@
-- * @∫ cos(x)\/x dx = Ci(x)@
-- * @∫ e^x\/x dx = Ei(x)@
-- * @∫ 1\/√(1−x²) dx = arcsin(x)@
-- * @∫ 1\/(1+x²) dx = arctan(x)@
-- * @∫ 1\/√(1+x²) dx = log(x + √(x²+1))@ (arcsinh)
-- * @∫ 1\/√(x²−1) dx = log(x + √(x²−1))@ (arccosh)
-- * @∫ x\/√(1−x²) dx = −√(1−x²)@
-- * @∫ √(1−x²) dx = x\/2 · √(1−x²) + arcsin(x)\/2@
--
-- All patterns are scoped to the exact forms shown; generalisations
-- (e.g. @sin(2x)\/x@) are not recognised here.
recognizeSpecialIntegral :: Expr -> String -> Maybe RischResult
recognizeSpecialIntegral f var =
  case f of
    Exp (Neg (Pow (Var v) (Const 2))) | v == var ->
      Just $ Elementary $ simplify $
        Mul (Div (Pow Pi (Const 0.5)) (Const 2)) (Erf (Var var))

    Div (Const 1) (Log (Var v)) | v == var ->
      Just $ Elementary (Li (Var var))

    Div (Sin (Var v)) (Var v') | v == var && v' == var ->
      Just $ Elementary (Si (Var var))

    Div (Cos (Var v)) (Var v') | v == var && v' == var ->
      Just $ Elementary (Ci (Var var))

    Div (Exp (Var v)) (Var v') | v == var && v' == var ->
      Just $ Elementary (Ei (Var var))

    Div (Const 1) (Pow (Sub (Const 1) (Pow (Var v) (Const 2))) (Const 0.5))
      | v == var ->
      Just $ Elementary (Arcsin (Var var))

    Div (Const 1) (Pow (Add (Const 1) (Pow (Var v) (Const 2))) (Const 0.5))
      | v == var ->
      Just $ Elementary $ simplify $
        Log (Add (Var var)
                 (Pow (Add (Pow (Var var) (Const 2)) (Const 1)) (Const 0.5)))

    Div (Const 1) (Pow (Sub (Pow (Var v) (Const 2)) (Const 1)) (Const 0.5))
      | v == var ->
      Just $ Elementary $ simplify $
        Log (Add (Var var)
                 (Pow (Sub (Pow (Var var) (Const 2)) (Const 1)) (Const 0.5)))

    Div (Var v) (Pow (Sub (Const 1) (Pow (Var v') (Const 2))) (Const 0.5))
      | v == var && v' == var ->
      Just $ Elementary $ simplify $
        Neg (Pow (Sub (Const 1) (Pow (Var var) (Const 2))) (Const 0.5))

    Div (Const 1) (Add (Const 1) (Pow (Var v) (Const 2)))
      | v == var ->
      Just $ Elementary (Arctan (Var var))

    Pow (Sub (Const 1) (Pow (Var v) (Const 2))) (Const 0.5)
      | v == var ->
      Just $ Elementary $ simplify $
        Add (Mul (Div (Var var) (Const 2))
                 (Pow (Sub (Const 1) (Pow (Var var) (Const 2))) (Const 0.5)))
            (Div (Arcsin (Var var)) (Const 2))

    _ -> Nothing

-- | Dispatch integration by integrand class.
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
      -- Integration by parts: ∫ g·log(h) dx = G·log(h) − ∫ G·(h'\/h) dx
      -- where G = ∫ g dx.
      case rischIntegrate (simplify g) var of
        Elementary gAntideriv ->
          let inner = simplify (Mul gAntideriv (Div (deriveBase h) h))
          in case rischIntegrate inner var of
               Elementary e -> Elementary (simplify (Sub (Mul gAntideriv (Log h)) e))
               other        -> other
        other -> other

    TrigCase ->
      case exprToTrigRatFun f var of
        Nothing -> NotImplemented
          "Trig integration only supports rational expressions in \
          \sin(var) and cos(var) directly"
        Just rf ->
          let theta = Mul I (Var var)
              field = addExtension (baseField var) (Exponential theta)
          in case integrateExponential rf field of
               ExponentialElementary e  -> Elementary (foldEuler (simplify e))
               ExponentialNonElementary -> NonElementary
               ExponentialError msg     -> NotImplemented msg

    PolynomialCase poly ->
      Elementary (simplify (integratePolynomial poly var))

    GeneralCase ->
      NotImplemented "General case — full tower not yet implemented"

-- | Convert a trig expression to a rational function in @t = e^(ix)@
-- via Euler's formula:
--
-- * @sin(x) = (t² − 1) \/ (2it)@
-- * @cos(x) = (t² + 1) \/ (2t)@
--
-- Returns 'Nothing' for any subexpression not in the grammar
-- @{sin(var), cos(var), Const, Add, Sub, Mul, Div, Neg}@.
exprToTrigRatFun :: Expr -> String -> Maybe (RatFun AlgNum)
exprToTrigRatFun expr var = go expr
  where
    tVar              = "t"
    tSquaredMinusOne  = Poly tVar [negate algOne, algZero, algOne]
    tSquaredPlusOne   = Poly tVar [algOne, algZero, algOne]
    twoT              = Poly tVar [algZero, fromQ 2]
    twoIT             = Poly tVar [algZero, fromQ 2 * algI]
    onePolyT          = Poly tVar [algOne]

    go (Sin (Var v)) | v == var = Just (ratFun tSquaredMinusOne twoIT)
    go (Cos (Var v)) | v == var = Just (ratFun tSquaredPlusOne twoT)
    go (Const c)                = Just (ratFun (Poly tVar [fromRational (toRational c)]) onePolyT)
    go (Add f g)                = addRat <$> go f <*> go g
    go (Sub f g)                = subRat <$> go f <*> go g
    go (Mul f g)                = mulRat <$> go f <*> go g
    go (Div f g)                = divRat <$> go f <*> go g
    go (Neg f)                  = negRat <$> go f
    go _                        = Nothing

-- | Classification of an integrand for dispatch.
data IntegrandClass
  = RationalCase (RatFun Double)
    -- ^ A rational function @p(x)\/q(x)@.
  | ExponentialCase Expr
    -- ^ @exp(f)@ where @f' = const@.
  | LogCase Expr Expr
    -- ^ @g · log(h)@: multiplier @g@ and log argument @h@.
    -- The @g = 1@ case covers plain @log(h)@.
  | TrigCase
    -- ^ An expression involving @sin@ and\/or @cos@.
  | PolynomialCase Expr
    -- ^ A polynomial in the integration variable.
  | GeneralCase
    -- ^ Does not fit any recognised pattern.
  deriving (Show)

-- | Classify an integrand by structural pattern matching.
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

-- | Test whether an 'Expr' is a polynomial in @var@.
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

-- | Test whether an 'Expr' is a rational function in @var@.
isRational :: Expr -> String -> Bool
isRational (Div f g) var = isPoly f var && isPoly g var
isRational f         var = isPoly f var

-- | Test whether an 'Expr' is of the form @exp(f)@ or @c · exp(f)@.
isExpForm :: Expr -> String -> Bool
isExpForm (Exp _)         _ = True
isExpForm (Mul _ (Exp _)) _ = True
isExpForm (Mul (Exp _) _) _ = True
isExpForm _               _ = False

-- | Test whether an 'Expr' involves @log@.
isLogForm :: Expr -> String -> Bool
isLogForm (Log _)         _ = True
isLogForm (Mul _ (Log _)) _ = True
isLogForm (Mul (Log _) _) _ = True
isLogForm _               _ = False

-- | Test whether an 'Expr' involves @sin@ or @cos@.
hasTrig :: Expr -> Bool
hasTrig (Sin _)   = True
hasTrig (Cos _)   = True
hasTrig (Add f g) = hasTrig f || hasTrig g
hasTrig (Mul f g) = hasTrig f || hasTrig g
hasTrig (Div f g) = hasTrig f || hasTrig g
hasTrig (Neg f)   = hasTrig f
hasTrig _         = False

-- | Extract the argument of an exponential subexpression.
extractExpArg :: Expr -> Maybe Expr
extractExpArg (Exp f)         = Just f
extractExpArg (Mul _ (Exp f)) = Just f
extractExpArg (Mul (Exp f) _) = Just f
extractExpArg _               = Nothing

-- | Extract @(multiplier, log-argument)@ from a log-form expression.
-- @log(f)@ is treated as multiplier @Const 1@.
extractLogParts :: Expr -> Maybe (Expr, Expr)
extractLogParts (Log f)         = Just (Const 1, f)
extractLogParts (Mul g (Log f)) = Just (g, f)
extractLogParts (Mul (Log f) g) = Just (g, f)
extractLogParts _               = Nothing

-- | Convert a rational-function 'Expr' to a 'RatFun Double'.
exprToRatFun :: Expr -> String -> RatFun Double
exprToRatFun (Div f g) var = ratFun (exprToPoly f var) (exprToPoly g var)
exprToRatFun f         var = ratFun (exprToPoly f var) (onePoly var)

-- | Convert a polynomial 'Expr' to a 'Poly Double'.
-- Non-polynomial subexpressions are treated as zero.
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

-- | Integrate a polynomial expression term by term.
integratePolynomial :: Expr -> String -> Expr
integratePolynomial f var =
  polyToExpr (polyDoubleToAlgNum (integratePoly (exprToPoly f var)))

-- | Integrate a 'Poly Double' term by term:
-- @∫ ∑ cₙ xⁿ dx = ∑ cₙ\/(n+1) x^(n+1)@.
integratePoly :: Poly Double -> Poly Double
integratePoly (Poly x []) = zeroPoly x
integratePoly (Poly x cs) =
  Poly x $ 0 : [ c / fromIntegral (n+1 :: Int)
                | (n, c) <- zip [0..] cs ]

-- | Lift a 'Poly Double' to 'Poly AlgNum'.
--
-- Ordinary polynomial and rational-function integration stays in
-- @Double@ throughout the classification and polynomial-integration
-- logic, converting to 'AlgNum' only at the boundary with
-- 'LimCalc.Risch.Primitive', which requires 'AlgNum' coefficients
-- to support complex numbers for trig integration.
polyDoubleToAlgNum :: Poly Double -> Poly AlgNum
polyDoubleToAlgNum (Poly x cs) = Poly x (map (fromRational . toRational) cs)

-- | Lift a 'RatFun Double' to 'RatFun AlgNum'. See 'polyDoubleToAlgNum'.
ratFunDoubleToAlgNum :: RatFun Double -> RatFun AlgNum
ratFunDoubleToAlgNum (RatFun p q) =
  RatFun (polyDoubleToAlgNum p) (polyDoubleToAlgNum q)