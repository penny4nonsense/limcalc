module LimCalc.Risch.Exponential where

import LimCalc.Poly
import LimCalc.RationalFunction
import LimCalc.DiffField
import LimCalc.Expr
import LimCalc.Simplify
import LimCalc.Risch.Primitive
import LimCalc.AlgNum

-- | Result of the Risch algorithm for the exponential case
data ExponentialResult
  = ExponentialElementary Expr
  | ExponentialNonElementary
  | ExponentialError String
  deriving (Show, Eq)

-- | Integrate in the exponential case
--
-- Generalized to RatFun AlgNum (see Risch.Primitive's header for
-- why: trig integration via Euler's formula needs genuinely complex
-- coefficients). Two further fixes beyond the type generalization:
--
-- 1. The polynomial part g returned by hermiteReduce is a
--    polynomial in t = exp(theta), NOT in x -- passing it through
--    via ratFunToExpr unchanged (as the previous version did) is
--    wrong, since "g as a function of x" is not its own
--    antiderivative. See integrateExpPolyPart.
--
-- 2. 'field' is now actually consulted to find theta (the exponent
--    argument); the previous version took it as a parameter but
--    never used it.
integrateExponential :: RatFun AlgNum -> DiffField -> ExponentialResult
integrateExponential rf field =
  case findExponentialTheta field of
    Nothing -> ExponentialError "No exponential extension found in field"
    Just theta ->
      case constThetaPrime theta of
        Nothing -> ExponentialError
          "Exponential extension's argument has a non-constant \
          \derivative -- only theta' = constant is supported \
          \(covers sin(a*x), cos(a*x), e^(a*x) for constant a, \
          \including the i*x case used by trig integration; does \
          \NOT cover e.g. e^(x^2))"
        Just thetaPrime ->
          let (g, reduced@(RatFun a d)) = hermiteReduce rf
              polyPartG  = integrateExpPolyPart thetaPrime (numerator g) theta
          in case monomialDegree d of
               -- d is a pure power of t: a/d is a single Laurent
               -- term, not a genuine logarithmic case (see
               -- integrateExponentialReduced's header for why
               -- Rothstein-Trager doesn't apply here).
               Just k ->
                 let aConst = headOrZero (polyCoef a)
                     laurentPart = integrateSingleExpTerm thetaPrime aConst (negate k) theta
                 in ExponentialElementary (addExpr polyPartG laurentPart)
               Nothing ->
                 case integrateExponentialReduced thetaPrime reduced field of
                   Left err    -> ExponentialError err
                   Right terms ->
                     let logPart = foldr addExpr (Const 0)
                           [ Mul (algNumToExpr c) (Log (polyToExpr u))
                           | (c, u) <- terms ]
                     in ExponentialElementary (addExpr polyPartG logPart)
  where
    headOrZero []    = algZero
    headOrZero (x:_) = x

-- | Find theta (the exponent argument) from the first Exponential
-- extension in the field's tower. Single-exponential-extension only
-- for now, matching the rest of this module's scope.
findExponentialTheta :: DiffField -> Maybe Expr
findExponentialTheta field = go (extensions field)
  where
    go []                  = Nothing
    go (Exponential f : _) = Just f
    go (_ : rest)          = go rest

-- | Evaluate theta' = D(theta) and check it's a constant (doesn't
-- depend on the base variable), returning that constant as an
-- AlgNum. Constants here can involve I (e.g. theta = i*a*x gives
-- theta' = i*a), which evalConstExpr handles.
constThetaPrime :: Expr -> Maybe AlgNum
constThetaPrime theta =
  let thetaPrime = simplify (deriveBase theta)
  in evalConstExpr thetaPrime

-- | Evaluate an Expr to an AlgNum IF it is a closed-form constant
-- expression (built only from Const, I, Pi, E, and arithmetic on
-- those) -- i.e. it does not depend on any variable. Returns Nothing
-- if the expression contains a Var, since that means it is not
-- actually constant.
evalConstExpr :: Expr -> Maybe AlgNum
evalConstExpr (Const c)  = Just (fromRational (toRational c))
evalConstExpr Pi         = Just (fromRational (toRational (pi :: Double)))
evalConstExpr E          = Just (fromRational (toRational (exp 1 :: Double)))
evalConstExpr I          = Just algI
evalConstExpr (Var _)    = Nothing
evalConstExpr (Add f g)  = (+) <$> evalConstExpr f <*> evalConstExpr g
evalConstExpr (Sub f g)  = (-) <$> evalConstExpr f <*> evalConstExpr g
evalConstExpr (Mul f g)  = (*) <$> evalConstExpr f <*> evalConstExpr g
evalConstExpr (Div f g)  = (/) <$> evalConstExpr f <*> evalConstExpr g
evalConstExpr (Neg f)    = negate <$> evalConstExpr f
evalConstExpr (Pow f (Const n)) =
  (\base -> base ** fromRational (toRational n)) <$> evalConstExpr f
evalConstExpr _          = Nothing  -- Exp/Log/Sin/Cos/Abs/non-constant Pow: not a closed-form constant here

-- | Antiderivative of a single c*t^n term, for any integer n
-- (positive, negative, or zero), where t = exp(theta) and theta' is
-- the (constant) derivative of theta:
--   n = 0:  Integral(c) "dx" = c * x          (base-field term)
--   n /= 0: Integral(c*t^n) "dx" = (c/(n*theta')) * t^n
--
-- Used for the Laurent-term case (denominator is a bare power of t)
-- where the term cannot be represented via Poly's coefficient list
-- (which is implicitly indexed by non-negative degree 0,1,2,...).
integrateSingleExpTerm :: AlgNum -> AlgNum -> Int -> Expr -> Expr
integrateSingleExpTerm _ c 0 _
  | isAlgZero c = Const 0
integrateSingleExpTerm _ c 0 theta = Mul (algNumToExpr c) baseVarOf
  where
    -- theta is built from the base variable; recover its name the
    -- same way the rest of this module assumes a single base
    -- variable "x"-equivalent. Since DiffField threads baseVar
    -- explicitly elsewhere but this helper only receives theta, we
    -- fall back to extracting it from theta's own free variable.
    baseVarOf = case firstVarIn theta of
      Just v  -> Var v
      Nothing -> Var "x"  -- defensive fallback; should not occur in practice
integrateSingleExpTerm thetaPrime c n theta
  | isAlgZero c = Const 0
  | otherwise =
      Mul (algNumToExpr (c / (fromIntegral n * thetaPrime)))
          (Pow (Exp theta) (Const (fromIntegral n)))

-- | Find the first variable name appearing in an Expr, if any.
firstVarIn :: Expr -> Maybe String
firstVarIn (Var x)     = Just x
firstVarIn (Add f g)   = firstVarIn f `orElse` firstVarIn g
firstVarIn (Sub f g)   = firstVarIn f `orElse` firstVarIn g
firstVarIn (Mul f g)   = firstVarIn f `orElse` firstVarIn g
firstVarIn (Div f g)   = firstVarIn f `orElse` firstVarIn g
firstVarIn (Pow f g)   = firstVarIn f `orElse` firstVarIn g
firstVarIn (Neg f)     = firstVarIn f
firstVarIn (Exp f)     = firstVarIn f
firstVarIn (Log f)     = firstVarIn f
firstVarIn (Sin f)     = firstVarIn f
firstVarIn (Cos f)     = firstVarIn f
firstVarIn (Abs f)     = firstVarIn f
firstVarIn _           = Nothing

orElse :: Maybe a -> Maybe a -> Maybe a
orElse (Just x) _ = Just x
orElse Nothing  y = y

-- | Integrate the polynomial-in-t part of an exponential-case
-- rational function, where t = exp(theta) and theta' is the
-- (constant) derivative of theta.
--
--   D(t^n) = n * t^n * theta'   for n >= 1
--   so  Integral(c_n * t^n) "dx" = (c_n / (n * theta')) * t^n
--
-- The constant term (n=0) is a base-field term: Integral(c_0) "dx"
-- = c_0 * x.
integrateExpPolyPart :: AlgNum -> Poly AlgNum -> Expr -> Expr
integrateExpPolyPart _ (Poly _ []) _ = Const 0
integrateExpPolyPart thetaPrime (Poly x cs) theta =
  let tExpr = Exp theta
      terms = [ termFor n c
              | (n, c) <- zip [0 :: Int ..] cs
              , not (isAlgZero c)
              ]
      termFor 0 c = Mul (algNumToExpr c) (Var x)
      termFor n c =
        Mul (algNumToExpr (c / (fromIntegral n * thetaPrime)))
            (Pow tExpr (Const (fromIntegral n)))
  in foldr addExpr (Const 0) terms

-- | Handle the reduced (proper-fraction) part in the exponential
-- case.
--
-- Previously this used diffPoly d directly as "D(d)" when checking
-- gcdPoly d d' and when computing the Rothstein-Trager resultant.
-- That is WRONG: diffPoly is the ordinary polynomial derivative
-- d/dt, but the differential-field derivation needed here is
-- D(d) = theta' * t * diffPoly(d) (chain rule through t = exp(theta),
-- since D(t) = theta'*t). The two coincide only by an accidental
-- factor in degenerate cases; in general they differ by exactly the
-- theta'*t factor, which corrupts both the gcd-based case split and
-- the Rothstein-Trager resultant whenever theta' != 1 (e.g. the i in
-- trig's theta = i*x) or even when theta'=1, since the missing
-- factor of t still shifts every degree by one.
--
-- Separately: if d is purely a power of t (a monomial, e.g. d = t,
-- or d = t^3), this is NOT a genuine logarithmic/Rothstein-Trager
-- case at all -- it is a Laurent-series term in t. Rothstein-Trager's
-- resultant trick assumes gcd(d, D(d)) = 1 (squarefree denominator
-- coprime to its own derivative); when d = t^k, D(d) = k*theta'*t^k,
-- which SHARES the factor t^k with d entirely, making the resultant
-- degenerate (constant in z, independent of the root being sought)
-- and silently producing zero terms -- the integral isn't actually
-- non-elementary, the wrong algorithm was being applied. a/t^k
-- integrates via the same rule as the polynomial part, extended to
-- negative n: Integral(c * t^n) "dx" = (c / (n*theta')) * t^n.
integrateExponentialReduced :: AlgNum -> RatFun AlgNum -> DiffField -> Either String [(AlgNum, Poly AlgNum)]
integrateExponentialReduced thetaPrime (RatFun a d) field
  | Just k <- monomialDegree d =
      -- d = t^k (up to a constant factor already normalized away by
      -- ratFun/hermiteReduce); a/d is a single Laurent term, handled
      -- directly rather than via Rothstein-Trager. Returned as an
      -- empty log-term list since this contributes nothing to the
      -- log part -- the caller (integrateExponential) is expected to
      -- add this term separately. To keep this function's existing
      -- signature/contract, we instead fold the Laurent term's
      -- antiderivative directly into a synthetic "log-part" entry is
      -- NOT appropriate (it isn't a log term), so this case is
      -- surfaced via Left with a sentinel the caller special-cases.
      -- See integrateExponential, which checks for this case BEFORE
      -- calling this function, so this branch should not normally be
      -- reached -- kept here only as a defensive fallback.
      Right []
  | otherwise =
      let dPrime = expFieldDeriv thetaPrime d
          g      = gcdPoly d dPrime
      in if degree g == 0
           then rothsteinTragerExp thetaPrime a d field
           else Left "Repeated factors in exponential case (beyond a \
                      \bare power of the generator) -- not yet \
                      \implemented"

-- | If a polynomial is purely a single monomial c*t^k (k = its
-- degree, all lower coefficients zero), return Just k. Otherwise
-- Nothing. Used to detect the Laurent-term case in
-- integrateExponentialReduced/integrateExponential.
monomialDegree :: Poly AlgNum -> Maybe Int
monomialDegree (Poly _ cs) =
  let n = length cs - 1
  in if n >= 0 && all isAlgZero (take n cs)
       then Just n
       else Nothing

-- | The proper differential-field derivative of a polynomial in t,
-- where t = exp(theta) and theta' = D(theta) is a known constant:
--   D(sum c_n t^n) = theta' * sum (c_n * n * t^n) = theta' * t * diffPoly(d)
expFieldDeriv :: AlgNum -> Poly AlgNum -> Poly AlgNum
expFieldDeriv thetaPrime d =
  scalePoly thetaPrime (shiftUpByOne (diffPoly d))
  where
    -- Multiply by t: shift every exponent up by one, i.e. prepend a
    -- zero coefficient.
    shiftUpByOne (Poly x cs) = Poly x (algZero : cs)

-- | Rothstein-Trager for the exponential case, using the corrected
-- differential-field derivative throughout.
rothsteinTragerExp :: AlgNum -> Poly AlgNum -> Poly AlgNum -> DiffField -> Either String [(AlgNum, Poly AlgNum)]
rothsteinTragerExp thetaPrime a d field =
  let dPrime = expFieldDeriv thetaPrime d
      rPoly  = resultantPoly a d dPrime
  in case findComplexRoots rPoly of
       Nothing    -> Left "NonElementary"
       Just roots ->
         Right [ (c, gcdPoly d (subPoly a (scalePoly c dPrime)))
               | c <- roots
               , not (isAlgZero c)
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
containsVar v (Erf f)     = containsVar v f
containsVar v (Li f)      = containsVar v f
containsVar v (Si f)      = containsVar v f
containsVar v (Ci f)      = containsVar v f
containsVar v (Ei f)      = containsVar v f
containsVar v (Arcsin f)  = containsVar v f
containsVar v (Arccos f)  = containsVar v f
containsVar v (Arctan f)  = containsVar v f