-- | Risch algorithm: exponential case.
--
-- Handles integration of rational functions in a single exponential
-- extension @t = exp(θ)@ where @θ' = D(θ)@ is a constant.
-- Covers @e^(ax)@, @sin(ax)@, and @cos(ax)@ (the last two via the
-- Euler substitution @t = e^(ix)@ used in
-- 'LimCalc.Integration.Risch').
--
-- = Algorithm outline
--
-- 1. /Hermite reduction/: split @f = g' + h@ where @g@ is a polynomial
--    in @t@ (the rational-integral part) and @h@ is a proper fraction.
--
-- 2. /Polynomial part/: integrate @g@ term by term using
--    @D(t^n) = n · θ' · t^n@, giving antiderivative
--    @∑ (cₙ\/(n·θ')) · t^n@.
--
-- 3. /Reduced part/: if the denominator is a pure power of @t@
--    (a Laurent term), integrate directly. Otherwise apply a
--    version of Rothstein-Trager using the correct differential-field
--    derivative @D(d) = θ' · t · d'(t)@ (not the ordinary polynomial
--    derivative).
--
-- = Complex roots
--
-- Unlike 'LimCalc.Integration.Risch.Primitive', this module uses
-- 'findComplexRoots' (not 'findRationalRoots') in
-- 'rothsteinTragerExp'. The field has already been extended to include
-- @i@ via the Euler substitution, so complex-coefficient logarithms
-- are legitimate here. The imaginary parts cancel when the Euler-form
-- output is folded back to @sin\/cos@ via
-- 'LimCalc.Core.Simplify.foldEuler'.
--
-- = Key corrections from the original implementation
--
-- * The differential-field derivative @D(d) = θ' · t · d'(t)@ is used
--   throughout, not the ordinary polynomial derivative @d'(t)@. These
--   coincide only when @θ' = 1@ and the monomial shift is ignored —
--   an accidental coincidence that masked the bug for simple cases.
--
-- * Laurent terms (denominator = pure power of @t@) are detected via
--   'monomialDegree' and integrated directly, bypassing Rothstein-Trager
--   (which degenerates on such inputs since @gcd(t^k, D(t^k)) = t^k ≠ 1@).
module LimCalc.Integration.Risch.Exponential
  ( -- * Integration
    integrateExponential
  , integrateExp
  , ExponentialResult (..)
    -- * Field analysis
  , findExponentialTheta
  , constThetaPrime
  , evalConstExpr
    -- * Term integration
  , integrateSingleExpTerm
  , integrateExpPolyPart
  , integrateExponentialReduced
    -- * Rothstein-Trager (exponential)
  , rothsteinTragerExp
  , expFieldDeriv
  , monomialDegree
    -- * Utilities
  , containsVar
  , firstVarIn
  , orElse
  ) where

import LimCalc.Algebra.Poly
import LimCalc.Algebra.RationalFunction
import LimCalc.Differentiation.DiffField
import LimCalc.Core.Expr
import LimCalc.Core.Simplify
import LimCalc.Integration.Risch.Primitive
import LimCalc.Algebra.AlgNum

-- | Result of the Risch exponential-case integration algorithm.
data ExponentialResult
  = ExponentialElementary Expr
    -- ^ The integral is elementary; the 'Expr' is the antiderivative.
  | ExponentialNonElementary
    -- ^ The integral is non-elementary.
  | ExponentialError String
    -- ^ An internal error occurred.
  deriving (Show, Eq)

-- | Integrate a rational function in the exponential case.
--
-- Requires that the 'DiffField' contains exactly one exponential
-- extension @t = exp(θ)@ with @θ' = D(θ)@ a constant. Applies
-- Hermite reduction, integrates the polynomial part via
-- 'integrateExpPolyPart', and handles the reduced part either
-- as a Laurent term or via 'integrateExponentialReduced'.
integrateExponential :: RatFun AlgNum -> DiffField -> ExponentialResult
integrateExponential rf field =
  case findExponentialTheta field of
    Nothing -> ExponentialError "No exponential extension found in field"
    Just theta ->
      case constThetaPrime theta of
        Nothing -> ExponentialError
          "Exponential extension's argument has a non-constant \
          \derivative -- only theta' = constant is supported"
        Just thetaPrime ->
          let (g, reduced@(RatFun a d)) = hermiteReduce rf
              polyPartG = integrateExpPolyPart thetaPrime (numerator g) theta
          in case monomialDegree d of
               Just k ->
                 let aConst      = headOrZero (polyCoef a)
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

-- | Find the argument @θ@ of the first exponential extension in a
-- differential field tower.
findExponentialTheta :: DiffField -> Maybe Expr
findExponentialTheta field = go (extensions field)
  where
    go []                   = Nothing
    go (Exponential f : _)  = Just f
    go (_ : rest)           = go rest

-- | Evaluate @D(θ)@ and verify it is a constant (independent of the
-- base variable), returning its value as an 'AlgNum'.
--
-- Constants may involve @i@ (e.g. @θ = i·a·x@ gives @θ' = i·a@),
-- which 'evalConstExpr' handles.
constThetaPrime :: Expr -> Maybe AlgNum
constThetaPrime theta =
  let thetaPrime = simplify (deriveBase theta)
  in evalConstExpr thetaPrime

-- | Evaluate an 'Expr' to an 'AlgNum' if it is a closed-form constant
-- (built from 'Const', 'I', 'Pi', 'E', and arithmetic on those only).
-- Returns 'Nothing' if the expression contains any 'Var'.
evalConstExpr :: Expr -> Maybe AlgNum
evalConstExpr (Const c)         = Just (fromRational (toRational c))
evalConstExpr Pi                = Just (fromRational (toRational (pi :: Double)))
evalConstExpr E                 = Just (fromRational (toRational (exp 1 :: Double)))
evalConstExpr I                 = Just algI
evalConstExpr (Var _)           = Nothing
evalConstExpr (Add f g)         = (+) <$> evalConstExpr f <*> evalConstExpr g
evalConstExpr (Sub f g)         = (-) <$> evalConstExpr f <*> evalConstExpr g
evalConstExpr (Mul f g)         = (*) <$> evalConstExpr f <*> evalConstExpr g
evalConstExpr (Div f g)         = (/) <$> evalConstExpr f <*> evalConstExpr g
evalConstExpr (Neg f)           = negate <$> evalConstExpr f
evalConstExpr (Pow f (Const n)) =
  (\base -> base ** fromRational (toRational n)) <$> evalConstExpr f
evalConstExpr _                 = Nothing

-- | Integrate a single term @c · t^n@ where @t = exp(θ)@ and @n@ is
-- any integer (positive, negative, or zero).
--
-- * @n = 0@: @∫ c dx = c · x@
-- * @n ≠ 0@: @∫ c · t^n dx = (c\/(n·θ')) · t^n@
integrateSingleExpTerm :: AlgNum -> AlgNum -> Int -> Expr -> Expr
integrateSingleExpTerm _ c 0 _
  | isAlgZero c = Const 0
integrateSingleExpTerm _ c 0 theta = Mul (algNumToExpr c) baseVarOf
  where
    baseVarOf = case firstVarIn theta of
      Just v  -> Var v
      Nothing -> Var "x"
integrateSingleExpTerm thetaPrime c n theta
  | isAlgZero c = Const 0
  | otherwise   =
      Mul (algNumToExpr (c / (fromIntegral n * thetaPrime)))
          (Pow (Exp theta) (Const (fromIntegral n)))

-- | Find the first variable name appearing in an 'Expr'.
firstVarIn :: Expr -> Maybe String
firstVarIn (Var x)   = Just x
firstVarIn (Add f g) = firstVarIn f `orElse` firstVarIn g
firstVarIn (Sub f g) = firstVarIn f `orElse` firstVarIn g
firstVarIn (Mul f g) = firstVarIn f `orElse` firstVarIn g
firstVarIn (Div f g) = firstVarIn f `orElse` firstVarIn g
firstVarIn (Pow f g) = firstVarIn f `orElse` firstVarIn g
firstVarIn (Neg f)   = firstVarIn f
firstVarIn (Exp f)   = firstVarIn f
firstVarIn (Log f)   = firstVarIn f
firstVarIn (Sin f)   = firstVarIn f
firstVarIn (Cos f)   = firstVarIn f
firstVarIn (Abs f)   = firstVarIn f
firstVarIn _         = Nothing

-- | @orElse x y@ returns @x@ if it is @Just@, otherwise @y@.
orElse :: Maybe a -> Maybe a -> Maybe a
orElse (Just x) _ = Just x
orElse Nothing  y = y

-- | Integrate the polynomial-in-@t@ part returned by Hermite reduction.
--
-- For each term @cₙ · t^n@:
--
-- * @n = 0@: @∫ c₀ dx = c₀ · x@
-- * @n ≥ 1@: @∫ cₙ · t^n dx = (cₙ\/(n·θ')) · t^n@
--
-- where @t = exp(θ)@ and @θ' = D(θ)@ is a constant.
integrateExpPolyPart :: AlgNum -> Poly AlgNum -> Expr -> Expr
integrateExpPolyPart _ (Poly _ []) _ = Const 0
integrateExpPolyPart thetaPrime (Poly x cs) theta =
  let tExpr = Exp theta
      terms  = [ termFor n c
               | (n, c) <- zip [0 :: Int ..] cs
               , not (isAlgZero c)
               ]
      termFor 0 c = Mul (algNumToExpr c) (Var x)
      termFor n c =
        Mul (algNumToExpr (c / (fromIntegral n * thetaPrime)))
            (Pow tExpr (Const (fromIntegral n)))
  in foldr addExpr (Const 0) terms

-- | Handle the reduced (proper-fraction) part in the exponential case.
--
-- If the denominator is a pure power of @t@ (detected via
-- 'monomialDegree'), returns an empty term list — the Laurent-term
-- integral is handled by the caller. Otherwise, verifies the
-- denominator is squarefree under the differential-field derivative
-- 'expFieldDeriv', then applies 'rothsteinTragerExp'.
integrateExponentialReduced :: AlgNum -> RatFun AlgNum -> DiffField
                            -> Either String [(AlgNum, Poly AlgNum)]
integrateExponentialReduced thetaPrime (RatFun a d) _field
  | Just _ <- monomialDegree d = Right []
  | otherwise =
      let dPrime = expFieldDeriv thetaPrime d
          g      = gcdPoly d dPrime
      in if degree g == 0
           then rothsteinTragerExp thetaPrime a d _field
           else Left "Repeated factors in exponential case (beyond a \
                      \bare power of the generator) -- not yet implemented"

-- | Test whether a polynomial is a pure monomial @c · t^k@ (degree @k@,
-- all lower coefficients zero). Returns @Just k@ if so, @Nothing@ otherwise.
--
-- Used to detect the Laurent-term case, where Rothstein-Trager is
-- inapplicable (since @gcd(t^k, D(t^k)) = t^k \neq 1@).
monomialDegree :: Poly AlgNum -> Maybe Int
monomialDegree (Poly _ cs) =
  let n = length cs - 1
  in if n >= 0 && all isAlgZero (take n cs)
       then Just n
       else Nothing

-- | The differential-field derivative of a polynomial @d(t)@ in the
-- exponential extension @t = exp(θ)@:
--
-- @D(∑ cₙ tⁿ) = θ' · t · d'(t)@
--
-- where @d'(t)@ is the ordinary polynomial derivative. Implemented
-- as @θ' · (0 : diffPoly(d))@ (prepending a zero coefficient to
-- shift the degree up by one, representing multiplication by @t@).
expFieldDeriv :: AlgNum -> Poly AlgNum -> Poly AlgNum
expFieldDeriv thetaPrime d =
  scalePoly thetaPrime (shiftUpByOne (diffPoly d))
  where
    shiftUpByOne (Poly x cs) = Poly x (algZero : cs)

-- | Rothstein-Trager algorithm for the exponential case, using the
-- correct differential-field derivative 'expFieldDeriv' throughout.
--
-- Uses 'findComplexRoots' (not 'findRationalRoots') since complex
-- roots are legitimate in the @i@-extended field used by trig
-- integration.
rothsteinTragerExp :: AlgNum -> Poly AlgNum -> Poly AlgNum -> DiffField
                   -> Either String [(AlgNum, Poly AlgNum)]
rothsteinTragerExp thetaPrime a d _field =
  let dPrime = expFieldDeriv thetaPrime d
      rPoly  = resultantPoly a d dPrime
  in case findComplexRoots rPoly of
       Nothing    -> Left "NonElementary"
       Just roots ->
         Right [ (c, gcdPoly d (subPoly a (scalePoly c dPrime)))
               | c <- roots
               , not (isAlgZero c)
               ]

-- | Directly integrate @exp(f)@ with respect to @var@.
--
-- Succeeds when @D(f)@ is a non-zero constant with respect to @var@,
-- giving @∫ exp(f) dx = exp(f) \/ f'@. Returns @Left@ for non-constant
-- or zero derivatives (non-elementary or degenerate cases).
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

-- | Test whether an 'Expr' contains a given variable name.
containsVar :: String -> Expr -> Bool
containsVar v (Var x)    = v == x
containsVar _ (Const _)  = False
containsVar _ Pi         = False
containsVar _ E          = False
containsVar _ I          = False
containsVar v (Add f g)  = containsVar v f || containsVar v g
containsVar v (Sub f g)  = containsVar v f || containsVar v g
containsVar v (Mul f g)  = containsVar v f || containsVar v g
containsVar v (Div f g)  = containsVar v f || containsVar v g
containsVar v (Pow f g)  = containsVar v f || containsVar v g
containsVar v (Neg f)    = containsVar v f
containsVar v (Exp f)    = containsVar v f
containsVar v (Log f)    = containsVar v f
containsVar v (Sin f)    = containsVar v f
containsVar v (Cos f)    = containsVar v f
containsVar v (Abs f)    = containsVar v f
containsVar v (Erf f)    = containsVar v f
containsVar v (Li f)     = containsVar v f
containsVar v (Si f)     = containsVar v f
containsVar v (Ci f)     = containsVar v f
containsVar v (Ei f)     = containsVar v f
containsVar v (Arcsin f) = containsVar v f
containsVar v (Arccos f) = containsVar v f
containsVar v (Arctan f) = containsVar v f