-- | Differential field tower for the Risch integration algorithm.
--
-- The Risch algorithm operates in a tower of differential field
-- extensions over the base field ℂ(x). Each extension adds a new
-- element θ with a known derivation Dθ expressed in terms of
-- previously introduced elements:
--
-- * /Primitive extension/ (@θ = log(f)@): @Dθ = Df\/f@
-- * /Exponential extension/ (@θ = exp(f)@): @Dθ = θ · Df@
-- * /Algebraic extension/ (@p(θ) = 0@): @Dθ = −(∂p\/∂x) \/ (∂p\/∂θ)@
--   via implicit differentiation
-- * /Special function extension/ (e.g. @θ = erf(f)@): @Dθ@ is
--   provided explicitly
--
-- = Two differentiation paths
--
-- limcalc provides two symbolic differentiation paths:
--
-- * 'deriveExpr' (this module): the /algebraic chain-rule path/.
--   Differentiates by structural recursion, treating each constructor
--   as a known primitive. Correct for iterated partial derivatives.
--   Used by 'LimCalc.Calculus.partialDiff'.
--
-- * @diff@ ('LimCalc.Calculus'): the /series extraction path/.
--   Expands @f(x + h)@ as a log-Puiseux series and reads off the
--   @h^1@ coefficient. Used by 'LimCalc.Calculus.diff'.
--
-- The series path generates Taylor coefficients by iterating
-- 'deriveBase', making it a consequence of the algebraic path rather
-- than an independent implementation.
module LimCalc.DiffField
  ( -- * Extensions
    Extension (..)
    -- * Differential field tower
  , DiffField (..)
  , baseField
  , addExtension
    -- * Derivation
  , derive
  , deriveExpr
  , deriveBase
    -- * Tower construction
  , buildTower
  , collectExts
    -- * Helpers
  , lookupExtDeriv
  , extensionDeriv
  , isBaseField
  , specialExt
  , lookupSpecial
  , knownSpecial
  ) where

import LimCalc.Expr
import LimCalc.Poly
import LimCalc.RationalFunction
import LimCalc.Simplify

-- | A single extension of the differential field, adding a new
-- element θ with a known derivation.
data Extension
  = Primitive Expr
    -- ^ Primitive (logarithmic) extension: @θ = log(f)@, @Dθ = Df\/f@.
  | Exponential Expr
    -- ^ Exponential extension: @θ = exp(f)@, @Dθ = θ · Df@.
  | Algebraic Expr
    -- ^ Algebraic extension: @p(θ, x) = 0@. The derivation is
    -- computed via implicit differentiation:
    -- @Dθ = −(∂p\/∂x) \/ (∂p\/∂θ)@.
    -- The 'Expr' is the minimal polynomial @p@ in which the extension
    -- variable appears as @Var \"t1\"@ (or @\"t2\"@, etc.) and the
    -- base variable as @Var \"x\"@.
  | Special
    -- ^ A special function extension with an explicitly provided
    -- derivative. Used for @erf@, @li@, @Si@, @Ci@, @Ei@.
      { specialName  :: String
        -- ^ Name of the special function (e.g. @\"erf\"@, @\"Si\"@).
      , specialArg   :: Expr
        -- ^ The argument of the special function.
      , specialDeriv :: Expr
        -- ^ The derivative @D(f(arg))@ as an 'Expr'.
      }
  deriving (Show, Eq)

-- | A differential field tower.
--
-- Represents the field @ℂ(x)(θ₁, θ₂, …, θₙ)@ where each @θᵢ@ is
-- an extension element with a known derivation. The base field is
-- @ℂ(x)@.
data DiffField = DiffField
  { baseVar    :: String
    -- ^ The base variable (typically @\"x\"@).
  , extensions :: [Extension]
    -- ^ Extensions in the tower, innermost first.
  , extNames   :: [String]
    -- ^ Variable names for the extensions (@\"t1\"@, @\"t2\"@, …).
  } deriving (Show, Eq)

-- | The base differential field @ℂ(x)@, with no extensions.
baseField :: String -> DiffField
baseField x = DiffField x [] []

-- | Add an extension to the top of the differential field tower.
--
-- The new extension is assigned the next available name
-- (@\"t1\"@, @\"t2\"@, etc.).
addExtension :: DiffField -> Extension -> DiffField
addExtension (DiffField x exts names) ext =
  let name = "t" ++ show (length exts + 1)
  in DiffField x (exts ++ [ext]) (names ++ [name])

-- | Compute the derivation of a field element. Alias for 'deriveExpr'.
derive :: DiffField -> Expr -> Expr
derive field expr = deriveExpr field expr

-- | Compute the derivative of an 'Expr' in a differential field tower.
--
-- Uses the algebraic chain rule, dispatching on each 'Expr'
-- constructor. For extension variables (@Var \"t1\"@, etc.), the
-- derivation is looked up from the tower via 'lookupExtDeriv'.
-- For the base variable, @D(x) = 1@. For constants, @D(c) = 0@.
deriveExpr :: DiffField -> Expr -> Expr
deriveExpr _     (Const _)         = Const 0
deriveExpr _     Pi                = Const 0
deriveExpr _     E                 = Const 0
deriveExpr field (Var x)
  | x == baseVar field             = Const 1
  | otherwise                      = lookupExtDeriv field x
deriveExpr field (Add f g)         =
  Add (deriveExpr field f) (deriveExpr field g)
deriveExpr field (Sub f g)         =
  Sub (deriveExpr field f) (deriveExpr field g)
deriveExpr field (Mul f g)         =
  Add (Mul (deriveExpr field f) g) (Mul f (deriveExpr field g))
deriveExpr field (Div f g)         =
  Div (Sub (Mul (deriveExpr field f) g) (Mul f (deriveExpr field g)))
      (Mul g g)
deriveExpr field (Pow f (Const n)) =
  Mul (Mul (Const n) (Pow f (Const (n-1)))) (deriveExpr field f)
deriveExpr field (Neg f)           =
  Neg (deriveExpr field f)
deriveExpr field (Exp f)           =
  Mul (Exp f) (deriveExpr field f)
deriveExpr field (Log f)           =
  Div (deriveExpr field f) f
deriveExpr field (Sin f)           =
  Mul (Cos f) (deriveExpr field f)
deriveExpr field (Cos f)           =
  Neg (Mul (Sin f) (deriveExpr field f))
deriveExpr field (Arcsin f)        =
  -- D(arcsin(f)) = Df / sqrt(1 - f²)
  Div (deriveExpr field f)
      (Pow (Sub (Const 1) (Mul f f)) (Const 0.5))
deriveExpr field (Arccos f)        =
  -- D(arccos(f)) = -Df / sqrt(1 - f²)
  Neg (Div (deriveExpr field f)
           (Pow (Sub (Const 1) (Mul f f)) (Const 0.5)))
deriveExpr field (Arctan f)        =
  -- D(arctan(f)) = Df / (1 + f²)
  Div (deriveExpr field f)
      (Add (Const 1) (Mul f f))
deriveExpr field (Abs f)           =
  Div (Mul f (deriveExpr field f)) (Abs f)
deriveExpr field (Erf f)           =
  -- D(erf(f)) = (2/√π) · e^(-f²) · Df
  Mul (Mul (Div (Const 2) (Pow Pi (Const 0.5)))
           (Exp (Neg (Mul f f))))
      (deriveExpr field f)
deriveExpr field (Li f)            =
  -- D(li(f)) = Df / log(f)
  Div (deriveExpr field f) (Log f)
deriveExpr field (Si f)            =
  -- D(Si(f)) = sin(f)/f · Df
  Mul (Div (Sin f) f) (deriveExpr field f)
deriveExpr field (Ci f)            =
  -- D(Ci(f)) = cos(f)/f · Df
  Mul (Div (Cos f) f) (deriveExpr field f)
deriveExpr field (Ei f)            =
  -- D(Ei(f)) = e^f/f · Df
  Mul (Div (Exp f) f) (deriveExpr field f)
deriveExpr _     _                 = Const 0

-- | Look up the derivation of an extension variable in the tower.
--
-- Returns @Const 0@ if the variable is not found in the tower
-- (treating it as a constant with respect to the derivation).
lookupExtDeriv :: DiffField -> String -> Expr
lookupExtDeriv field@(DiffField x exts names) name =
  case lookup name (zip names exts) of
    Nothing  -> Const 0
    Just ext -> extensionDeriv x ext name

-- | Compute the derivation of an extension element θ.
--
-- Dispatches on the extension type:
--
-- * 'Primitive': @Dθ = Df\/f@ where @f@ is the argument of @log@.
-- * 'Exponential': @Dθ = θ · Df@ where @f@ is the argument of @exp@.
-- * 'Algebraic': implicit differentiation, @Dθ = −(∂p\/∂x)\/(∂p\/∂θ)@.
-- * 'Special': returns the pre-computed @specialDeriv@.
extensionDeriv :: String -> Extension -> String -> Expr
extensionDeriv _ (Primitive f) _theta =
  Div (deriveBase f) f
extensionDeriv _ (Exponential f) theta =
  Mul (Var theta) (deriveBase f)
extensionDeriv baseX (Algebraic p) theta =
  -- Implicit differentiation: p(θ, x) = 0
  -- ∂p/∂θ · Dθ + ∂p/∂x = 0  =>  Dθ = −(∂p/∂x) / (∂p/∂θ)
  let dpdx     = simplify (deriveExpr (baseField baseX) p)
      dpdtheta = simplify (deriveExpr (baseField theta) p)
  in simplify (Neg (Div dpdx dpdtheta))
extensionDeriv _ (Special _name _arg deriv) _theta =
  deriv

-- | Differentiate an 'Expr' in the base field @ℂ(x)@ only,
-- with no extension variables.
--
-- Equivalent to @'deriveExpr' ('baseField' \"x\")@. Used by
-- 'LimCalc.Risch' and 'LimCalc.SymExpand' to compute derivatives
-- of Taylor series coefficients.
deriveBase :: Expr -> Expr
deriveBase expr = deriveExpr (baseField "x") expr

-- | Build a differential field tower for an expression by collecting
-- all extensions it requires.
buildTower :: String -> Expr -> DiffField
buildTower x expr = foldr addExt (baseField x) (collectExts expr)
  where
    addExt ext field = addExtension field ext

-- | Collect all extensions required to express an 'Expr' in the
-- differential field tower.
--
-- Returns one 'Extension' per transcendental or special-function
-- subexpression. Duplicates may appear and are not deduplicated here.
collectExts :: Expr -> [Extension]
collectExts (Const _)  = []
collectExts (Var _)    = []
collectExts Pi         = []
collectExts E          = []
collectExts I          = []
collectExts (Add f g)  = collectExts f ++ collectExts g
collectExts (Sub f g)  = collectExts f ++ collectExts g
collectExts (Mul f g)  = collectExts f ++ collectExts g
collectExts (Div f g)  = collectExts f ++ collectExts g
collectExts (Pow f g)  = collectExts f ++ collectExts g
collectExts (Neg f)    = collectExts f
collectExts (Abs f)    = collectExts f
collectExts (Exp f)    = Exponential f : collectExts f
collectExts (Log f)    = Primitive f   : collectExts f
collectExts (Sin f)    =
  -- sin(f) = Im(exp(i·f)), so add exp(i·f) as extension
  Exponential (Mul I f) : collectExts f
collectExts (Cos f)    =
  Exponential (Mul I f) : collectExts f
collectExts (Erf f)    = specialExt "erf" f : collectExts f
collectExts (Li f)     = specialExt "li"  f : collectExts f
collectExts (Si f)     = specialExt "Si"  f : collectExts f
collectExts (Ci f)     = specialExt "Ci"  f : collectExts f
collectExts (Ei f)     = specialExt "Ei"  f : collectExts f
collectExts (Arcsin f) = collectExts f
collectExts (Arccos f) = collectExts f
collectExts (Arctan f) = collectExts f

-- | Construct a 'Special' extension for a named special function.
--
-- Looks up the function in 'knownSpecial'. Errors if the name is not
-- registered — this should not happen for the five names used in
-- 'collectExts'.
specialExt :: String -> Expr -> Extension
specialExt name arg = case lookupSpecial name arg of
  Just ext -> ext
  Nothing  -> error ("specialExt: unknown special function " ++ name)

-- | Test whether an 'Expr' belongs to the base field @ℂ(x)@
-- (no transcendental or special-function subexpressions).
isBaseField :: String -> Expr -> Bool
isBaseField x (Const _) = True
isBaseField x (Var v)   = v == x
isBaseField _ Pi        = True
isBaseField _ E         = True
isBaseField x (Add f g) = isBaseField x f && isBaseField x g
isBaseField x (Sub f g) = isBaseField x f && isBaseField x g
isBaseField x (Mul f g) = isBaseField x f && isBaseField x g
isBaseField x (Div f g) = isBaseField x f && isBaseField x g
isBaseField x (Pow f g) = isBaseField x f && isBaseField x g
isBaseField x (Neg f)   = isBaseField x f
isBaseField _ _         = False

-- | Registry of known special functions and their derivative rules.
--
-- Each entry is a pair @(name, constructor)@ where @constructor arg@
-- builds a 'Special' extension for the given argument expression.
knownSpecial :: [(String, Expr -> Extension)]
knownSpecial =
  [ ("erf", \arg -> Special "erf" arg
      (Mul (Div (Const 2) (Pow Pi (Const 0.5)))
           (Exp (Neg (Mul arg arg)))))
  , ("li",  \arg -> Special "li" arg
      (Div (Const 1) (Log arg)))
  , ("Si",  \arg -> Special "Si" arg
      (Div (Sin arg) arg))
  , ("Ci",  \arg -> Special "Ci" arg
      (Div (Cos arg) arg))
  , ("Ei",  \arg -> Special "Ei" arg
      (Div (Exp arg) arg))
  ]

-- | Look up a special function by name in 'knownSpecial'.
lookupSpecial :: String -> Expr -> Maybe Extension
lookupSpecial name arg =
  fmap (\f -> f arg) (lookup name knownSpecial)