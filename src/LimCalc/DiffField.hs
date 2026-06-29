module LimCalc.DiffField where

import LimCalc.Expr
import LimCalc.Poly
import LimCalc.RationalFunction

-- | A single extension of the differential field
-- Each extension adds a new element θ to the field with a known derivation
data Extension
  = Primitive Expr        -- ^ θ = log(f), Dθ = Df/f
  | Exponential Expr      -- ^ θ = exp(f), Dθ = θ·Df
  | Algebraic Expr        -- ^ minimal polynomial p(θ) = 0
  | Special               -- ^ special function with known derivative
      { specialName  :: String  -- ^ "erf", "li", "Si", etc.
      , specialArg   :: Expr    -- ^ the argument
      , specialDeriv :: Expr    -- ^ D(f(x)) as an Expr
      }
  deriving (Show, Eq)

-- | The differential field tower
-- Base field: ℂ(x) — rational functions in the base variable
-- Each extension adds a new transcendental or algebraic element
data DiffField = DiffField
  { baseVar    :: String        -- ^ base variable (usually "x")
  , extensions :: [Extension]   -- ^ tower of extensions, innermost first
  , extNames   :: [String]      -- ^ variable names for extensions ("t1", "t2", ...)
  } deriving (Show, Eq)

-- | Empty field — just ℂ(x)
baseField :: String -> DiffField
baseField x = DiffField x [] []

-- | Add an extension to the tower
addExtension :: DiffField -> Extension -> DiffField
addExtension (DiffField x exts names) ext =
  let name = "t" ++ show (length exts + 1)
  in DiffField x (exts ++ [ext]) (names ++ [name])

-- | The derivation of a field element
-- Given the tower, compute D(expr) using the chain rule
derive :: DiffField -> Expr -> Expr
derive field expr = deriveExpr field expr

-- | Compute the derivative of an expression in the differential field
deriveExpr :: DiffField -> Expr -> Expr
deriveExpr field (Const _)   = Const 0
deriveExpr field Pi          = Const 0
deriveExpr field E           = Const 0
deriveExpr field (Var x)
  | x == baseVar field = Const 1
  | otherwise          = lookupExtDeriv field x
deriveExpr field (Add f g)   =
  Add (deriveExpr field f) (deriveExpr field g)
deriveExpr field (Sub f g)   =
  Sub (deriveExpr field f) (deriveExpr field g)
deriveExpr field (Mul f g)   =
  Add (Mul (deriveExpr field f) g) (Mul f (deriveExpr field g))
deriveExpr field (Div f g)   =
  Div (Sub (Mul (deriveExpr field f) g) (Mul f (deriveExpr field g)))
      (Mul g g)
deriveExpr field (Pow f (Const n)) =
  Mul (Mul (Const n) (Pow f (Const (n-1)))) (deriveExpr field f)
deriveExpr field (Neg f)     =
  Neg (deriveExpr field f)
deriveExpr field (Exp f)     =
  Mul (Exp f) (deriveExpr field f)
deriveExpr field (Log f)     =
  Div (deriveExpr field f) f
deriveExpr field (Sin f)     =
  Mul (Cos f) (deriveExpr field f)
deriveExpr field (Cos f)     =
  Neg (Mul (Sin f) (deriveExpr field f))
deriveExpr field (Abs f)     =
  Div (Mul f (deriveExpr field f)) (Abs f)
deriveExpr field _           = Const 0

-- | Look up the derivation of an extension variable
lookupExtDeriv :: DiffField -> String -> Expr
lookupExtDeriv (DiffField x exts names) name =
  case lookup name (zip names exts) of
    Nothing  -> Const 0  -- not in tower, treat as constant
    Just ext -> extensionDeriv ext name

-- | The derivation of an extension element θ
extensionDeriv :: Extension -> String -> Expr
extensionDeriv (Primitive f) theta =
  -- Dθ = Df/f
  Div (deriveBase f) f
extensionDeriv (Exponential f) theta =
  -- Dθ = θ · Df
  Mul (Var theta) (deriveBase f)
extensionDeriv (Algebraic p) theta =
  -- Dθ via implicit differentiation: Dp(θ) = 0
  -- p'(θ)·Dθ + ∂p/∂x = 0 → Dθ = -(∂p/∂x)/p'(θ)
  Const 0  -- TODO: implement properly
extensionDeriv (Special name arg deriv) theta =
  -- D(f(arg)) = deriv · D(arg)
  deriv

-- | Derive in the base field ℂ(x) only
deriveBase :: Expr -> Expr
deriveBase expr = deriveExpr (baseField "x") expr

-- | Build a differential field tower from an Expr
-- Analyzes the expression and builds the appropriate tower
buildTower :: String -> Expr -> DiffField
buildTower x expr = foldr addExt (baseField x) (collectExts expr)
  where
    addExt ext field = addExtension field ext

-- | Collect extensions needed for an expression
collectExts :: Expr -> [Extension]
collectExts (Const _)   = []
collectExts (Var _)     = []
collectExts Pi          = []
collectExts E           = []
collectExts I           = []
collectExts (Add f g)   = collectExts f ++ collectExts g
collectExts (Sub f g)   = collectExts f ++ collectExts g
collectExts (Mul f g)   = collectExts f ++ collectExts g
collectExts (Div f g)   = collectExts f ++ collectExts g
collectExts (Pow f g)   = collectExts f ++ collectExts g
collectExts (Neg f)     = collectExts f
collectExts (Abs f)     = collectExts f
collectExts (Exp f)     = Exponential f : collectExts f
collectExts (Log f)     = Primitive f   : collectExts f
collectExts (Sin f)     =
  -- sin(f) = Im(exp(i·f)), so add exp(i·f) as extension
  Exponential (Mul I f) : collectExts f
collectExts (Cos f)     =
  -- cos(f) = Re(exp(i·f))
  Exponential (Mul I f) : collectExts f

-- | Check if an expression is in the base field ℂ(x)
isBaseField :: String -> Expr -> Bool
isBaseField x (Const _)     = True
isBaseField x (Var v)       = v == x
isBaseField x Pi            = True
isBaseField x E             = True
isBaseField x (Add f g)     = isBaseField x f && isBaseField x g
isBaseField x (Sub f g)     = isBaseField x f && isBaseField x g
isBaseField x (Mul f g)     = isBaseField x f && isBaseField x g
isBaseField x (Div f g)     = isBaseField x f && isBaseField x g
isBaseField x (Pow f g)     = isBaseField x f && isBaseField x g
isBaseField x (Neg f)       = isBaseField x f
isBaseField x _             = False

-- | Known special functions and their derivatives
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

-- | Look up a special function by name
lookupSpecial :: String -> Expr -> Maybe Extension
lookupSpecial name arg =
  fmap (\f -> f arg) (lookup name knownSpecial)