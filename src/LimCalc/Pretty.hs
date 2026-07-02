-- | Pretty-printing of 'Expr' trees.
--
-- Produces human-readable mathematical notation with minimal
-- parentheses, using standard precedence rules. The output is
-- intended for display in REPLs, test output, and paper examples —
-- not for round-tripping back to 'Expr' (no parser exists yet).
--
-- = Precedence levels
--
-- @
-- 0  —  top level, 'Add', 'Sub'
-- 1  —  'Mul', 'Div'
-- 2  —  unary 'Neg'
-- 3  —  'Pow'
-- 4  —  atoms: 'Const', 'Var', function application
-- @
--
-- Associativity is handled by passing different precedences to the
-- left and right sub-expressions:
--
-- * @Sub@: right operand gets precedence 1 (not 0), so
--   @a − (b − c)@ parenthesises the right operand correctly.
-- * @Div@: right operand gets precedence 2, so @a \/ (b \/ c)@
--   parenthesises correctly.
--
-- = Numeric rendering
--
-- 'Const' values are rendered via 'prettyDouble':
--
-- * Integer-valued doubles print without a decimal point (@2.0@ → @\"2\"@).
-- * Doubles that snap to a small rational (@|q| ≤ 100@) print as
--   @p\/q@ (@0.333…@ → @\"1\/3\"@).
-- * Other doubles print via 'show'.
module LimCalc.Pretty
  ( -- * Pretty-printing
    prettyExpr
  , prettyPrec
    -- * Numeric rendering
  , prettyDouble
  , snapToRational
    -- * Utilities
  , parensIf
  ) where

import LimCalc.Core.Expr
import Data.Ratio (Ratio, numerator, denominator, (%))

-- | Pretty-print an 'Expr' with minimal parentheses.
prettyExpr :: Expr -> String
prettyExpr = prettyPrec 0

-- | Pretty-print with an ambient precedence level.
-- Wraps the result in parentheses when the expression's own
-- precedence is lower than the ambient level.
prettyPrec :: Int -> Expr -> String
prettyPrec _ (Const d)          = prettyDouble d
prettyPrec _ Pi                 = "π"
prettyPrec _ E                  = "e"
prettyPrec _ I                  = "i"
prettyPrec _ (Var x)            = x

prettyPrec p (Add a b)          = parensIf (p > 0) $
  prettyPrec 0 a ++ " + " ++ prettyPrec 0 b

prettyPrec p (Sub a b)          = parensIf (p > 0) $
  prettyPrec 0 a ++ " - " ++ prettyPrec 1 b

prettyPrec p (Mul a b)          = parensIf (p > 1) $
  prettyPrec 1 a ++ " * " ++ prettyPrec 1 b

prettyPrec p (Div a b)          = parensIf (p > 1) $
  prettyPrec 1 a ++ " / " ++ prettyPrec 2 b

prettyPrec _ (Pow b (Const 0.5)) = "sqrt(" ++ prettyPrec 0 b ++ ")"
prettyPrec p (Pow b e)           = parensIf (p > 3) $
  prettyPrec 4 b ++ "^" ++ prettyPrec 4 e

prettyPrec p (Neg (Neg x))      = prettyPrec p x
prettyPrec _ (Neg x)            = "-" ++ prettyPrec 4 x

prettyPrec _ (Abs x)            = "|" ++ prettyPrec 0 x ++ "|"
prettyPrec _ (Exp x)            = "exp(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Log x)            = "log(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Sin x)            = "sin(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Cos x)            = "cos(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Arcsin x)         = "arcsin(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Arccos x)         = "arccos(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Arctan x)         = "arctan(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Erf x)            = "erf(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Li x)             = "li(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Si x)             = "Si(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Ci x)             = "Ci(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Ei x)             = "Ei(" ++ prettyPrec 0 x ++ ")"

-- | Wrap a string in parentheses if the condition holds.
parensIf :: Bool -> String -> String
parensIf True  s = "(" ++ s ++ ")"
parensIf False s = s

-- | Render a 'Double' as a human-readable string.
--
-- * Integer-valued: @show n@ (no decimal point).
-- * Snaps to @p\/q@ with @|q| ≤ 100@: @\"p\/q\"@.
-- * Otherwise: @show d@.
prettyDouble :: Double -> String
prettyDouble d
  | d == fromIntegral n = show n
  | otherwise           = case snapToRational d of
      Just r  -> show (numerator r) ++ "/" ++ show (denominator r)
      Nothing -> show d
  where
    n = round d :: Int

-- | Try to express a 'Double' as @p\/q@ with @|q| ≤ 100@,
-- returning 'Nothing' if no such representation exists within
-- @1e-10@ absolute tolerance.
snapToRational :: Double -> Maybe (Ratio Integer)
snapToRational d = go 1
  where
    go q
      | q > 100   = Nothing
      | otherwise =
          let p = round (d * fromIntegral q) :: Integer
              r = p % q
          in if abs (fromRational r - d) < 1e-10
               then Just r
               else go (q + 1)