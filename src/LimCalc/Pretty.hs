module LimCalc.Pretty where

import LimCalc.Expr
import Data.Ratio (Ratio, numerator, denominator, (%))

-- | Pretty-print an Expr as a human-readable string with minimal parentheses.
prettyExpr :: Expr -> String
prettyExpr expr = prettyPrec 0 expr

-- Precedence levels:
--   0 = top level / Add / Sub
--   1 = Mul / Div
--   2 = Neg (unary minus)
--   3 = Pow
--   4 = atoms (Const, Var, function application)

prettyPrec :: Int -> Expr -> String
prettyPrec _ (Const d)   = prettyDouble d
prettyPrec _ Pi          = "π"
prettyPrec _ E           = "e"
prettyPrec _ I           = "i"
prettyPrec _ (Var x)     = x

prettyPrec p (Add a b)   = parensIf (p > 0) $
  prettyPrec 0 a ++ " + " ++ prettyPrec 0 b

prettyPrec p (Sub a b)   = parensIf (p > 0) $
  prettyPrec 0 a ++ " - " ++ prettyPrec 1 b

prettyPrec p (Mul a b)   = parensIf (p > 1) $
  prettyPrec 1 a ++ " * " ++ prettyPrec 1 b

prettyPrec p (Div a b)   = parensIf (p > 1) $
  prettyPrec 1 a ++ " / " ++ prettyPrec 2 b

prettyPrec _ (Pow b (Const 0.5)) = "sqrt(" ++ prettyPrec 0 b ++ ")"
prettyPrec p (Pow b e)   = parensIf (p > 3) $
  prettyPrec 4 b ++ "^" ++ prettyPrec 4 e

prettyPrec p (Neg (Neg x)) = prettyPrec p x
prettyPrec _ (Neg x)     = "-" ++ prettyPrec 4 x

prettyPrec _ (Abs x)     = "|" ++ prettyPrec 0 x ++ "|"
prettyPrec _ (Exp x)     = "exp(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Log x)     = "log(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Sin x)     = "sin(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Cos x)     = "cos(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Arcsin x)  = "arcsin(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Arccos x)  = "arccos(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Arctan x)  = "arctan(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Erf x)     = "erf(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Li x)      = "li(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Si x)      = "Si(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Ci x)      = "Ci(" ++ prettyPrec 0 x ++ ")"
prettyPrec _ (Ei x)      = "Ei(" ++ prettyPrec 0 x ++ ")"

parensIf :: Bool -> String -> String
parensIf True  s = "(" ++ s ++ ")"
parensIf False s = s

-- | Render a Double as p/q if it snaps cleanly, else as decimal.
prettyDouble :: Double -> String
prettyDouble d
  | d == fromIntegral n = show n
  | otherwise = case snapToRational d of
      Just r  -> show (numerator r) ++ "/" ++ show (denominator r)
      Nothing -> show d
  where n = round d :: Int

-- | Try to express a Double as p/q with |q| <= 100.
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