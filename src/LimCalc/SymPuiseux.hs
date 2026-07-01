-- | Symbolic Puiseux series with 'Expr' coefficients.
--
-- 'SymPuiseuxSeries' is a Puiseux series whose coefficients are
-- symbolic 'Expr' values rather than numeric 'AlgNum' values. This
-- is the representation used by the symbolic differentiation path
-- in 'LimCalc.SymExpand': expanding @f(x + h)@ symbolically yields
-- a series whose coefficients are expressions in @x@, from which the
-- derivative (the @h^1@ coefficient) can be read off as an 'Expr'.
--
-- Unlike 'LimCalc.Puiseux.LogPuiseuxSeries', this type supports only
-- pure power terms (no @log(h)@ terms), since the symbolic expansion
-- path is used for differentiation only and the functions it handles
-- are all analytic at generic points.
--
-- Coefficients are /not/ simplified automatically; callers apply
-- 'LimCalc.Simplify.simplify' to the extracted coefficient as needed.
module LimCalc.SymPuiseux
  ( -- * Types
    SymPuiseuxTerm (..)
  , SymPuiseuxSeries (..)
    -- * Constructors
  , zeroSym
    -- * Normalisation
  , normalizeSym
    -- * Arithmetic
  , addSymSeries
  , scaleSymSeries
  , mulSymSeries
  , combineLikeSym
    -- * Accessors
  , symCoeffAt
  , symLeadingTerm
    -- * Structural operations
  , truncateSymSeries
  , shiftSymExponents
  , removeSymTerm
  ) where

import Data.List (sortBy)
import Data.Ord (comparing)
import LimCalc.Expr

-- | A single term in a symbolic Puiseux series: @c(x) · h^p@.
--
-- The coefficient 'symCoeff' is a symbolic 'Expr' (typically a
-- function of the base variable @x@). The exponent 'symExp' is
-- rational, allowing fractional powers.
data SymPuiseuxTerm = SymPuiseuxTerm
  { symExp   :: Rational
    -- ^ Exponent of @h@.
  , symCoeff :: Expr
    -- ^ Coefficient as a symbolic expression.
  } deriving (Show, Eq)

-- | A symbolic Puiseux series: a finite list of 'SymPuiseuxTerm'
-- values in ascending exponent order.
--
-- Canonical form: terms are sorted by 'symExp' in ascending order.
-- Terms with equal exponents are combined by 'combineLikeSym'.
-- Unlike 'LimCalc.Puiseux.LogPuiseuxSeries', zero coefficients are
-- not stripped automatically (since symbolic zero-detection would
-- require a full simplifier call).
newtype SymPuiseuxSeries = SymPuiseuxSeries
  { symTerms :: [SymPuiseuxTerm]
    -- ^ Terms in ascending exponent order.
  } deriving (Show, Eq)

-- | The zero symbolic series (empty term list).
zeroSym :: SymPuiseuxSeries
zeroSym = SymPuiseuxSeries []

-- | Sort terms into ascending exponent order.
normalizeSym :: SymPuiseuxSeries -> SymPuiseuxSeries
normalizeSym (SymPuiseuxSeries ts) =
  SymPuiseuxSeries $ sortBy (comparing symExp) ts

-- | Merge two symbolic series, combining coefficients of terms with
-- equal exponents by addition.
--
-- Both inputs must already be in canonical order. The output is
-- normalised but coefficients are not simplified.
combineLikeSym :: SymPuiseuxSeries -> SymPuiseuxSeries -> SymPuiseuxSeries
combineLikeSym (SymPuiseuxSeries s1) (SymPuiseuxSeries s2) =
  normalizeSym $ SymPuiseuxSeries $ mergePlus s1 s2
  where
    mergePlus [] ys = ys
    mergePlus xs [] = xs
    mergePlus (x:xs) (y:ys)
      | symExp x == symExp y =
          SymPuiseuxTerm (symExp x) (Add (symCoeff x) (symCoeff y))
            : mergePlus xs ys
      | symExp x < symExp y  = x : mergePlus xs (y:ys)
      | otherwise            = y : mergePlus (x:xs) ys

-- | Add two symbolic series.
addSymSeries :: SymPuiseuxSeries -> SymPuiseuxSeries -> SymPuiseuxSeries
addSymSeries = combineLikeSym

-- | Scale every term's coefficient by a symbolic expression.
scaleSymSeries :: Expr -> SymPuiseuxSeries -> SymPuiseuxSeries
scaleSymSeries c (SymPuiseuxSeries ts) =
  SymPuiseuxSeries $ map (\t -> t { symCoeff = Mul c (symCoeff t) }) ts

-- | Multiply two symbolic series (Cauchy product).
--
-- Computes the cross-product of all term pairs, sorts by exponent,
-- then combines like terms by summing their coefficients.
mulSymSeries :: SymPuiseuxSeries -> SymPuiseuxSeries -> SymPuiseuxSeries
mulSymSeries (SymPuiseuxSeries s1) (SymPuiseuxSeries s2) =
  let raw = [ SymPuiseuxTerm (symExp t1 + symExp t2)
                              (Mul (symCoeff t1) (symCoeff t2))
            | t1 <- s1, t2 <- s2 ]
      sorted = normalizeSym (SymPuiseuxSeries raw)
  in foldTerms sorted
  where
    foldTerms (SymPuiseuxSeries []) = SymPuiseuxSeries []
    foldTerms (SymPuiseuxSeries (t:ts)) =
      let (same, rest) = span (\t' -> symExp t' == symExp t) ts
          combined = foldl (\acc t' -> Add acc (symCoeff t')) (symCoeff t) same
          SymPuiseuxSeries tl = foldTerms (SymPuiseuxSeries rest)
      in SymPuiseuxSeries (SymPuiseuxTerm (symExp t) combined : tl)

-- | Extract the coefficient of @h^n@.
--
-- If multiple terms share the exponent @n@ (which can happen before
-- 'combineLikeSym' is applied), their coefficients are summed. Returns
-- @Const 0@ if no term has exponent @n@.
symCoeffAt :: Rational -> SymPuiseuxSeries -> Expr
symCoeffAt n (SymPuiseuxSeries ts) =
  case filter (\t -> symExp t == n) ts of
    []     -> Const 0
    [t]    -> symCoeff t
    (t:ts) -> foldl (\acc t' -> Add acc (symCoeff t')) (symCoeff t) ts

-- | Return the leading term (smallest exponent), or 'Nothing' if the
-- series is empty.
symLeadingTerm :: SymPuiseuxSeries -> Maybe SymPuiseuxTerm
symLeadingTerm (SymPuiseuxSeries [])    = Nothing
symLeadingTerm (SymPuiseuxSeries (t:_)) = Just t

-- | Retain only the first @n@ terms (by position in the sorted list).
truncateSymSeries :: Int -> SymPuiseuxSeries -> SymPuiseuxSeries
truncateSymSeries n (SymPuiseuxSeries ts) = SymPuiseuxSeries (take n ts)

-- | Add @delta@ to every term's exponent.
shiftSymExponents :: Rational -> SymPuiseuxSeries -> SymPuiseuxSeries
shiftSymExponents delta (SymPuiseuxSeries ts) =
  SymPuiseuxSeries
    [ SymPuiseuxTerm (symExp t + delta) (symCoeff t) | t <- ts ]

-- | Remove all terms with a given exponent.
removeSymTerm :: Rational -> SymPuiseuxSeries -> SymPuiseuxSeries
removeSymTerm e (SymPuiseuxSeries ts) =
  SymPuiseuxSeries $ filter (\t -> symExp t /= e) ts