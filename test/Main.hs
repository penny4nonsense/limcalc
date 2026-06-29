module Main where

import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Map.Strict as Map
import LimCalc.Expr
import LimCalc.Puiseux
import LimCalc.Expand
import LimCalc.Calculus
import LimCalc.Simplify
import LimCalc.Limit
import LimCalc.Types

main :: IO ()
main = defaultMain tests

tests :: TestTree
tests = testGroup "limcalc"
  [ seriesTests
  , expandTests
  , derivativeTests
  , simplifyTests
  , limitTests
  ]

-- | Helper: expand at a point
pt :: Double -> Map.Map String Double
pt x = Map.fromList [("x", x)]

-- | Helper: get Right value
unRight :: Either a b -> b
unRight (Right x) = x
unRight (Left _)  = error "Expected Right"

------------------------------------------------------------------------
-- Series arithmetic
------------------------------------------------------------------------

seriesTests :: TestTree
seriesTests = testGroup "PuiseuxSeries arithmetic"
  [ testCase "addSeries combines like terms" $
      let s1 = PuiseuxSeries [PuiseuxTerm 1 2.0, PuiseuxTerm 2 3.0]
          s2 = PuiseuxSeries [PuiseuxTerm 1 4.0, PuiseuxTerm 3 1.0]
          r  = addSeries s1 s2
      in terms r @?= [ PuiseuxTerm 1 6.0
                     , PuiseuxTerm 2 3.0
                     , PuiseuxTerm 3 1.0 ]

  , testCase "addSeries cancels zero terms" $
      let s1 = PuiseuxSeries [PuiseuxTerm 1 1.0]
          s2 = PuiseuxSeries [PuiseuxTerm 1 (-1.0)]
      in terms (addSeries s1 s2) @?= []

  , testCase "mulSeries Cauchy product" $
      let s1 = PuiseuxSeries [PuiseuxTerm 0 1.0, PuiseuxTerm 1 1.0]
          s2 = PuiseuxSeries [PuiseuxTerm 0 1.0, PuiseuxTerm 1 1.0]
          r  = mulSeries s1 s2
      in terms r @?= [ PuiseuxTerm 0 1.0
                     , PuiseuxTerm 1 2.0
                     , PuiseuxTerm 2 1.0 ]

  , testCase "mulSeries with fractional exponents" $
      let s1 = PuiseuxSeries [PuiseuxTerm (1/2) 1.0]
          s2 = PuiseuxSeries [PuiseuxTerm (1/2) 1.0]
          r  = mulSeries s1 s2
      in terms r @?= [PuiseuxTerm 1 1.0]
  ]

------------------------------------------------------------------------
-- Expansion engine
------------------------------------------------------------------------

expandTests :: TestTree
expandTests = testGroup "Expansion engine"
  [ testCase "expand Const" $
      expand (Const 3.0) (pt 0) "x"
        @?= Right (PuiseuxSeries [PuiseuxTerm 0 3.0])

  , testCase "expand Var at x0=0" $
      expand (Var "x") (pt 0) "x"
        @?= Right (PuiseuxSeries [PuiseuxTerm 1 1.0])

  , testCase "expand Var at x0=2" $
      expand (Var "x") (pt 2) "x"
        @?= Right (PuiseuxSeries [PuiseuxTerm 0 2.0, PuiseuxTerm 1 1.0])

  , testCase "expand sin(x) at x0=0 leading term" $
      let Right s = expand (Sin (Var "x")) (pt 0) "x"
          Just lt = leadingTermNZ s
      in (pExp lt, coeff lt) @?= (1, 1.0)

  , testCase "expand sin(x) at x0=pi/2 constant term" $
      let Right s = expand (Sin (Var "x")) (pt (pi/2)) "x"
      in abs (constantTerm s - 1.0) < 1e-10 @?= True

  , testCase "expand exp(x) at x0=0 constant term" $
      let Right s = expand (Exp (Var "x")) (pt 0) "x"
      in constantTerm s @?= 1.0

  , testCase "expand 1/x at x0=0 is pole" $
      let Right s = expand (Div (Const 1) (Var "x")) (pt 0) "x"
          Just lt = leadingTermNZ s
      in pExp lt @?= (-1)

  , testCase "expand sin^(1/2)(x) at x0=0 leading exponent" $
      let Right s = expand (Pow (Sin (Var "x")) (Const 0.5)) (pt 0) "x"
          Just lt = leadingTermNZ s
      in pExp lt @?= (1/2)

  , testCase "expand sin^(1/2)(x) at x0=0 leading coeff" $
      let Right s = expand (Pow (Sin (Var "x")) (Const 0.5)) (pt 0) "x"
          Just lt = leadingTermNZ s
      in coeff lt @?= 1.0

  , testCase "expand log(x) at x0<=0 returns error" $
      case expand (Log (Var "x")) (pt 0) "x" of
        Left (Undefined _) -> return ()
        _                  -> assertFailure "Expected Undefined error"

  , testCase "expand unknown variable returns error" $
      case expand (Var "y") (pt 0) "x" of
        Left (Unknown _) -> return ()
        _                -> assertFailure "Expected Unknown error"
  ]

------------------------------------------------------------------------
-- Derivative tests
------------------------------------------------------------------------

derivativeTests :: TestTree
derivativeTests = testGroup "Derivatives"
  [ testCase "d/dx sin(x) at x=0 is 1" $
      unRight (derivative (Sin (Var "x")) (pt 0) "x") @?= 1.0

  , testCase "d/dx sin(x) at x=pi/2 is ~0" $
      abs (unRight (derivative (Sin (Var "x")) (pt (pi/2)) "x")) < 1e-10
        @?= True

  , testCase "d/dx x^2 at x=3 is 6" $
      unRight (derivative (Pow (Var "x") (Const 2)) (pt 3) "x") @?= 6.0

  , testCase "d/dx exp(x) at x=0 is 1" $
      unRight (derivative (Exp (Var "x")) (pt 0) "x") @?= 1.0

  , testCase "d/dx log(x) at x=1 is 1" $
      unRight (derivative (Log (Var "x")) (pt 1) "x") @?= 1.0

  , testCase "d/dx cos(x) at x=0 is 0" $
      abs (unRight (derivative (Cos (Var "x")) (pt 0) "x")) < 1e-10
        @?= True

  , testCase "d/dx sin^(1/2)(x) at x=1 matches analytic" $
      let result   = unRight (derivative (Pow (Sin (Var "x")) (Const 0.5)) (pt 1) "x")
          expected = cos 1 / (2 * sqrt (sin 1))
      in abs (result - expected) < 1e-10 @?= True

  , testCase "symbolic d/dx sin(x) is cos(x)" $
      fmap simplify (diff (Sin (Var "x")) "x")
        @?= Right (Cos (Var "x"))

  , testCase "symbolic d/dx cos(x) is -sin(x)" $
      fmap simplify (diff (Cos (Var "x")) "x")
        @?= Right (Neg (Sin (Var "x")))

  , testCase "symbolic d/dx exp(x) is exp(x)" $
      fmap simplify (diff (Exp (Var "x")) "x")
        @?= Right (Exp (Var "x"))

  , testCase "symbolic d/dx log(x) is 1/x" $
      fmap simplify (diff (Log (Var "x")) "x")
        @?= Right (Div (Const 1.0) (Var "x"))

  , testCase "symbolic d/dx x^2 is 2x" $
      fmap simplify (diff (Pow (Var "x") (Const 2)) "x")
        @?= Right (Mul (Const 2.0) (Var "x"))

  , testCase "symbolic d/dx x^3 is 3x^2" $
      fmap simplify (diff (Pow (Var "x") (Const 3)) "x")
        @?= Right (Mul (Const 3.0) (Pow (Var "x") (Const 2.0)))

  , testCase "product rule: d/dx (x*sin(x))" $
      fmap simplify (diff (Mul (Var "x") (Sin (Var "x"))) "x")
        @?= Right (Add (Mul (Var "x") (Cos (Var "x"))) (Sin (Var "x")))

  , testCase "chain rule: d/dx sin(x^2)" $
      fmap simplify (diff (Sin (Pow (Var "x") (Const 2))) "x")
        @?= Right (Mul (Const 2.0) (Mul (Cos (Pow (Var "x") (Const 2.0))) (Var "x")))
  ]

------------------------------------------------------------------------
-- Simplifier tests
------------------------------------------------------------------------

simplifyTests :: TestTree
simplifyTests = testGroup "Simplifier"
  [ testCase "simplify 0 + x = x" $
      simplify (Add (Const 0) (Var "x")) @?= Var "x"

  , testCase "simplify x + 0 = x" $
      simplify (Add (Var "x") (Const 0)) @?= Var "x"

  , testCase "simplify 1 * x = x" $
      simplify (Mul (Const 1) (Var "x")) @?= Var "x"

  , testCase "simplify x * 1 = x" $
      simplify (Mul (Var "x") (Const 1)) @?= Var "x"

  , testCase "simplify 0 * x = 0" $
      simplify (Mul (Const 0) (Var "x")) @?= Const 0

  , testCase "simplify x - x = 0" $
      simplify (Sub (Var "x") (Var "x")) @?= Const 0

  , testCase "simplify x / x = 1" $
      simplify (Div (Var "x") (Var "x")) @?= Const 1

  , testCase "simplify x^0 = 1" $
      simplify (Pow (Var "x") (Const 0)) @?= Const 1

  , testCase "simplify x^1 = x" $
      simplify (Pow (Var "x") (Const 1)) @?= Var "x"

  , testCase "simplify --x = x" $
      simplify (Neg (Neg (Var "x"))) @?= Var "x"

  , testCase "simplify 2 * 3 = 6" $
      simplify (Mul (Const 2) (Const 3)) @?= Const 6

  , testCase "simplify 2 + 3 = 5" $
      simplify (Add (Const 2) (Const 3)) @?= Const 5
  ]

------------------------------------------------------------------------
-- Limit tests
------------------------------------------------------------------------

limitTests :: TestTree
limitTests = testGroup "Limits"
  [ testCase "lim_{x->0} sin(x)/x = 1" $
      limit (Div (Sin (Var "x")) (Var "x")) "x" 0
        @?= Exists 1.0

  , testCase "lim_{x->0} x^2 = 0" $
      limit (Pow (Var "x") (Const 2)) "x" 0
        @?= Exists 0.0

  , testCase "lim_{x->2} x^2 = 4" $
      limit (Pow (Var "x") (Const 2)) "x" 2
        @?= Exists 4.0

  , testCase "lim_{x->0} 1/x is pole" $
      case limit (Div (Const 1) (Var "x")) "x" 0 of
        Pole _ -> return ()
        _      -> assertFailure "Expected Pole"

  , testCase "lim_{x->0} exp(x) = 1" $
      limit (Exp (Var "x")) "x" 0
        @?= Exists 1.0

  , testCase "lim_{x->0} (1-cos(x))/x^2 = 1/2" $
      let f = Div (Sub (Const 1) (Cos (Var "x")))
                  (Pow (Var "x") (Const 2))
          Exists v = limit f "x" 0
      in abs (v - 0.5) < 1e-10 @?= True
  ]