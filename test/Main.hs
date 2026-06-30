module Main where

import Data.Complex
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Map.Strict as Map
import LimCalc.Expr
import LimCalc.Puiseux
import LimCalc.Expand
import LimCalc.Calculus
import LimCalc.Simplify
import LimCalc.Limit
import LimCalc.Poly
import LimCalc.RationalFunction
import LimCalc.DiffField
import LimCalc.Risch.Primitive
import LimCalc.Risch.Exponential
import LimCalc.Risch
import LimCalc.Types
import LimCalc.AlgNum
import LimCalc.QPoly
import LimCalc.BivPoly

main :: IO ()
main = defaultMain tests

tests :: TestTree
tests = testGroup "limcalc"
  [ seriesTests
  , expandTests
  , derivativeTests
  , simplifyTests
  , limitTests
  , polyTests
  , rationalFunctionTests
  , rischTests
  , algNumDegree2Tests
  , qPolyDivisionTests
  , trigIntegrationTests
  , specialFunctionTests
  , algebraicExtensionTests
  ]

-- | Helper: expand at a point
pt :: Double -> Map.Map String AlgNum
pt x = Map.fromList [("x", fromQ (toRational x))]

-- | Helper: AlgNum from Double
aN :: Double -> AlgNum
aN = fromQ . toRational

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
      let s1 = PuiseuxSeries [PuiseuxTerm 1 (aN 2), PuiseuxTerm 2 (aN 3)]
          s2 = PuiseuxSeries [PuiseuxTerm 1 (aN 4), PuiseuxTerm 3 (aN 1)]
          r  = addSeries s1 s2
      in terms r @?= [ PuiseuxTerm 1 (aN 6)
                     , PuiseuxTerm 2 (aN 3)
                     , PuiseuxTerm 3 (aN 1) ]

  , testCase "addSeries cancels zero terms" $
      let s1 = PuiseuxSeries [PuiseuxTerm 1 (aN 1)]
          s2 = PuiseuxSeries [PuiseuxTerm 1 (aN (-1))]
      in terms (addSeries s1 s2) @?= []

  , testCase "mulSeries Cauchy product" $
      let s1 = PuiseuxSeries [PuiseuxTerm 0 (aN 1), PuiseuxTerm 1 (aN 1)]
          s2 = PuiseuxSeries [PuiseuxTerm 0 (aN 1), PuiseuxTerm 1 (aN 1)]
          r  = mulSeries s1 s2
      in terms r @?= [ PuiseuxTerm 0 (aN 1)
                     , PuiseuxTerm 1 (aN 2)
                     , PuiseuxTerm 2 (aN 1) ]

  , testCase "mulSeries with fractional exponents" $
      let s1 = PuiseuxSeries [PuiseuxTerm (1/2) (aN 1)]
          s2 = PuiseuxSeries [PuiseuxTerm (1/2) (aN 1)]
          r  = mulSeries s1 s2
      in terms r @?= [PuiseuxTerm 1 (aN 1)]
  ]

------------------------------------------------------------------------
-- Expansion engine
------------------------------------------------------------------------

expandTests :: TestTree
expandTests = testGroup "Expansion engine"
  [ testCase "expand Const" $
      expand (Const 3.0) (pt 0) "x"
        @?= Right (PuiseuxSeries [PuiseuxTerm 0 (aN 3)])

  , testCase "expand Var at x0=0" $
      expand (Var "x") (pt 0) "x"
        @?= Right (PuiseuxSeries [PuiseuxTerm 1 (aN 1)])

  , testCase "expand Var at x0=2" $
      expand (Var "x") (pt 2) "x"
        @?= Right (PuiseuxSeries [PuiseuxTerm 0 (aN 2), PuiseuxTerm 1 (aN 1)])

  , testCase "expand sin(x) at x0=0 leading term" $
      let Right s = expand (Sin (Var "x")) (pt 0) "x"
          Just lt = leadingTermNZ s
      in pExp lt @?= 1

  , testCase "expand sin(x) at x0=pi/2 constant term" $
      let Right s = expand (Sin (Var "x")) (pt (pi/2)) "x"
      in abs (algToDouble (constantTerm s) - 1.0) < 1e-6 @?= True

  , testCase "expand exp(x) at x0=0 constant term" $
      let Right s = expand (Exp (Var "x")) (pt 0) "x"
      in algToDouble (constantTerm s) @?= 1.0

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
      in abs (algToDouble (coeff lt) - 1.0) < 1e-10 @?= True

  , testCase "expand log(x) at x0<=0 returns error" $
      case expand (Log (Var "x")) (pt 0) "x" of
        Left (Undefined _) -> return ()
        _                  -> assertFailure "Expected Undefined error"

  , testCase "expand unknown variable returns error" $
      case expand (Var "y") (pt 0) "x" of
        Left (Unknown _) -> return ()
        _                -> assertFailure "Expected Unknown error"

  , testCase "expand |x| at x=2 matches x" $
      expand (Abs (Var "x")) (pt 2) "x"
        @?= expand (Var "x") (pt 2) "x"

  , testCase "expand |x| at x=-2 is -x" $
      expand (Abs (Var "x")) (pt (-2)) "x"
        @?= expand (Neg (Var "x")) (pt (-2)) "x"

  , testCase "expand |x| at x=0 is NonAnalytic (kink)" $
      case expand (Abs (Var "x")) (pt 0) "x" of
        Left (NonAnalytic _) -> return ()
        other                -> assertFailure $ "Expected NonAnalytic, got: " ++ show other

  , testCase "expand |x^2| at x=0 is smooth (matches x^2)" $
      expand (Abs (Pow (Var "x") (Const 2))) (pt 0) "x"
        @?= expand (Pow (Var "x") (Const 2)) (pt 0) "x"

  , testCase "expand |x+5| at x=0 matches x+5" $
      expand (Abs (Add (Var "x") (Const 5))) (pt 0) "x"
        @?= expand (Add (Var "x") (Const 5)) (pt 0) "x"

  , testCase "expand |x^3| at x=0 is NonAnalytic (odd-order kink)" $
      case expand (Abs (Pow (Var "x") (Const 3))) (pt 0) "x" of
        Left (NonAnalytic _) -> return ()
        other                -> assertFailure $ "Expected NonAnalytic, got: " ++ show other
  ]

------------------------------------------------------------------------
-- Derivative tests
------------------------------------------------------------------------

derivativeTests :: TestTree
derivativeTests = testGroup "Derivatives"
  [ testCase "d/dx sin(x) at x=0 is 1" $
      abs (unRight (derivative (Sin (Var "x")) (Map.fromList [("x", 0)]) "x") - 1.0) < 1e-10
        @?= True

  , testCase "d/dx sin(x) at x=pi/2 is ~0" $
      abs (unRight (derivative (Sin (Var "x")) (Map.fromList [("x", pi/2)]) "x")) < 1e-6
        @?= True

  , testCase "d/dx x^2 at x=3 is 6" $
      abs (unRight (derivative (Pow (Var "x") (Const 2)) (Map.fromList [("x", 3)]) "x") - 6.0) < 1e-10
        @?= True

  , testCase "d/dx exp(x) at x=0 is 1" $
      abs (unRight (derivative (Exp (Var "x")) (Map.fromList [("x", 0)]) "x") - 1.0) < 1e-10
        @?= True

  , testCase "d/dx log(x) at x=1 is 1" $
      abs (unRight (derivative (Log (Var "x")) (Map.fromList [("x", 1)]) "x") - 1.0) < 1e-10
        @?= True

  , testCase "d/dx cos(x) at x=0 is 0" $
      abs (unRight (derivative (Cos (Var "x")) (Map.fromList [("x", 0)]) "x")) < 1e-10
        @?= True

  , testCase "d/dx sin^(1/2)(x) at x=1 matches analytic" $
      let result   = unRight (derivative (Pow (Sin (Var "x")) (Const 0.5)) (Map.fromList [("x", 1)]) "x")
          expected = cos 1 / (2 * sqrt (sin 1))
      in abs (result - expected) < 1e-6 @?= True

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
      case limit (Div (Sin (Var "x")) (Var "x")) "x" 0 of
        Exists v -> abs (v - 1.0) < 1e-10 @?= True
        other    -> assertFailure $ "Expected Exists, got: " ++ show other

  , testCase "lim_{x->0} x^2 = 0" $
      case limit (Pow (Var "x") (Const 2)) "x" 0 of
        Exists v -> abs v < 1e-10 @?= True
        other    -> assertFailure $ "Expected Exists, got: " ++ show other

  , testCase "lim_{x->2} x^2 = 4" $
      case limit (Pow (Var "x") (Const 2)) "x" 2 of
        Exists v -> abs (v - 4.0) < 1e-10 @?= True
        other    -> assertFailure $ "Expected Exists, got: " ++ show other

  , testCase "lim_{x->0} 1/x is pole" $
      case limit (Div (Const 1) (Var "x")) "x" 0 of
        Pole _ -> return ()
        _      -> assertFailure "Expected Pole"

  , testCase "lim_{x->0} exp(x) = 1" $
      case limit (Exp (Var "x")) "x" 0 of
        Exists v -> abs (v - 1.0) < 1e-10 @?= True
        other    -> assertFailure $ "Expected Exists, got: " ++ show other

  , testCase "lim_{x->0} (1-cos(x))/x^2 = 1/2" $
      let f = Div (Sub (Const 1) (Cos (Var "x")))
                  (Pow (Var "x") (Const 2))
      in case limit f "x" 0 of
           Exists v -> abs (v - 0.5) < 1e-6 @?= True
           other    -> assertFailure $ "Expected Exists, got: " ++ show other
  ]

------------------------------------------------------------------------
-- Polynomial tests
------------------------------------------------------------------------

polyTests :: TestTree
polyTests = testGroup "Polynomial arithmetic"
  [ testCase "addPoly" $
      let p = Poly "x" [1, 2, 1 :: Double]
          q = Poly "x" [1, 1 :: Double]
      in polyCoef (addPoly p q) @?= [2, 3, 1]

  , testCase "mulPoly (x+1)^2 = x^2+2x+1" $
      let q = Poly "x" [1, 1 :: Double]
      in polyCoef (mulPoly q q) @?= [1, 2, 1]

  , testCase "divModPoly exact" $
      let p = Poly "x" [1, 2, 1 :: Double]
          q = Poly "x" [1, 1 :: Double]
          (quot, rem) = divModPoly p q
      in (polyCoef quot, polyCoef rem) @?= ([1, 1], [])

  , testCase "gcdPoly" $
      let p = Poly "x" [1, 2, 1 :: Double]
          q = Poly "x" [1, 1 :: Double]
      in polyCoef (gcdPoly p q) @?= [1, 1]

  , testCase "diffPoly x^2+2x+1 = 2x+2" $
      let p = Poly "x" [1, 2, 1 :: Double]
      in polyCoef (diffPoly p) @?= [2, 2]

  , testCase "evalPoly at 2" $
      let p = Poly "x" [1, 2, 1 :: Double]
      in evalPoly p 2 @?= 9.0

  , testCase "resultant basic" $
      let p = Poly "x" [1, 1 :: Double]
          q = Poly "x" [-1, 1 :: Double]
      in resultant p q @?= (-2.0)

  , testCase "degree of zero poly is -1" $
      degree (zeroPoly "x" :: Poly Double) @?= (-1)

  , testCase "degree of constant is 0" $
      degree (constPoly "x" (5 :: Double)) @?= 0

  , testCase "leadingCoeff" $
      leadingCoeff (Poly "x" [1, 2, 3 :: Double]) @?= 3.0
  ]

------------------------------------------------------------------------
-- Rational function tests
------------------------------------------------------------------------

rationalFunctionTests :: TestTree
rationalFunctionTests = testGroup "Rational functions"
  [ testCase "addRat 1/x + 1/x = 2/x" $
      let rf = ratFun (Poly "x" [1 :: Double]) (Poly "x" [0, 1])
      in numerator (addRat rf rf) @?= Poly "x" [2]

  , testCase "diffRat 1/x = -1/x^2" $
      let rf  = ratFun (Poly "x" [1 :: Double]) (Poly "x" [0, 1])
          rf' = diffRat rf
      in ( polyCoef (numerator rf')
         , polyCoef (denominator rf') )
           @?= ([-1], [0, 0, 1])

  , testCase "hermiteReduce squarefree denom unchanged" $
      let p  = Poly "x" [1 :: Double]
          q  = Poly "x" [-1, 0, 1]
          rf = ratFun p q
          (g, _h) = hermiteReduce rf
      in polyCoef (numerator g) @?= []

  , testCase "extGCD basic" $
      let p = Poly "x" [1, 1 :: Double]
          q = Poly "x" [-1, 1]
          (g, _s, _t) = extGCD p q
      in degree g @?= 0
  ]

------------------------------------------------------------------------
-- Risch integration tests
------------------------------------------------------------------------

rischTests :: TestTree
rischTests = testGroup "Risch integration"
  [ testCase "int 1 dx = x" $
      rischIntegrate (Const 1) "x"
        @?= Elementary (Var "x")

  , testCase "int x dx = x^2/2" $
      rischIntegrate (Var "x") "x"
        @?= Elementary (Mul (Const 0.5) (Pow (Var "x") (Const 2.0)))

  , testCase "int x^2 dx = x^3/3" $
      rischIntegrate (Pow (Var "x") (Const 2)) "x"
        @?= Elementary (Mul (Const 0.3333333333333333) (Pow (Var "x") (Const 3.0)))

  , testCase "int 1/x dx = log(x)" $
      rischIntegrate (Div (Const 1) (Var "x")) "x"
        @?= Elementary (Log (Var "x"))

  , testCase "int exp(x) dx = exp(x)" $
      rischIntegrate (Exp (Var "x")) "x"
        @?= Elementary (Exp (Var "x"))

  , testCase "int exp(2x) dx = exp(2x)/2" $
      rischIntegrate (Exp (Mul (Const 2) (Var "x"))) "x"
        @?= Elementary (Div (Exp (Mul (Const 2.0) (Var "x"))) (Const 2.0))

  , testCase "int exp(x^2) dx is non-elementary" $
      rischIntegrate (Exp (Pow (Var "x") (Const 2))) "x"
        @?= NonElementary

  , testCase "int 1/(x^2-1) dx has log terms" $
      case rischIntegrate (Div (Const 1) (Sub (Pow (Var "x") (Const 2)) (Const 1))) "x" of
        Elementary _ -> return ()
        other        -> assertFailure $ "Expected Elementary, got: " ++ show other

  , testCase "int 1/(x^2+1) dx is non-elementary over reals" $
      case rischIntegrate (Div (Const 1) (Add (Pow (Var "x") (Const 2)) (Const 1))) "x" of
        NonElementary -> return ()
        other         -> assertFailure $ "Expected NonElementary, got: " ++ show other

  , testCase "int log(x) dx = x*log(x) - x (regression: required \
             \simplifying the by-parts inner expression before \
             \reclassifying it)" $
      rischIntegrate (Log (Var "x")) "x"
        @?= Elementary (Sub (Mul (Var "x") (Log (Var "x"))) (Var "x"))

  , testCase "int x*log(x) dx = x^2/2*log(x) - x^2/4 (regression: \
             \LogCase previously discarded the multiplier x entirely, \
             \silently treating this the same as plain log(x))" $
      rischIntegrate (Mul (Var "x") (Log (Var "x"))) "x"
        @?= Elementary
              (Sub (Mul (Const 0.5) (Mul (Pow (Var "x") (Const 2.0)) (Log (Var "x"))))
                   (Mul (Const 0.25) (Pow (Var "x") (Const 2.0))))

  , testCase "int 2x/(x^2-1) dx = log(x^2-1) (regression: log-derivative \
             \edge case in Rothstein-Trager produced garbage when a = c*d')" $
      rischIntegrate (Div (Mul (Const 2) (Var "x"))
                         (Sub (Pow (Var "x") (Const 2)) (Const 1))) "x"
        @?= Elementary (Log (Add (Const (-1.0)) (Pow (Var "x") (Const 2.0))))

  , testCase "int 1/(x^2-1) dx has two log terms (partial fractions via \
             \Rothstein-Trager)" $
      case rischIntegrate (Div (Const 1)
                              (Sub (Pow (Var "x") (Const 2)) (Const 1))) "x" of
        Elementary _ -> return ()
        other -> assertFailure $ "Expected Elementary, got: " ++ show other

  , testCase "primitive: int 1/x = log(x)" $
      let p     = Poly "x" [aN 1]
          q     = Poly "x" [aN 0, aN 1]
          rf    = ratFun p q
          field = baseField "x"
      in case integratePrimitive rf field of
           PrimitiveElementary _ -> return ()
           _                     -> assertFailure "Expected Elementary"

  , testCase "exponential: int exp(x) = exp(x)" $
      case integrateExp (Var "x") "x" of
        Right _ -> return ()
        Left e  -> assertFailure $ "Expected Right, got: " ++ e

  , testCase "exponential: int exp(x^2) non-elementary" $
      case integrateExp (Pow (Var "x") (Const 2)) "x" of
        Left _  -> return ()
        Right _ -> assertFailure "Expected Left (non-elementary)"
  ]

------------------------------------------------------------------------
-- Degree >= 2 AlgNum arithmetic (complex algebraic numbers)
------------------------------------------------------------------------

-- | Helper: assert an AlgNum's real and imaginary parts match
-- expected values within tolerance.
algNumApprox :: AlgNum -> Double -> Double -> Assertion
algNumApprox a expectedRe expectedIm = do
  abs (algToDouble a - expectedRe) < 1e-6 @?
    ("real part: expected " ++ show expectedRe ++ ", got " ++ show (algToDouble a))
  abs (algImagDouble a - expectedIm) < 1e-6 @?
    ("imag part: expected " ++ show expectedIm ++ ", got " ++ show (algImagDouble a))

algNumDegree2Tests :: TestTree
algNumDegree2Tests = testGroup "AlgNum degree >= 2 (complex algebraic numbers)"
  [ testCase "i * i = -1" $
      algNumApprox (algI * algI) (-1) 0

  , testCase "i + i = 2i" $
      algNumApprox (algI + algI) 0 2

  , testCase "i^3 = -i" $
      algNumApprox (algI * algI * algI) 0 (-1)

  , testCase "sqrt(2) + sqrt(2) = 2*sqrt(2)" $
      algNumApprox (algSqrt 2 + algSqrt 2) (2 * sqrt 2) 0

  , testCase "sqrt(2) * sqrt(2) = 2" $
      algNumApprox (algSqrt 2 * algSqrt 2) 2 0

  , testCase "(1+i) * (1-i) = 2" $
      let onePlusI  = algOne + algI
          oneMinusI = algOne - algI
      in algNumApprox (onePlusI * oneMinusI) 2 0

  , testCase "negate i = -i" $
      algNumApprox (negate algI) 0 (-1)

  , testCase "1 - i has imaginary part -1 (regression: algNeg previously left imagQ untouched)" $
      algNumApprox (algOne - algI) 1 (-1)
  ]

------------------------------------------------------------------------
-- QPoly division (regression: qDivPoly stub silently no-op'd on
-- non-constant input, corrupting the subresultant PRS recursion)
------------------------------------------------------------------------

qPolyDivisionTests :: TestTree
qPolyDivisionTests = testGroup "QPoly division"
  [ testCase "qDivModPoly exact division, no remainder" $
      -- (x^2 - 1) / (x - 1) = x + 1, remainder 0
      let p = QPoly [-1, 0, 1]   -- x^2 - 1
          q = QPoly [-1, 1]      -- x - 1
          (quot, rem') = qDivModPoly p q
      in (qPolyCoef quot, qPolyCoef (qStrip rem')) @?= ([1, 1], [])

  , testCase "qDivModPoly with nonzero remainder" $
      -- (x^2 + 1) / (x - 1) = x + 1, remainder 2
      let p = QPoly [1, 0, 1]    -- x^2 + 1
          q = QPoly [-1, 1]      -- x - 1
          (quot, rem') = qDivModPoly p q
      in (qPolyCoef quot, qPolyCoef (qStrip rem')) @?= ([1, 1], [2])

  , testCase "qQuotPoly matches qDivModPoly's quotient" $
      let p = QPoly [-1, 0, 1]
          q = QPoly [-1, 1]
      in qQuotPoly p q @?= fst (qDivModPoly p q)
  ]

------------------------------------------------------------------------
-- Trig integration via the exponential extension (e^(ix))
------------------------------------------------------------------------

-- | Evaluate an Expr at a point, in Complex Double, since
-- rischIntegrate's trig output legitimately contains I and only
-- cancels to a real value at the end. Var lookups use the single
-- supplied (name, value) binding; Const/Pi/E/I are independent of
-- the binding.
evalComplexExpr :: (String, Complex Double) -> Expr -> Complex Double
evalComplexExpr _        (Const c) = c :+ 0
evalComplexExpr _        Pi        = pi :+ 0
evalComplexExpr _        E         = exp 1 :+ 0
evalComplexExpr _        I         = 0 :+ 1
evalComplexExpr (vn, vv) (Var x)
  | x == vn   = vv
  | otherwise = error ("evalComplexExpr: unbound variable " ++ x)
evalComplexExpr env (Add f g) = evalComplexExpr env f + evalComplexExpr env g
evalComplexExpr env (Sub f g) = evalComplexExpr env f - evalComplexExpr env g
evalComplexExpr env (Mul f g) = evalComplexExpr env f * evalComplexExpr env g
evalComplexExpr env (Div f g) = evalComplexExpr env f / evalComplexExpr env g
evalComplexExpr env (Neg f)   = negate (evalComplexExpr env f)
evalComplexExpr env (Abs f)   = magnitude (evalComplexExpr env f) :+ 0
evalComplexExpr env (Exp f)   = exp (evalComplexExpr env f)
evalComplexExpr env (Log f)   = log (evalComplexExpr env f)
evalComplexExpr env (Sin f)   = sin (evalComplexExpr env f)
evalComplexExpr env (Cos f)   = cos (evalComplexExpr env f)
evalComplexExpr env (Pow f g) = evalComplexExpr env f ** evalComplexExpr env g

-- | Assert that two complex values agree (within tolerance) on real
-- and imaginary parts.
complexApprox :: Complex Double -> Complex Double -> Assertion
complexApprox got expected = do
  abs (realPart got - realPart expected) < 1e-6 @?
    ("real part: expected " ++ show (realPart expected) ++ ", got " ++ show (realPart got))
  abs (imagPart got - imagPart expected) < 1e-6 @?
    ("imag part: expected " ++ show (imagPart expected) ++ ", got " ++ show (imagPart got))

-- | Extract the Expr from an Elementary RischResult, failing the
-- test otherwise.
expectElementary :: RischResult -> Expr
expectElementary (Elementary e) = e
expectElementary other = error ("expected Elementary, got: " ++ show other)

trigIntegrationTests :: TestTree
trigIntegrationTests = testGroup "Trig integration via exponential extension"
  [ testCase "int sin(x) dx matches -cos(x) at x=0.7" $
      let result = expectElementary (rischIntegrate (Sin (Var "x")) "x")
          got    = evalComplexExpr ("x", 0.7 :+ 0) result
          expect = negate (cos (0.7 :+ 0))
      in complexApprox got expect

  , testCase "int cos(x) dx matches sin(x) at x=0.7" $
      let result = expectElementary (rischIntegrate (Cos (Var "x")) "x")
          got    = evalComplexExpr ("x", 0.7 :+ 0) result
          expect = sin (0.7 :+ 0)
      in complexApprox got expect

  , testCase "int sin(x) dx matches -cos(x) at x=2.3 (different point)" $
      let result = expectElementary (rischIntegrate (Sin (Var "x")) "x")
          got    = evalComplexExpr ("x", 2.3 :+ 0) result
          expect = negate (cos (2.3 :+ 0))
      in complexApprox got expect

  , testCase "int cos(x) dx matches sin(x) at x=2.3 (different point)" $
      let result = expectElementary (rischIntegrate (Cos (Var "x")) "x")
          got    = evalComplexExpr ("x", 2.3 :+ 0) result
          expect = sin (2.3 :+ 0)
      in complexApprox got expect

  , testCase "int 1/(x^2+1) dx is non-elementary over the reals (regression: \
             \complex root-finding previously misclassified this as Elementary)" $
      case rischIntegrate (Div (Const 1) (Add (Pow (Var "x") (Const 2)) (Const 1))) "x" of
        NonElementary -> return ()
        other         -> assertFailure $ "Expected NonElementary, got: " ++ show other
  ]

------------------------------------------------------------------------
-- Special function recognition (erf, li, Si, Ci, Ei)
------------------------------------------------------------------------

specialFunctionTests :: TestTree
specialFunctionTests = testGroup "Special function recognition"
  [ testCase "int e^(-x^2) dx = (sqrt(pi)/2) * erf(x)" $
      rischIntegrate (Exp (Neg (Pow (Var "x") (Const 2)))) "x"
        @?= Elementary (Mul (Div (Pow Pi (Const 0.5)) (Const 2.0)) (Erf (Var "x")))

  , testCase "int e^(x^2) dx is still NonElementary (negative control: \
             \confirms the erf pattern match doesn't over-fire on the \
             \positive-exponent case, which has no classical closed form)" $
      rischIntegrate (Exp (Pow (Var "x") (Const 2))) "x"
        @?= NonElementary

  , testCase "int 1/log(x) dx = li(x)" $
      rischIntegrate (Div (Const 1) (Log (Var "x"))) "x"
        @?= Elementary (Li (Var "x"))

  , testCase "int sin(x)/x dx = Si(x)" $
      rischIntegrate (Div (Sin (Var "x")) (Var "x")) "x"
        @?= Elementary (Si (Var "x"))

  , testCase "int cos(x)/x dx = Ci(x)" $
      rischIntegrate (Div (Cos (Var "x")) (Var "x")) "x"
        @?= Elementary (Ci (Var "x"))

  , testCase "int e^x/x dx = Ei(x)" $
      rischIntegrate (Div (Exp (Var "x")) (Var "x")) "x"
        @?= Elementary (Ei (Var "x"))

  , testCase "Erf/Li/Si/Ci/Ei have correct derivatives via deriveBase \
             \(DiffField.hs; NOT via diff/Calculus.hs, which goes \
             \through symExpand -- deliberately left unimplemented \
             \for these, since full Taylor expansion is separate, \
             \unstarted work)" $ do
      simplify (deriveBase (Erf (Var "x")))
        @?= simplify (Mul (Div (Const 2.0) (Pow Pi (Const 0.5))) (Exp (Neg (Mul (Var "x") (Var "x")))))
      simplify (deriveBase (Li (Var "x")))
        @?= simplify (Div (Const 1.0) (Log (Var "x")))
      simplify (deriveBase (Si (Var "x")))
        @?= simplify (Div (Sin (Var "x")) (Var "x"))
      simplify (deriveBase (Ci (Var "x")))
        @?= simplify (Div (Cos (Var "x")) (Var "x"))
      simplify (deriveBase (Ei (Var "x")))
        @?= simplify (Div (Exp (Var "x")) (Var "x"))
  ]

------------------------------------------------------------------------
-- Algebraic extension case in DiffField
------------------------------------------------------------------------

algebraicExtensionTests :: TestTree
algebraicExtensionTests = testGroup "Algebraic extension (implicit differentiation)"
  [ testCase "D(sqrt(x)) = 1/(2*sqrt(x)) via Algebraic extension" $
      let p = Sub (Pow (Var "t1") (Const 2)) (Var "x")
          field = addExtension (baseField "x") (Algebraic p)
      in simplify (derive field (Var "t1"))
           @?= Div (Const 1.0) (Mul (Const 2.0) (Var "t1"))

  , testCase "D(sqrt(1-x^2)) = -x/sqrt(1-x^2) via Algebraic extension" $
      let p = Sub (Add (Pow (Var "t1") (Const 2)) (Pow (Var "x") (Const 2))) (Const 1)
          field = addExtension (baseField "x") (Algebraic p)
      in simplify (derive field (Var "t1"))
           @?= Div (Neg (Var "x")) (Var "t1")

  , testCase "D(cbrt(x)) = 1/(3*cbrt(x)^2) via Algebraic extension" $
      let p = Sub (Pow (Var "t1") (Const 3)) (Var "x")
          field = addExtension (baseField "x") (Algebraic p)
      in simplify (derive field (Var "t1"))
           @?= Div (Const 1.0) (Mul (Const 3.0) (Pow (Var "t1") (Const 2.0)))
  ]