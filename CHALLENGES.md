# Technical Challenges in Building limcalc

This document records the significant technical challenges encountered during
development of limcalc, intended as reference material when writing the
accompanying paper. Challenges are grouped roughly by which session they
arose in, but the lessons are broadly applicable.

---

## 1. The Zero Detection Problem in AlgNum

The most persistent challenge in the AlgNum refactor. Algebraic numbers are
represented as minimal polynomial plus isolating rectangle, with arithmetic
via resultants. The problem: determining whether an AlgNum is zero requires
checking if the minimal polynomial's root (identified by the rectangle) is
actually zero. This is trivial for exact rational inputs but breaks down when:

- `algAdd (fromQ 1) (fromQ (-1))` produces `AlgNum(x, rect)` — the minimal
  polynomial `x` correctly identifies 0 as the root, but the isolating
  rectangle after arithmetic may not tightly contain 0
- `stripZeros` in `PuiseuxSeries` checks `abs (coeff t) > ε` but the
  `algToDouble` approximation via rectangle midpoint accumulates error
- The `|| True` placeholder in `hasRoot` — we deferred proper root isolation
  (winding number computation) and this caused `refineToRoot` to converge to
  wrong roots

**Lesson for the paper:** Exact algebraic number arithmetic requires careful
separation between the algebraic structure (minimal polynomial) and the
numerical approximation (isolating rectangle). The two can diverge under
compound arithmetic operations, requiring periodic re-isolation.

---

## 2. Infinite Loops in Series Computation

Several places where naive implementation caused non-termination:

- `binomialSeries` — `iterate (mulSeries w)` created an infinite list; without
  truncation at each step the series grew unboundedly. Fixed by truncating
  inside the iterate.
- `subresultantPRS` in BivPoly — the `qDivPoly` placeholder returning `p`
  unchanged for non-trivial cases caused the pseudo-division remainder to never
  shrink, looping forever.
- `resultant` in Poly — the original subresultant PRS implementation using
  integer exponentiation `^` crashed with "Negative exponent" for certain
  inputs. Required switching to `**` with `fromIntegral` throughout.
- `mulResultantQ` in BivPoly — fell through to `bivResultant` for the general
  case which uses the looping `subresultantPRS`. Required special-casing the
  linear polynomial case.

**Lesson for the paper:** Symbolic computation engines are particularly
susceptible to non-termination because the natural mathematical definition of
many operations (infinite series, recursive GCD) is infinite. Truncation and
depth bounds must be explicit architectural decisions, not afterthoughts.

---

## 3. The Puiseux Series Depth/Truncation Tradeoff

The `depth` parameter (`depth = 8`) is a global constant controlling how many
terms to compute. This creates tension:

- Too few terms: `sin^(1/2)(x)` requires enough terms in the sin expansion
  before taking the power, otherwise coefficients are wrong
- Too many terms: `binomialSeries` computing `w^n` for `n` up to `depth`
  makes the computation exponential in depth for composed functions
- Truncation must happen inside `iterate` not outside — discovered this the
  hard way when `binomialSeries` looped

The composition `sin^(1/2)(sin(x^2))` for example requires the engine to
handle nested compositions where each level potentially needs more terms from
the level below.

**Lesson for the paper:** Adaptive depth — computing more terms when
cancellation occurs in leading coefficients — is likely necessary for a
production system. The fixed depth is a pragmatic compromise.

---

## 4. The Circular Import Problem

AlgNum needed BivPoly for resultant computation, and BivPoly needed QPoly
types. The naive organization created a circular dependency. Resolution
required extracting `GaussianQ`, `QPoly`, and all their operations into a
separate `LimCalc.QPoly` module that neither AlgNum nor BivPoly imports from
each other.

**Lesson for the paper:** In a mathematical library with deep
interdependencies between algebraic structures, module architecture must be
designed around the mathematical dependency order — base types first, then
structures over those types, then algorithms using those structures. Circular
dependencies almost always indicate a missing abstraction layer.

---

## 5. The `exponent` Name Clash

Naming the exponent field of `PuiseuxTerm` as `exponent` caused an ambiguity
error — Haskell's Prelude exports a function called `exponent` for floating
point numbers. Required renaming to `pExp`. Similarly `properFraction` in
`RationalFunction` clashed with Haskell's `properFraction` from `RealFrac`.

**Lesson for the paper:** Mathematical terminology and Haskell Prelude names
frequently collide. A systematic naming convention (prefixes like `p` for
Puiseux, `q` for QPoly) reduces these collisions and improves code
readability.

---

## 6. The `Poly a` Parameterization Refactor

Originally `Poly` used `Double` coefficients. When we decided to parameterize
to `Poly a`, the ripple effect touched every downstream module —
`RationalFunction`, `DiffField`, `Risch.Primitive`, `Risch.Exponential`,
`Risch`, and the test suite. Key difficulties:

- `pseudoDivMod` required `Fractional a` but was typed with `Num a` — the
  type error only appeared when downstream code tried to use it
- `RatFun` becoming `RatFun a` required updating every pattern match and type
  signature in Risch modules
- The test suite used `Poly "x" [1, 2, 1]` which became ambiguous — needed
  explicit `:: Double` annotations throughout

**Lesson for the paper:** Parameterizing over coefficient fields is
mathematically natural (field extensions are central to Risch) but
architecturally expensive. Starting with a parameterized type from the
beginning would have saved significant refactoring. The Haskell type system
made the refactor tractable — every type error was a precise pointer to a
site needing updating.

---

## 7. Floating Point Accumulation in AlgNum via Double

The current AlgNum implementation uses Double approximation internally for
transcendental functions — `algSin`, `algCos`, `algExp`, `algLog` all go
through `fromRational . toRational . sin . algToDouble`. This means:

- `algExp algZero` correctly returns approximately `1.0`
- But `algToDouble (expTaylor algZero)` accumulates error across 8 Taylor terms
- The constant term of `expand (Exp (Var "x")) (pt 0) "x"` comes out as
  `1.000000238` not `1.0`

The root cause is that `fromRational . toRational` converts a Double to a
rational approximation that is exact in ℚ but not equal to the true value.
Composition of many such approximations drifts.

**Lesson for the paper:** The AlgNum layer requires genuinely exact
transcendental function evaluation for the Taylor series coefficients — in
practice this means keeping coefficients as exact rationals (for polynomial
and rational functions) and only approximating at the final evaluation step.
The current architecture conflates these two concerns.

---

## 8. The `hasRoot` Implementation for Isolating Rectangles

Proper root isolation in the complex plane requires the argument principle —
counting zeros by integrating the logarithmic derivative around the rectangle
boundary. This is the correct approach but requires complex contour
integration, which is itself a non-trivial computation. Our simplified
implementations:

- First attempt: `|| True` — always returned True, causing `refineToRoot` to
  always pick the upper half, converging to the wrong root
- Second attempt: sign change on boundary corners — missed roots on the real
  axis when all corner evaluations had the same sign
- Third attempt: sign change on real axis only — works for real roots but
  incorrect for complex roots

**Lesson for the paper:** The isolating rectangle representation for complex
algebraic numbers requires the argument principle for correctness.
Approximations based on boundary sign changes are insufficient in general.
This is a known difficulty in exact real arithmetic systems and motivates
alternative representations (Thom encoding, Cauchy index).

---

## 9. The Dropbox/OneDrive Build Directory Conflict

A practical engineering challenge: Haskell's cabal build performs atomic file
rename operations during the build, which conflict with cloud sync tools
(Dropbox, OneDrive) that hold file locks. The `--builddir` flag redirecting
build artifacts outside the sync directory resolved this, but required adding
the flag to every cabal invocation.

**Lesson for the paper:** Not relevant to the paper directly, but worth noting
in a README — build tools and cloud sync are frequently incompatible.

---

## 10. Symbolic Differentiation Without a Rules Engine

The key design challenge that motivated the entire project: how to compute
`d/dx f(x)` symbolically without pattern matching on the AST to apply
derivative rules.

The solution — expand `f(x + h)` as a `SymPuiseuxSeries` with `Expr`
coefficients, then read off the `h^1` coefficient — works elegantly but
required:

- A separate `SymPuiseuxSeries` type with `Expr` coefficients rather than
  `Double`
- Symbolic Taylor series for sin, cos, exp, log where coefficients are `Expr`
  nodes (`Sin(Var "x")`, `Cos(Var "x")` etc.) rather than numeric values
- Fixing `symCoeffAt` to sum all terms at a given exponent rather than
  returning the first — the product rule bug where `x·sin(x)` differentiated
  to only `x·cos(x)` missing the `sin(x)` term

**Lesson for the paper:** The cleanest formulation of the limit-first approach
uses two separate but parallel series types — `PuiseuxSeries AlgNum` for
numeric computation and `SymPuiseuxSeries` for symbolic differentiation. The
duality between these two mirrors the classical duality between numerical and
symbolic computation.

---

## 11. The `simplifyMul` Ping-Pong Infinite Loop

The "pull constant left" rules in `simplifyMul` — `x * (c * z) → c * (x * z)`
and `(c * y) * z → c * (y * z)` — created an infinite loop when both operands
had leading constants. `Const 3 * (Const 2 * x)` would rewrite to
`Const 2 * (Const 3 * x)` which rewrote back indefinitely, producing an
infinite `Expr` tree. The tree itself was finite at any given step but
`simplify`'s fixed-point check (`expr' == expr`) diverged because equality on
an infinite structure never terminates. Discovered when `partialDiff` started
calling `simplify` on `deriveExpr` output, which for the first time produced
expressions like `3*(2*y^1*1)`. Fixed by guarding both rules with `notConst`
— only pull a constant left when the other operand is not itself a constant.

**Lesson for the paper:** Fixed-point simplification loops are particularly
dangerous because the failure mode is non-termination rather than a wrong
answer. Rewriting rules that move terms between positions can oscillate
indefinitely. Guards based on structural properties of the operands (not just
the top-level constructor) are necessary to ensure termination.

---

## 12. The Rothstein-Trager Log-Derivative Edge Case

When the numerator `a` of a rational function `a/d` is a scalar multiple of
`d'` (the denominator's derivative), the Rothstein-Trager resultant polynomial
degenerates — the resultant's root at `z=c` causes `gcd(d, a - c*d') =
gcd(d, 0) = d`, returning the entire denominator as the factor rather than the
correct irreducible pieces. This silently produced garbage: `∫2x/(x²-1)dx`
returned a sum of `log(1)` terms (effectively zero) instead of `log(x²-1)`.
The fix is a pre-check (`logDerivativeCheck`) that detects `a = c*d'` before
entering the resultant machinery and returns `c*log(d)` directly.

**Lesson for the paper:** The Rothstein-Trager algorithm has a classical
degenerate case when the integrand is a logarithmic derivative. This case
arises naturally in practice (any rational function whose numerator equals a
scalar multiple of the denominator's derivative) and must be handled as a
pre-check rather than relying on the general algorithm to produce a meaningful
result.

---

## 13. The Durand-Kerner Precision Collapse in Multi-Factor Partial Fractions

Partial fractions with three or more linear factors required factoring the
denominator via Rothstein-Trager's resultant, whose roots are found
numerically by Durand-Kerner. For `1/((x-1)(x-2)(x-3))`, the resultant has a
repeated rational root at `z=0.5` (appearing twice), but Durand-Kerner
returned `0.5000000011...` and `0.5000000003...` — close but not identical.
When these imprecise values were used in `gcd(q, p - c*d')`, the GCD
computation degenerated because `p - 0.5000000011*d'` doesn't exactly vanish
at the right factors. All three GCDs returned the same wrong factor. Fixed by
`snapToRational`: after Durand-Kerner, snap each root to the nearest rational
with denominator ≤ 100. Since the Rothstein-Trager resultant has rational
coefficients when the input is rational, its rational roots are exactly
representable as simple fractions.

**Lesson for the paper:** Algorithms that are exact in theory but use
numerical root-finding internally require a rationalization step at the
boundary between numerical and symbolic computation. The proximity of
numerical roots to exact rationals can be exploited via continued fraction
approximation or bounded-denominator rational reconstruction.

---

## 14. The Circular Dependency Between RationalFunction and Risch.Primitive

`partialFractions` naturally belongs in `RationalFunction.hs` (it operates on
`RatFun` values), but its implementation requires `rothsteinTrager` from
`Risch.Primitive`, which itself imports `RationalFunction`. Moving
`partialFractions` to `Risch.Primitive` resolved the cycle but is
architecturally awkward — it places a presentation-layer function (partial
fraction decomposition) in an algorithmic module. The underlying cause is that
partial fraction decomposition over ℚ requires polynomial factorization, which
in our engine is accomplished via Rothstein-Trager as a side effect. A cleaner
resolution would factor polynomial factorization out into a separate module
that both `RationalFunction` and `Risch.Primitive` import.

**Lesson for the paper:** Partial fraction decomposition and Risch integration
are more tightly coupled than they initially appear — the same resultant/GCD
machinery that computes integration coefficients implicitly factors the
denominator polynomial. This coupling should be made explicit in the module
architecture rather than papering over it with a misplaced function.

---

## 15. The Algebraic Extension Derivation and the symExpand Inconsistency

`DiffField.deriveExpr` correctly implements symbolic differentiation via the
algebraic chain rule for all `Expr` constructors including `Erf`, `Si`, `Ei`,
`Ci`, `Arcsin`, etc. But `diff` in `Calculus.hs` goes through `symExpand`
(the series-expansion path), which for these constructors returned `Unknown`
or required implementing Taylor series with logarithmic terms. The result was
two parallel differentiation paths with different capabilities — `deriveBase`
worked correctly for all functions, `diff` didn't. Resolved by implementing
`erfSymTaylor`, `siSymTaylor`, `eiSymTaylor`, `ciSymTaylor` via iterated
`deriveBase`, making `diff` consistent with `deriveBase` for all implemented
functions. `Li` remains asymmetric — `deriveBase (Li x)` gives `1/log(x)`
correctly but `diff (Li x) "x"` returns `Unknown` because `li(x) = Ei(log(x))`
requires a doubly-logarithmic Puiseux extension.

**Lesson for the paper:** A system with two differentiation paths (algebraic
chain rule and series extraction) will inevitably diverge unless the series
path is systematically derived from the algebraic path. The resolution —
generating Taylor coefficients by iterating `deriveBase` — makes the series
path a consequence of the algebraic path rather than an independent
implementation, eliminating the inconsistency for all but the most singular
cases.

---

## 16. The Recursive Partial Fractions Decomposition Failure

The first two attempts at decomposing compound factors in `partialFractions`
both failed subtly. First attempt: recurse on `c*dᵢ'/dᵢ` (the
Rothstein-Trager log-derivative form) — but this is always itself a
log-derivative, so RT returns it as a single whole-denominator term, making
the recursion a no-op. Second attempt: call `partialFractions(1/dᵢ)` and
scale each result by `c` — wrong coefficients, because the sub-decomposition
of `1/dᵢ` is independent of the original numerator and the scaling is not a
simple multiplication. The correct solution: use RT only to *identify* the
irreducible factors (including recursive splitting of compound factors via
`partialFractions(1/dᵢ)`), then apply the residue formula
`A = p(r)/(q/(x-r))(r)` for each linear factor using the *original* `p` and
`q`. The factors come from the recursive structure; the coefficients come
entirely from the original problem.

**Lesson for the paper:** Partial fraction decomposition has two logically
separate subproblems — factoring the denominator and computing the numerator
coefficients. Conflating them (trying to extract coefficients from the same
Rothstein-Trager run that factors the denominator) produces incorrect results.
The residue formula for the coefficients must be applied to the original
rational function, not to the sub-problems generated during factoring.
