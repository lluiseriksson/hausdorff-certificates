# Mathematical statements used by this repository

Notation: `mu` is a finite positive Borel measure, `b_n = int v^n dmu(v)`,
`H_N = (b_{i+j})_{0<=i,j<=N}`, `L_N = (b_{i+j} - b_{i+j+1})`,
`L_N^theta = (theta b_{i+j} - b_{i+j+1})`, `S_N = (b_{i+j+1})`.

## 1. Hausdorff criterion and the finite certificates

**Theorem (Hausdorff).** A real sequence `(b_n)_{n>=0}` is the moment
sequence of a positive measure on `[0,1]` iff it is completely monotone:
`(-1)^k Delta^k b_n >= 0` for all `n, k >= 0`, where `Delta b_n = b_{n+1} -
b_n`. Equivalently, iff `H_N >= 0` and `L_N >= 0` for every `N`.

*Why the matrices:* for any polynomial `p(v) = sum c_i v^i`,
`c^T H_N c = int p(v)^2 dmu >= 0` and `c^T L_N c = int p(v)^2 (1-v) dmu >= 0`;
conversely positivity of all `H_N, L_N` solves the moment problem on `[0,1]`
(classical; see e.g. Shohat--Tamarkin, *The Problem of Moments*).

The difference table `(-1)^k Delta^k b_n = int v^n (1-v)^k dmu` gives the
same information localised at `(n, k)`; the implementation uses it as a
secondary, cheaply localised falsifier (`cm_check` in the certificates).

**What a finite run means.** `H_N, L_N >= 0` for the tested `N` is
consistency only. A certified negative eigenvalue at any single `N` refutes
the existence of the representing measure -- that is the falsifier this
repository is built around.

## 2. Support test (the theta certificate)

**Proposition.** Let `mu >= 0` on `[0,1]` and `theta in (0,1]`. Then
`supp(mu) subset [0, theta]` iff `H_N >= 0` and `L_N^theta >= 0` for all `N`.

*Proof.* If the support bound holds, `c^T L_N^theta c = int p(v)^2 (theta -
v) dmu >= 0`. Conversely, positivity of `(b_n)` and of `(theta b_n -
b_{n+1})` says both `(b_n)` and `(theta b_n - b_{n+1})` are Stieltjes-type
positive sequences; rescaling `w = v/theta` turns `(b_n theta^{-n})` into a
Hausdorff sequence on `[0,1]`, whose representing measure pushed back has
support in `[0, theta]`. (Diagonal congruence: `L_N^theta` for `b` equals
`D L_N D` for the rescaled sequence with `D = diag(theta^i)` up to a
positive factor, so the two formulations certify each other.) QED

## 3. Yang--Mills-side dictionary: coercivity as a support bound

Let `K = K^T` be finite-dimensional with eigenvalues `lambda_j`, fix
`x0 > 0`, and set

```
M_n(K; x0) = Tr (K + x0 I)^{-(n+1)} = sum_j (lambda_j + x0)^{-(n+1)},
b_n        = x0^n M_n = sum_j u_j v_j^n,
u_j = 1/(lambda_j + x0),   v_j = x0/(lambda_j + x0).
```

So `b_n = int v^n dmu` with `mu = sum u_j delta_{v_j}`, and

```
K >= c I   <=>   lambda_min >= c   <=>   v_max <= theta := x0/(c + x0)
           <=>   supp(mu) subset (0, theta].
```

Hence `H_N >= 0` and `L_N^theta >= 0` (all `N`) certify the claimed
coercivity constant; a certified negative eigenvalue of `L_N^theta` at any
`N` refutes it. Since `v_max` determines `lambda_min = x0(1 - v_max)/v_max`,
the hierarchy also *estimates* the coercivity constant from below. All of
this is computed exactly over `Q` (`moments.resolvent_trace_moments`).

This is the "cheap finite screen" for candidate precisions/quadratic forms
on the P4 front: refute wrong constants for pennies before attempting an
analytic proof. It says nothing about uniformity in volume or scale.

## 4. Riemann-side dictionary: single-point resolvent-Hausdorff criterion

For `x > 1/4` define `S_Xi(x) = (1/(2 sqrt x)) (xi'/xi)(1/2 + sqrt x)`.
Under RH, `S_Xi(x) = sum_{gamma > 0} 1/(gamma^2 + x)` (ordinates of the
nontrivial zeros), which is the Stieltjes transform of the positive measure
`nu = sum delta_{gamma^2}`. Fix one `x0 > 1/4` and set

```
b_n(x0) = x0^n (-1)^n S_Xi^{(n)}(x0) / n!
        = int_0^1 v^n dmu_{x0}(v),   mu_{x0} = sum u_gamma delta_{v_gamma},
u_gamma = 1/(gamma^2 + x0),   v_gamma = x0/(gamma^2 + x0),
```

so under RH `(b_n)` is a Hausdorff moment sequence and every `H_N, L_N >= 0`.
Conversely (source criterion, knowledge-tree card 8): complete monotonicity
of `(b_n(x0))` at a *single* `x0` yields, via the Hausdorff theorem, a
holomorphic Stieltjes extension of `S_Xi` to the slit plane
`C \\ (-infty, 0]`, which excludes zeros of `Xi` off the real axis and hence
implies RH. **Novelty status of this equivalence is an open verification
gate** (knowledge-tree `verification_queue` V3): moment-problem criteria for
RH already exist in the literature; only the packaging (one-point resolvent
derivatives + slit-plane extension + reuse of the same certificates on the
YM side) is a candidate for novelty. Nothing in this repository depends on
the novelty claim.

An interval-certified negative eigenvalue of `H_N` or `L_N` built from
correct enclosures of `b_n(x0)` would refute RH. Finitely many PSD
verdicts prove nothing.

## 5. Enclosures for the zeta pipeline

### 5.1 Truncated part

The first `Z` ordinates are taken from `mpmath.zetazero` at 60 dps and
inflated to intervals `gamma_k +- 10^{-52}`; all subsequent arithmetic is
outward-rounded interval arithmetic at >= 220 dps. (Accuracy of
`mpmath.zetazero` within the inflation radius is verification gate V-B in
SOURCES.md.)

### 5.2 Tail bound

Let `T` be a certified lower bound for `gamma_Z`, `f(t) = (t^2 +
x0)^{-(n+1)}` (positive, decreasing, `f(infty) = 0`), and `N(t)` the zero
counting function, so `N(T) >= Z`. Stieltjes integration by parts twice:

```
sum_{gamma > T} f(gamma) = -f(T) N(T) + int_T^infty (-f'(t)) N(t) dt
    <= f(T) (N_up(T) - Z) + int_T^infty f(t) N_up'(t) dt
```

with any `C^1` upper bound `N_up >= N`. We use Trudgian (2014):

```
|N(t) - (t/2pi) log(t/(2pi e)) - 7/8| <= C1 log t + C2 log log t + C3,
(C1, C2, C3) = (0.112, 0.278, 2.510),   t >= e,
```

(verification gate V-A) so `N_up'(t) = (1/2pi) log(t/2pi) + C1/t +
C2/(t log t)`. Bounding `f(t) <= t^{-m}`, `m = 2(n+1)`, the integral has
elementary closed-form upper bounds:

```
int_T^infty t^{-m} log(t/c) dt = T^{1-m} [ log(T/c)/(m-1) + 1/(m-1)^2 ],
int_T^infty C1 t^{-m-1} dt     = C1 T^{-m}/m,
int_T^infty C2 t^{-m-1}/log t dt <= (C2/log T) T^{-m}/m,
```

all evaluated in interval arithmetic (`zeta.tail_upper_bound`). The final
enclosure is `[trunc.lo, (trunc + [0, tail]).hi]`.

### 5.3 Lemma S (monotone splitting) -- why consistency does not pay the tail

**Lemma S.** Let `mu = mu_1 + mu_2` with `mu_2 >= 0` supported in `[0, s]`,
`s <= 1`. Then for every `N`:

```
H_N(mu) >= H_N(mu_1),    L_N(mu) >= L_N(mu_1),
L_N^theta(mu) >= L_N^theta(mu_1)  whenever theta >= s.
```

*Proof.* Each matrix is linear in the measure, and the `mu_2` block is
itself a moment matrix of a positive measure against the nonnegative
weights `1`, `(1 - v)`, `(theta - v)` respectively on `[0, s]`, hence PSD;
a PSD summand preserves the Loewner order. QED

Consequence: PSD certificates computed on the tight truncated-measure
enclosures (`include_tail=False`, applied with `mu_1 = mu_trunc`, `mu_2 =
mu_tail`, `s = x0/(T^2+x0)`) transfer to the full zero measure. The tail
bound is only needed for the *refutation* direction, where the full
enclosure must contain the true moments. The reference artifacts implement
exactly this split; `zeta_x0_1_full__H_N4` documents why the refutation
enclosure is width-limited at small `n` (`width(b_0) ~ log T /(2 pi T)`,
about `4.3e-3` with 60 zeros) and therefore shallow until many more zeros
are added.

### 5.4 Cross-check

`zeta.s_xi_moments_via_derivatives` recomputes `b_n(x0)` from Cauchy-integral
derivatives of `S_Xi` through `xi'/xi(s) = 1/s + 1/(s-1) - (1/2) log pi +
(1/2) psi(s/2) + zeta'(s)/zeta(s)` -- floating quadrature, **no certificate
semantics**. That these values land inside the rigorous zero-sum + tail
enclosures for every tested `n` (artifact
`zeta_x0_1_derivative_crosscheck.json`, `all_inside = true`) simultaneously
sanity-checks the zero list, the interval plumbing and the direction of the
tail bound.

## 6. Certification semantics

* Exact backend: LDL^T with symmetric diagonal pivoting is a complete PSD
  decision over `Q`; a PSD matrix whose running Schur complement has a zero
  diagonal must have the whole corresponding row vanish, otherwise a 2x2
  indefinite block exists, and the Schur-complement witness `w` lifts to the
  original matrix via `v = (-A11^{-1} A12 w, w)` with `v^T A v = w^T S w`
  (implemented as one back-substitution and *re-verified against the
  original matrix* before the certificate is emitted).
* Interval backend: completed interval Cholesky with positive pivot lower
  bounds implies every symmetric member of the enclosure is PD (Rump-style
  verification); a rational `v` with `sup v^T [A] v < 0` implies every
  member fails PSD. Overlap of neither kind is reported `INCONCLUSIVE`.
