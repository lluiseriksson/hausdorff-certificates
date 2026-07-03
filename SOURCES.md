# Sources and verification gates

This repository follows the knowledge-tree discipline: any external claim a
certificate leans on is listed here with an explicit gate. Exact-backend
artifacts depend on **nothing** below; only the zeta pipeline
("demonstration tier") does.

## Gate V-A: Trudgian's explicit Riemann-von Mangoldt bound

Used in `zeta.tail_upper_bound`:

    |N(t) - (t/2pi) log(t/(2pi e)) - 7/8| <= 0.112 log t + 0.278 log log t + 2.510,  t >= e.

Reference: T. S. Trudgian, "An improved upper bound for the argument of the
Riemann zeta-function on the critical line II", J. Number Theory 134 (2014)
280-292. **Gate:** verify the three constants and the validity range
against the published paper before treating any zeta refutation as final.
(Sharper constants, e.g. Hasanalizade-Shen-Wong 2021, can be dropped into
`zeta.TRUDGIAN_C*` without touching anything else; artifacts must then be
regenerated.)

## Gate V-B: accuracy of mpmath's zetazero

`data/zeta_zeros_60_dps60.json` stores the first 60 ordinates computed by
`mpmath.zetazero` at 60 dps; enclosures inflate them by +-1e-52. **Gate:**
cross-check the stored ordinates against an independent table (e.g.
Odlyzko's tables of zeros) to at least 52 decimal places, or recompute at
higher precision, before treating any zeta refutation as final. `make
zeros` recomputes the cache; the file is committed so CI is deterministic.

## Gate V-C (upstream, informational): novelty of the single-point criterion

The RH <=> one-point-Hausdorff packaging is knowledge-tree
`verification_queue` item V3 (`batch-2026-07-03-ym-unblock`); this
repository does not depend on its novelty, only on the classical Hausdorff
theorem (Shohat-Tamarkin, *The Problem of Moments*, AMS 1943) and on the
elementary dictionaries proved in MATH.md.

## Classical background (no gate)

- Hausdorff moment problem / complete monotonicity: Shohat-Tamarkin (1943).
- Verified interval Cholesky: S. M. Rump, "Verification methods: rigorous
  results using floating-point arithmetic", Acta Numerica 19 (2010).
