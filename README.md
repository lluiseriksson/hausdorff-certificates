# hausdorff-certificates

Finite, independently checkable positive-semidefiniteness certificates for
Hausdorff moment sequences:

```
H_N  = ( b_{i+j} )              >= 0        (moment matrix of  d mu)
L_N  = ( b_{i+j} - b_{i+j+1} )  >= 0        (moment matrix of (1-v) d mu)
L_N^theta = ( theta b_{i+j} - b_{i+j+1} ) >= 0   (support-in-[0,theta] test)
```

Two backends:

* **exact-rational** -- pure `fractions.Fraction`; PSD is *decided*, and the
  certificate is either an exact `P^T A P = L D L^T` factorisation with
  `D >= 0`, or a rational witness `v` with `v^T A v < 0` (a finite
  falsifier).  No floating point anywhere.
* **interval** -- enclosure inputs `[lo, hi]` (mpmath); a completed interval
  Cholesky certifies that *every* matrix in the enclosure is positive
  definite; a rational witness whose interval Rayleigh quotient has negative
  upper bound certifies that *every* matrix in the enclosure fails PSD.
  Anything else is reported `INCONCLUSIVE` -- never silently rounded.

Certificates are deterministic JSON files (sorted keys, no timestamps) that
`python -m hausdorff_certificates.verify` re-validates **without trusting
the generator**. CI regenerates all reference artifacts and fails on any
byte-level drift.

## Scope contract (read this first)

This repository certifies **finite matrices only**.

* `PSD_CERTIFIED` for every tested `N` is a *consistency* statement. It can
  never prove the Riemann Hypothesis, coercivity of an infinite-volume
  operator, or anything of that kind.
* `NOT_PSD_CERTIFIED` is a genuine finite falsifier for the corresponding
  claim (a Hausdorff representation on `[0,1]`, a support bound, a claimed
  coercivity constant) -- for the zeta pipeline, *modulo* the two
  literature/computation gates listed in [SOURCES.md](SOURCES.md).
* Nothing here touches `hRpoly`, the M4/M5 milestones, or the Clay problem.
  Distance to any of those after running everything in this repo: unchanged.

## Position in the Eriksson programme

This is the implementation satellite of the
[physmath-knowledge-tree](https://github.com/lluiseriksson/physmath-knowledge-tree)
node `bridge.hausdorff_moment_certificates` (edges
`edge.hausdorff_certificates.riemann`, `edge.hausdorff_certificates.ym_rg_activity`;
batch `batch-2026-07-03-ym-unblock`, Bridge Card 8). Reproducible results
feed back to the tree's Reproducible Run Ledger. It is deliberately **not**
hosted in [THE-ERIKSSON-PROGRAMME](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME):
nothing in this repo enters the Lean core's dependency tree.

## Install and run

```bash
pip install -e ".[test]"     # only runtime dependency: mpmath
pytest                       # 43 tests incl. one slow zeta smoke test
make artifacts               # regenerate artifacts/ deterministically
make verify                  # independent re-validation of every certificate
python -m hausdorff_certificates.manifest artifacts/manifest.json  # artifact digest + hash check
```

## The two dictionaries

**Riemann side** (single-point resolvent-Hausdorff criterion; MATH.md 2).
For fixed `x0 > 1/4`, `b_n(x0) = x0^n (-1)^n S_Xi^(n)(x0) / n!` where
`S_Xi(x) = sum_{gamma>0} 1/(gamma^2 + x)`. Under RH this is a Hausdorff
moment sequence, so every `H_N, L_N >= 0`; an interval-certified negative
eigenvalue at any `N` would refute RH. Verifying finitely many `N` proves
nothing -- it is an exact, reproducible hierarchy sensitive to the whole
spectrum.

**Yang--Mills side** (finite coercivity tests; MATH.md 3).
For a candidate precision `K`, `b_n = x0^n Tr (K + x0 I)^{-(n+1)}` has the
same structure, and `K >= c I` is *equivalent* to `supp(mu) subset (0,
theta]` with `theta = x0/(c + x0)`. So `H_N >= 0` together with
`L_N^theta >= 0` is a cheap finite screen for claimed coercivity constants
before investing in an analytic proof, and a certified negative eigenvalue
refutes the claimed constant exactly.

## Reference artifacts (`artifacts/`)

| artifact | backend | verdict | what it demonstrates |
|---|---|---|---|
| `hilbert_lebesgue__{H,L}` | exact | PSD | Lebesgue moments; H is the Hilbert matrix |
| `two_atoms__{H,L}` | exact | PSD | atomic measure, Hausdorff pair |
| `two_atoms_support_true` | exact | PSD | true support claim theta = 3/4 |
| `two_atoms_support_false` | exact | NOT_PSD | false claim theta = 1/2; exact witness |
| `laplacian5__H`, `laplacian5_coercivity_pass` | exact | PSD | `K >= 1*I` for tridiag(-1,3,-1) |
| `laplacian5_coercivity_fail` | exact | NOT_PSD | `K >= 2*I` refuted; exact witness |
| `corrupted_moments` | exact | NOT_PSD | falsifier demo + CM-table localisation |
| `zeta_x0_1_trunc__{H,L}_N14` | interval | PSD | 60 zeros, truncated measure; transfers to the full measure by Lemma S |
| `zeta_x0_1_full__H_N4` | interval | INCONCLUSIVE | honest width economics of the refutation-grade enclosure |
| `zeta_x0_1_derivative_crosscheck.json` | -- | all_inside = true | non-rigorous `S_Xi` derivatives land inside the rigorous enclosures (validates zero list *and* tail direction) |

Epistemic tiers: exact artifacts are unconditional; zeta artifacts are
*demonstration tier* -- rigorous modulo the [SOURCES.md](SOURCES.md) gates
(Trudgian's explicit `N(T)` constants; `mpmath.zetazero` accuracy within the
stated inflation radius).

## Layout

```
STATUS.md      current satellite heartbeat, blockers, and next small step
src/hausdorff_certificates/
  rational.py    exact LDL^T decision + witness lifting (pure Fractions)
  intervals.py   interval Cholesky, certified Rayleigh witnesses (mpmath.iv)
  hankel.py      H_N, L_N, L_N^theta, Stieltjes S_N, CM difference table
  moments.py     measures -> moments; exact resolvent traces of K; theta dictionary
  zeta.py        zero enclosures, Trudgian tail bound, Lemma-S truncation,
                 non-rigorous derivative cross-check
  certify.py     certificate format (deterministic JSON) + high-level API
  verify.py      independent re-checker (library + CLI)
scripts/generate_artifacts.py   deterministic reference run
tests/           43 tests incl. tamper detection, byte-determinism, manifest digest
MATH.md          precise statements and proofs of everything used
```

## License

MIT. See [CITATION.cff](CITATION.cff) for citing.
