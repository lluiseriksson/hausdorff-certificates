# Changelog

## Unreleased

Correct the infinite Hausdorff matrix criterion from `H_N, L_N >= 0` to
`S_N, L_N >= 0`, document the corresponding support criterion, and add the
exact counterexample `b_n = (-1)^n` as a regression test. Existing finite
certificate APIs and artifact formats are unchanged.

## 0.1.0 (2026-07-03)

Initial release: exact-rational and interval PSD certification for
Hausdorff moment sequences (H_N, L_N, L_N^theta, CM table); exact resolvent
-trace moments for coercivity screening; zeta pipeline with Trudgian tail
bound, Lemma-S truncation transfer and non-rigorous derivative cross-check;
independent verifier; deterministic reference artifacts; CI determinism gate.
