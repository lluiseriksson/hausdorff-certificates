# Reproducibility

## Determinism contract

- Certificates contain no timestamps, no environment data, no floats
  formatted by locale; JSON is written with sorted keys and `indent=1`.
- All exact-backend content is pure `fractions.Fraction`: identical bytes on
  any platform/Python >= 3.10.
- All interval-backend content is mpmath (pure Python arbitrary precision):
  identical bytes given the same parameters, which are pinned in
  `scripts/generate_artifacts.py` and recorded in `artifacts/manifest.json`.
- The zeta zero cache `data/zeta_zeros_60_dps60.json` is committed; `make
  zeros` regenerates it deterministically (`mpmath.zetazero`, 60 dps).

## Reproducing the reference run

```bash
pip install -e ".[test]"
make artifacts        # writes artifacts/ (including manifest.json)
make verify           # independent re-validation of every certificate
git diff --stat artifacts/   # must be empty against the committed set
```

CI (`.github/workflows/ci.yml`) performs exactly this on Python 3.10-3.12:
tests, regeneration into a scratch directory, byte-diff against the
committed artifacts, independent verification, artifact upload.

## What "verified" means here

`python -m hausdorff_certificates.verify FILE` rebuilds the matrix from the
stored moments and re-checks the stored evidence only:

- `ldlt`: exact reconstruction `P^T A P == L D L^T`, `D >= 0`;
- `negative_witness`: exact `v^T A v < 0` (and equality with the stored value);
- `interval_cholesky`: re-run of the factorisation at the stored precision;
- `interval_negative_witness`: interval Rayleigh upper bound < 0.

Trust in the generator is not required; tampering with any moment or factor
makes verification fail (covered by `test_certificate_tamper_detected`).
