# Status heartbeat

Last audited: 2026-07-04 13:48 UTC

Default branch `main` is at
`ce2ca84d47792a0e2edb3f25128aa171330ca44c`. The latest observed GitHub
Actions run on `main` was workflow `ci`, run `28707084762`, and completed
successfully on 2026-07-04.

## Mother-facing digest

This repository is a finite-certificate satellite. It produces deterministic
JSON certificates for finite Hausdorff-moment matrix checks, then revalidates
those certificates without trusting the generator.

Consumable API and theorem names:

- `certify_exact`, `certify_interval`, `verify_obj`
- `hankel_H`, `hankel_L`, `hankel_L_theta`
- `hausdorff_pair_exact`, `coercivity_bundle_exact`
- `theta_for_coercivity`, `lambda_min_lower_bound_from_theta`
- `zeta_moment_enclosures`, `tail_upper_bound`, `lemma_s_tail_bound`

Consumable command:

- `python -m hausdorff_certificates.manifest artifacts/manifest.json`

Consumable files:

- `MATH.md` contains the finite Hausdorff, support-bound, coercivity, zeta-tail,
  and certificate-verification statements used by this repo.
- `SOURCES.md` lists the explicit external gates for the zeta demonstration
  tier.
- `REPRODUCIBILITY.md` defines the deterministic artifact and verifier
  contract.
- `artifacts/manifest.json` names the committed reference certificate batch.

Current reference certificate roles:

- exact unconditional examples: `hilbert_lebesgue__{H,L}`,
  `two_atoms__{H,L}`, `two_atoms_support_true`,
  `two_atoms_support_false`, `laplacian5__H`,
  `laplacian5_coercivity_pass`, `laplacian5_coercivity_fail`,
  `corrupted_moments`
- zeta demonstration-tier examples: `zeta_x0_1_trunc__{H,L}_N14`,
  `zeta_x0_1_full__H_N4`, `zeta_x0_1_derivative_crosscheck.json`

## Exact current blockers

- Gate V-A in `SOURCES.md`: independently check the Trudgian constants and
  validity range before treating zeta-tail certificates as final.
- Gate V-B in `SOURCES.md`: independently cross-check the committed
  `mpmath.zetazero` ordinates before treating zeta artifacts as final.

These blockers do not affect the exact-rational artifacts.

## Next small step

Add a CLI smoke test for `python -m hausdorff_certificates.verify` on one
positive exact certificate and one committed negative witness certificate.
