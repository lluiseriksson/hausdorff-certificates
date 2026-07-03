#!/usr/bin/env python3
"""Regenerate the reference artifacts, deterministically.

Usage:  python scripts/generate_artifacts.py [--out artifacts] [--skip-zeta]

Every certificate is written with sorted keys, fixed formatting and no
timestamps; `manifest.json` records the sha256 of each file.  CI regenerates
the whole set and diffs it against the committed one (determinism gate),
then re-validates every certificate with the independent verifier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import hausdorff_certificates as hc
from hausdorff_certificates import certify as C
from hausdorff_certificates import moments as M
from hausdorff_certificates import zeta as Z

ZEROS_FILE = os.path.join(os.path.dirname(__file__), "..", "data", "zeta_zeros_60_dps60.json")


def save(cert: C.Certificate, outdir: str, manifest: dict) -> None:
    path = os.path.join(outdir, cert.name + ".cert.json")
    digest = cert.save(path)
    manifest["files"][os.path.basename(path)] = digest
    print(f"  {cert.name:38s} {cert.verdict:18s} sha256={digest[:12]}...")


def save_json(name: str, obj: dict, outdir: str, manifest: dict) -> None:
    path = os.path.join(outdir, name)
    data = json.dumps(obj, sort_keys=True, indent=1) + "\n"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(data)
    manifest["files"][name] = hashlib.sha256(data.encode("utf-8")).hexdigest()
    print(f"  {name:38s} {'(data)':18s} sha256={manifest['files'][name][:12]}...")


def gen_exact(outdir: str, manifest: dict) -> None:
    print("[1/5] Lebesgue on [0,1] (Hilbert matrix), exact backend")
    b = M.moments_lebesgue(18)
    prov = {"source": "Lebesgue measure on [0,1]: b_n = 1/(n+1); H_N is the Hilbert matrix"}
    for cert in C.hausdorff_pair_exact("hilbert_lebesgue", b, 8, provenance=prov):
        save(cert, outdir, manifest)

    print("[2/5] Two-atom measure: support test (true and false claims), exact backend")
    atoms = [(F(1, 2), F(1, 3)), (F(1, 2), F(2, 3))]
    b = M.moments_from_atoms(atoms, 14)
    prov = {"source": "mu = (1/2) delta_{1/3} + (1/2) delta_{2/3}; sup(supp) = 2/3"}
    for cert in C.hausdorff_pair_exact("two_atoms", b, 6, provenance=prov):
        save(cert, outdir, manifest)
    save(
        C.certify_exact(
            "two_atoms_support_true", b, "hankel_L_theta", 6, theta=F(3, 4),
            provenance={**prov, "claim": "supp(mu) subset [0, 3/4] (TRUE)"},
        ),
        outdir, manifest,
    )
    save(
        C.certify_exact(
            "two_atoms_support_false", b, "hankel_L_theta", 6, theta=F(1, 2),
            provenance={**prov, "claim": "supp(mu) subset [0, 1/2] (FALSE: atom at 2/3)",
                        "expected": "NOT_PSD_CERTIFIED with an exact rational witness"},
        ),
        outdir, manifest,
    )

    print("[3/5] Coercivity certificates for K = tridiag(-1, 3, -1) (5x5), exact backend")
    n = 5
    K = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        K[i][i] = F(3)
        if i + 1 < n:
            K[i][i + 1] = K[i + 1][i] = F(-1)
    x0 = F(1)
    _, b = M.resolvent_trace_moments(K, x0, 14)
    lam_min = "3 - 2 cos(pi/6) = 3 - sqrt(3) ~ 1.2679"
    base_prov = {
        "source": "K = tridiag(-1,3,-1), 5x5; b_n = x0^n Tr (K+x0 I)^{-(n+1)}, x0 = 1",
        "lambda_min": lam_min,
        "dictionary": "coercivity K >= c I  <=>  supp subset (0, theta], theta = x0/(c+x0)",
    }
    save(C.certify_exact("laplacian5__H", b, "hankel_H", 6, provenance=base_prov), outdir, manifest)
    save(
        C.certify_exact(
            "laplacian5_coercivity_pass", b, "hankel_L_theta", 6,
            theta=M.theta_for_coercivity(F(1), x0),
            provenance={**base_prov, "claim": "K >= 1*I (TRUE, lambda_min ~ 1.2679)"},
        ),
        outdir, manifest,
    )
    save(
        C.certify_exact(
            "laplacian5_coercivity_fail", b, "hankel_L_theta", 6,
            theta=M.theta_for_coercivity(F(2), x0),
            provenance={**base_prov, "claim": "K >= 2*I (FALSE)",
                        "expected": "NOT_PSD_CERTIFIED: exact witness refutes the claimed constant"},
        ),
        outdir, manifest,
    )

    print("[4/5] Corrupted moment sequence: falsifier demo, exact backend")
    b = M.moments_lebesgue(14)
    bad = list(b)
    bad[5] += F(1, 100)
    save(
        C.certify_exact(
            "corrupted_moments", bad, "hankel_L", 6,
            provenance={"source": "Lebesgue moments with b_5 corrupted by +1/100",
                        "expected": "NOT_PSD_CERTIFIED; cm_check.violations locate (k, n) failures"},
        ),
        outdir, manifest,
    )


def gen_zeta(outdir: str, manifest: dict) -> None:
    print("[5/5] Riemann pipeline (demonstration tier): zeros -> b_n -> H_N, L_N")
    gammas = Z.load_or_compute_zeros(ZEROS_FILE, 60, dps=60)
    x0, n_max = 1, 33
    s_tail = Z.support_bound_of_tail(gammas, x0, dps=60)
    prov_common = {
        "pipeline": "first 60 zeta-zero ordinates (mpmath zetazero, 60 dps, "
                    "inflated +-1e-52) -> b_n(x0) = x0^n sum (gamma^2+x0)^{-(n+1)}, x0 = 1",
        "epistemic_status": "demonstration tier: rigorous modulo (i) Trudgian (2014) "
                            "N(T) constants (SOURCES.md gate V-A) and (ii) mpmath zetazero "
                            "accuracy within the inflation radius (gate V-B)",
        "knowledge_tree_node": "bridge.hausdorff_moment_certificates",
    }
    # (a) consistency certificates on the TRUNCATED measure (monotone splitting)
    encl_t = Z.zeta_moment_enclosures(gammas, x0, n_max, dps=220, inflate_exp=52, include_tail=False)
    prov_t = {
        **prov_common,
        "mode": "truncated measure mu_trunc (60 zeros, no tail)",
        "transfer": "Lemma S (MATH.md): mu = mu_trunc + mu_tail with mu_tail >= 0 "
                    f"supported in [0, s], s <= {s_tail}; hence H_N(mu) >= H_N(mu_trunc) "
                    "and L_N(mu) >= L_N(mu_trunc); PSD below transfers to the full measure.",
        "meaning": "consistency with RH at (x0, N); a PSD outcome can never prove RH",
    }
    save(C.certify_interval("zeta_x0_1_trunc__H_N14", encl_t, "hankel_H", 14, dps=220,
                            provenance=prov_t, digits=70), outdir, manifest)
    save(C.certify_interval("zeta_x0_1_trunc__L_N14", encl_t, "hankel_L", 14, dps=220,
                            provenance=prov_t, digits=70), outdir, manifest)
    # (b) full-measure enclosures (tail included): the refutation-grade object
    encl_f = Z.zeta_moment_enclosures(gammas, x0, 13, dps=220, inflate_exp=52, include_tail=True)
    prov_f = {
        **prov_common,
        "mode": "full measure: truncated interval sum + Trudgian tail upper bound",
        "meaning": "any NOT_PSD verdict on these enclosures would refute RH "
                   "(modulo the SOURCES.md gates); PSD verdicts remain consistency only",
        "expected": "INCONCLUSIVE at this depth: the tail bound gives b_0 a width "
                    "~4.3e-3 (only 60 zeros, T ~ 163), which drowns the ~1e-12 "
                    "Cholesky pivots from N = 2 on. The refutation-grade enclosure "
                    "narrows like log(T)/T as more zeros are added; the consistency "
                    "direction does not need it (Lemma S). This artifact documents "
                    "the honest width economics of the falsifier.",
    }
    save(C.certify_interval("zeta_x0_1_full__H_N4", encl_f, "hankel_H", 4, dps=220,
                            provenance=prov_f, digits=70), outdir, manifest)
    # (c) NON-RIGOROUS cross-check: derivative-based b_n inside the enclosures
    deriv = Z.s_xi_moments_via_derivatives(x0, 6, dps=40, radius="0.5")
    import mpmath
    rows = []
    all_inside = True
    for n, d in enumerate(deriv):
        e = encl_f[n]
        inside = bool(e.lo <= d <= e.hi)
        all_inside &= inside
        rows.append({
            "n": n,
            "b_n_from_derivatives": mpmath.nstr(d, 20),
            "enclosure": [mpmath.nstr(e.lo, 20), mpmath.nstr(e.hi, 20)],
            "inside": inside,
        })
    save_json("zeta_x0_1_derivative_crosscheck.json", {
        "description": "b_n(1) from Cauchy-integral derivatives of S_Xi via "
                       "xi'/xi (floating quadrature, NON-RIGOROUS) compared against "
                       "the rigorous zero-sum + tail enclosures; 'inside' validates "
                       "both the zero list and the tail bound direction",
        "all_inside": all_inside,
        "rows": rows,
        "epistemic_status": "demonstration tier; no certificate semantics",
    }, outdir, manifest)
    if not all_inside:
        raise SystemExit("cross-check FAILED: derivative moments left the enclosures")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "artifacts"))
    ap.add_argument("--skip-zeta", action="store_true")
    args = ap.parse_args()
    outdir = os.path.abspath(args.out)
    os.makedirs(outdir, exist_ok=True)
    manifest = {
        "format": "hausdorff-certificates-manifest/1",
        "package_version": hc.__version__,
        "parameters": {
            "zeta": {"zeros": 60, "zeros_dps": 60, "x0": "1", "inflate_exp": 52,
                     "work_dps": 220, "serialize_digits": 70},
        },
        "files": {},
    }
    gen_exact(outdir, manifest)
    if not args.skip_zeta:
        gen_zeta(outdir, manifest)
    path = os.path.join(outdir, "manifest.json")
    data = json.dumps(manifest, sort_keys=True, indent=1) + "\n"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(data)
    print(f"manifest: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
