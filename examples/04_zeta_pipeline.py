"""Bridge Card 8 pipeline (demonstration tier):
zeros -> b_n(x0) -> H_N, L_N -> interval certification.

Consistency direction uses the truncated measure (Lemma S transfers PSD to
the full measure); the refutation-grade full enclosure includes the
Trudgian tail bound.  See MATH.md and SOURCES.md for the epistemic status.
"""
import os
from hausdorff_certificates import certify_interval
from hausdorff_certificates.zeta import (
    load_or_compute_zeros, zeta_moment_enclosures, s_xi_moments_via_derivatives,
)

here = os.path.dirname(__file__)
gammas = load_or_compute_zeros(os.path.join(here, "..", "data", "zeta_zeros_60_dps60.json"), 60, dps=60)

encl = zeta_moment_enclosures(gammas, x0=1, n_max=29, dps=220, inflate_exp=52, include_tail=False)
for kind, N in (("hankel_H", 12), ("hankel_L", 12)):
    cert = certify_interval(f"zeta_{kind}_{N}", encl, kind, N, dps=220, digits=70)
    print(kind, f"N={N} (truncated measure, Lemma S) ->", cert.verdict)

full = zeta_moment_enclosures(gammas, x0=1, n_max=6, dps=220, inflate_exp=52, include_tail=True)
deriv = s_xi_moments_via_derivatives(1, 4, dps=40)  # slow-ish, non-rigorous
for n, d in enumerate(deriv):
    e = full[n]
    print(f"b_{n}: derivative-based {float(d):.8g}  in  [{float(e.lo):.8g}, {float(e.hi):.8g}] ->",
          bool(e.lo <= d <= e.hi))
