"""Yang-Mills-side dictionary: finite coercivity tests for a precision K.

K = tridiag(-1, 3, -1) has lambda_min = 3 - sqrt(3) ~ 1.268.  The claim
K >= 1*I passes; the claim K >= 2*I is refuted with an exact witness."""
from fractions import Fraction as F
from hausdorff_certificates import (
    coercivity_bundle_exact, resolvent_trace_moments, theta_for_coercivity,
)

n = 5
K = [[F(0)] * n for _ in range(n)]
for i in range(n):
    K[i][i] = F(3)
    if i + 1 < n:
        K[i][i + 1] = K[i + 1][i] = F(-1)

_, b = resolvent_trace_moments(K, F(1), 14)
for c_claim in (F(1), F(2)):
    theta = theta_for_coercivity(c_claim, F(1))
    certs = coercivity_bundle_exact(f"K_c{c_claim}", b, 6, theta)
    for cert in certs:
        print(f"claim K >= {c_claim} I | {cert.matrix_kind:>15} -> {cert.verdict}")
        if cert.verdict == "NOT_PSD_CERTIFIED":
            print("   exact witness v =", cert.certificate["v"])
            print("   v^T A v =", cert.certificate["value"], "< 0")
