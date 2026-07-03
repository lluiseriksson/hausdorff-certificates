"""Lebesgue on [0,1]: b_n = 1/(n+1), H_N is the Hilbert matrix.  Exact PSD
certificate with independently checkable LDL^T factors."""
from hausdorff_certificates import moments_lebesgue, hausdorff_pair_exact

b = moments_lebesgue(18)
for cert in hausdorff_pair_exact("hilbert", b, 8):
    print(cert.name, "->", cert.verdict)
    print(" first pivots:", cert.certificate["D"][:4])
