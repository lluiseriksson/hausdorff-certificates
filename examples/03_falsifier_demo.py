"""What a refutation looks like: corrupt one moment of a valid Hausdorff
sequence and watch both detectors fire (Hankel witness + CM table)."""
from fractions import Fraction as F
from hausdorff_certificates import certify_exact, moments_lebesgue

b = moments_lebesgue(14)
b[5] += F(1, 100)  # corruption
cert = certify_exact("corrupted", b, "hankel_L", 6)
print("verdict:", cert.verdict)
print("witness value:", cert.certificate["value"])
print("CM violations (k, n):", cert.cm_check["violations"])
