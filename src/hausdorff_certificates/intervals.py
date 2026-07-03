"""Interval-arithmetic PSD certification (mpmath.iv backend).

Inputs are *enclosures*: symmetric matrices whose entries are intervals
``[lo, hi]`` guaranteed to contain the true (unknown) real entries.

* ``interval_cholesky_pd`` -- if interval Cholesky runs to completion with
  every pivot's lower bound > 0, then **every** symmetric matrix inside the
  enclosure is positive definite (Rump's classical verification argument).
  In particular the true matrix is PD.
* ``certified_negative_rayleigh`` -- given a rational witness vector v, if
  the interval evaluation of v^T [A] v has upper bound < 0, then **every**
  matrix inside the enclosure satisfies v^T A v < 0; in particular the true
  matrix is not PSD.  This is the finite falsifier.
* ``jacobi_eigen_min`` -- a small pure-Python Jacobi eigensolver used only
  as a floating-point *heuristic* to propose witness vectors; it plays no
  role in the certification itself.

Only dependency: mpmath (pure Python, arbitrary precision).
"""

from __future__ import annotations

from contextlib import contextmanager
from fractions import Fraction
from typing import List, Optional, Sequence, Tuple

import mpmath
from mpmath import iv


@contextmanager
def workdps_iv(dps: int):
    """Set BOTH the mp and iv context precisions (mpmath's ``workdps`` only
    touches ``mp``; the interval context has its own independent precision,
    and forgetting it silently degrades every interval to 53 bits)."""
    old_mp, old_iv = mpmath.mp.dps, iv.dps
    mpmath.mp.dps = dps
    iv.dps = dps
    try:
        yield
    finally:
        mpmath.mp.dps = old_mp
        iv.dps = old_iv


class Interval:
    """A thin, picklable wrapper: closed interval with mpf endpoints."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        if hi is None:
            hi = lo
        self.lo = mpmath.mpf(lo)
        self.hi = mpmath.mpf(hi)
        if self.lo > self.hi:
            raise ValueError("empty interval")

    @classmethod
    def from_fraction(cls, fr: Fraction, dps: int) -> "Interval":
        with mpmath.workdps(dps):
            lo = mpmath.mpf(fr.numerator) / mpmath.mpf(fr.denominator)
        # widen by one ulp on each side to stay safe
        eps = mpmath.mpf(10) ** (-(dps - 2))
        rad = abs(lo) * eps + mpmath.mpf(10) ** (-(dps + 10))
        return cls(lo - rad, lo + rad)

    def to_iv(self):
        return iv.mpf([self.lo, self.hi])

    def mid(self):
        return (self.lo + self.hi) / 2

    def __repr__(self):
        return f"[{mpmath.nstr(self.lo, 20)}, {mpmath.nstr(self.hi, 20)}]"


IntervalMatrix = List[List[Interval]]

#: precision used for elementary Interval endpoint operations (outward safe)
INTERVAL_OPS_DPS = 300


def interval_sub(x: "Interval", y: "Interval") -> "Interval":
    """x - y with outward rounding at high precision."""
    with workdps_iv(INTERVAL_OPS_DPS):
        r = iv.mpf([x.lo, x.hi]) - iv.mpf([y.lo, y.hi])
        return Interval(mpmath.mpf(r.a), mpmath.mpf(r.b))


def interval_scale_fraction(c: Fraction, x: "Interval") -> "Interval":
    """c * x for a rational scalar c, with outward rounding."""
    c = Fraction(c)
    with workdps_iv(INTERVAL_OPS_DPS):
        civ = iv.mpf(c.numerator) / iv.mpf(c.denominator)
        r = civ * iv.mpf([x.lo, x.hi])
        return Interval(mpmath.mpf(r.a), mpmath.mpf(r.b))


def _to_iv_matrix(a: IntervalMatrix):
    return [[x.to_iv() for x in row] for row in a]


def interval_cholesky_pd(a: IntervalMatrix, dps: int = 60) -> Tuple[bool, Optional[List]]:
    """Attempt an interval Cholesky factorisation of the enclosure ``a``.

    Returns ``(True, pivot_lower_bounds)`` when the factorisation completes
    with all pivots strictly positive (their lower bounds), which certifies
    that every symmetric matrix in the enclosure is positive definite.
    Returns ``(False, None)`` otherwise -- which is *inconclusive* for PSD
    (the enclosure may simply be too wide), never a refutation.
    """
    n = len(a)
    with workdps_iv(dps):
        A = _to_iv_matrix(a)
        L = [[iv.mpf(0) for _ in range(n)] for _ in range(n)]
        pivots = []
        for j in range(n):
            t = A[j][j]
            for k in range(j):
                t = t - L[j][k] * L[j][k]
            if not (t.a > 0):
                return False, None
            pivots.append(mpmath.mpf(t.a))
            L[j][j] = iv.sqrt(t)
            for i in range(j + 1, n):
                s = A[i][j]
                for k in range(j):
                    s = s - L[i][k] * L[j][k]
                L[i][j] = s / L[j][j]
        return True, pivots


def certified_negative_rayleigh(
    a: IntervalMatrix, v: Sequence[Fraction], dps: int = 60
) -> Tuple[bool, Optional[Tuple[str, str]]]:
    """Evaluate v^T [A] v in interval arithmetic for a rational vector v.

    Returns ``(True, (lo, hi))`` with hi < 0 when the quadratic form is
    certifiably negative for every matrix in the enclosure.
    """
    n = len(a)
    with workdps_iv(dps):
        A = _to_iv_matrix(a)
        viv = []
        for x in v:
            fx = Fraction(x)
            viv.append(iv.mpf(fx.numerator) / iv.mpf(fx.denominator))
        acc = iv.mpf(0)
        for i in range(n):
            for j in range(n):
                acc = acc + viv[i] * A[i][j] * viv[j]
        lo, hi = mpmath.mpf(acc.a), mpmath.mpf(acc.b)
        if hi < 0:
            return True, (mpmath.nstr(lo, 25), mpmath.nstr(hi, 25))
        return False, (mpmath.nstr(lo, 25), mpmath.nstr(hi, 25))


# ----------------------------------------------------------------------
# floating-point helper (heuristic only, not part of any certificate)
# ----------------------------------------------------------------------

def jacobi_eigen_min(a_mid: List[List[float]], sweeps: int = 60) -> Tuple[float, List[float]]:
    """Smallest eigenvalue and eigenvector of a symmetric float matrix
    via cyclic Jacobi rotations (pure Python; small n only)."""
    n = len(a_mid)
    a = [row[:] for row in a_mid]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = 0.0
        for p in range(n - 1):
            for q in range(p + 1, n):
                off += a[p][q] * a[p][q]
        if off < 1e-32:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-300:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
                t = (1.0 if theta >= 0 else -1.0) / (abs(theta) + (theta * theta + 1.0) ** 0.5)
                c = 1.0 / (t * t + 1.0) ** 0.5
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k][p], V[k][q]
                    V[k][p] = c * vkp - s * vkq
                    V[k][q] = s * vkp + c * vkq
    idx = min(range(n), key=lambda i: a[i][i])
    vec = [V[k][idx] for k in range(n)]
    return a[idx][idx], vec


def propose_rational_witness(a: IntervalMatrix, max_den: int = 10**6) -> List[Fraction]:
    """Heuristic witness: rationalised eigenvector of the midpoint matrix
    for its smallest eigenvalue.  Certification happens elsewhere."""
    mid = [[float(x.mid()) for x in row] for row in a]
    _, vec = jacobi_eigen_min(mid)
    scale = max(abs(x) for x in vec) or 1.0
    return [Fraction(x / scale).limit_denominator(max_den) for x in vec]
