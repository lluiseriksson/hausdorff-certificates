"""Hausdorff-moment certificate matrices.

For a sequence ``b = (b_0, b_1, ..., b_M)`` the classical Hausdorff moment
problem on [0, 1] is solvable (i.e. b_n = \\int_0^1 v^n d\\mu with \\mu >= 0)
iff the sequence is completely monotone. An exact infinite matrix
characterization uses, for every N,

    S_N = ( b_{i+j+1}          )_{0<=i,j<=N}  >= 0        (moments of v d\\mu)
    L_N = ( b_{i+j} - b_{i+j+1})_{0<=i,j<=N}  >= 0        (moments of (1-v) d\\mu).

Then H_N = S_N + L_N is automatically positive semidefinite. Positivity of
H_N and L_N alone is necessary but not sufficient: b_n = (-1)^n makes both
families positive semidefinite although it is represented by delta_{-1}, not
by a positive measure on [0, 1].

For an abstract sequence, the generalized upper-support characterization is

    S_N >= 0   and   L_N^\\theta = ( \\theta b_{i+j} - b_{i+j+1} ) >= 0    for all N.

Indeed theta H_N = S_N + L_N^theta. If a positive representing measure is
already known, L_N^theta alone characterizes support in [0, theta]. The H/L
and H/L-theta convenience bundles in this package remain useful finite
necessary screens, but they are not complete abstract moment decisions.

All builders work for both exact ``Fraction`` sequences and ``Interval``
sequences; the arithmetic is whatever the elements support.
"""

from __future__ import annotations

from fractions import Fraction
from typing import List, Sequence

from .intervals import Interval, interval_scale_fraction, interval_sub


def _sub(x, y):
    if isinstance(x, Interval) or isinstance(y, Interval):
        return interval_sub(x, y)
    return x - y


def _scale(c, x):
    if isinstance(x, Interval):
        return interval_scale_fraction(Fraction(c), x)
    return c * x


def hankel_H(b: Sequence, N: int) -> List[List]:
    """H_N = (b_{i+j}), size (N+1)x(N+1); requires len(b) >= 2N+1."""
    if len(b) < 2 * N + 1:
        raise ValueError(f"need {2*N+1} moments for H_{N}, got {len(b)}")
    return [[b[i + j] for j in range(N + 1)] for i in range(N + 1)]


def hankel_L(b: Sequence, N: int) -> List[List]:
    """L_N = (b_{i+j} - b_{i+j+1}); requires len(b) >= 2N+2."""
    if len(b) < 2 * N + 2:
        raise ValueError(f"need {2*N+2} moments for L_{N}, got {len(b)}")
    return [[_sub(b[i + j], b[i + j + 1]) for j in range(N + 1)] for i in range(N + 1)]


def hankel_L_theta(b: Sequence, N: int, theta) -> List[List]:
    """L_N^theta = (theta * b_{i+j} - b_{i+j+1}), the support-in-[0,theta] test."""
    if len(b) < 2 * N + 2:
        raise ValueError(f"need {2*N+2} moments for L^theta_{N}, got {len(b)}")
    th = theta if isinstance(theta, (Fraction, int)) else Fraction(theta)
    return [
        [_sub(_scale(th, b[i + j]), b[i + j + 1]) for j in range(N + 1)]
        for i in range(N + 1)
    ]


def shifted_hankel_S(b: Sequence, N: int) -> List[List]:
    """S_N = (b_{i+j+1}), the lower-support localizing matrix."""
    if len(b) < 2 * N + 2:
        raise ValueError(f"need {2*N+2} moments for S_{N}, got {len(b)}")
    return [[b[i + j + 1] for j in range(N + 1)] for i in range(N + 1)]


def complete_monotonicity_table(b: Sequence, max_order: int | None = None):
    """Finite differences  (-1)^k Delta^k b_n  for  n + k <= len(b)-1.

    For a Hausdorff sequence every entry equals \\int v^n (1-v)^k d\\mu >= 0.
    Returns a list of rows ``table[k][n]``.  A negative entry (for exact
    input) or an entry whose interval upper bound is negative (for interval
    input) is a localized falsification witness (n, k).
    """
    M = len(b) - 1
    if max_order is None:
        max_order = M
    rows = [list(b)]
    for k in range(1, max_order + 1):
        prev = rows[-1]
        # Delta applied k times, sign-flipped:  (-1)^k Delta^k b_n
        # recursion: c^{(k)}_n = c^{(k-1)}_n - c^{(k-1)}_{n+1}
        rows.append([_sub(prev[n], prev[n + 1]) for n in range(len(prev) - 1)])
    return rows


def cm_violations(table) -> List[tuple]:
    """Indices (k, n) whose entry is certifiably negative."""
    out = []
    for k, row in enumerate(table):
        for n, x in enumerate(row):
            if isinstance(x, Interval):
                if x.hi < 0:
                    out.append((k, n))
            else:
                if x < 0:
                    out.append((k, n))
    return out
