"""Moment-sequence constructors.

Three sources of Hausdorff moment sequences ``b_n = \\int_0^1 v^n d\\mu``:

1. explicit atomic/absolutely-continuous measures (exact, for tests),
2. a coercive symmetric rational matrix K via exact resolvent traces
   (the Yang--Mills-side adapter),
3. the zeta-zero measure (interval enclosures; see :mod:`.zeta`).

Yang--Mills-side dictionary (MATH.md, Bridge Card 8):

    M_n(K; x0) = Tr (K + x0 I)^{-(n+1)} = sum_j (lambda_j + x0)^{-(n+1)}
    b_n        = x0^n * M_n = sum_j u_j v_j^n,
                 u_j = 1/(lambda_j + x0),  v_j = x0/(lambda_j + x0) in (0, 1).

Coercivity  lambda_min >= c  is equivalent to  supp subset (0, theta] with
theta = x0/(c + x0); the finite certificates are  H_N >= 0  and
L_N^theta >= 0  (see :func:`hausdorff_certificates.hankel.hankel_L_theta`).
"""

from __future__ import annotations

from fractions import Fraction
from typing import List, Sequence, Tuple

from .rational import Matrix, as_fraction_matrix, check_symmetric, solve_linear_system


def moments_from_atoms(atoms: Sequence[Tuple[Fraction, Fraction]], n_max: int) -> List[Fraction]:
    """b_n for a finite atomic measure  mu = sum_j u_j delta_{v_j}.

    ``atoms`` is a list of (weight u_j, position v_j) with u_j >= 0 and
    v_j in [0, 1] for a genuine Hausdorff measure (not enforced -- the whole
    point of the certificates is to detect violations).
    """
    b = []
    for n in range(n_max + 1):
        b.append(sum(Fraction(u) * (Fraction(v) ** n) for (u, v) in atoms))
    return b


def moments_lebesgue(n_max: int) -> List[Fraction]:
    """b_n = 1/(n+1)  (Lebesgue measure on [0,1]; H_N is the Hilbert matrix)."""
    return [Fraction(1, n + 1) for n in range(n_max + 1)]


def resolvent_trace_moments(
    K_in: Sequence[Sequence], x0, n_max: int
) -> Tuple[List[Fraction], List[Fraction]]:
    """Exact  M_n(K; x0) = Tr (K + x0 I)^{-(n+1)}  and  b_n = x0^n M_n  over Q.

    ``K_in`` must be symmetric with rational entries and  K + x0 I  must be
    invertible (it is whenever K is coercive and x0 > 0).  Returns
    ``(M, b)`` with ``len == n_max + 1``.
    """
    K = as_fraction_matrix(K_in)
    check_symmetric(K)
    x0 = Fraction(x0)
    n = len(K)
    A = [[K[i][j] + (x0 if i == j else 0) for j in range(n)] for i in range(n)]
    # R = A^{-1} exactly: solve A X = I column by column
    eye = [[Fraction(1) if i == j else Fraction(0) for i in range(n)] for j in range(n)]
    cols = solve_linear_system(A, eye)  # columns of A^{-1}
    R: Matrix = [[cols[j][i] for j in range(n)] for i in range(n)]
    # powers R^{p}, p = 1..n_max+1, take traces
    M: List[Fraction] = []
    P = [row[:] for row in R]
    for p in range(1, n_max + 2):
        M.append(sum(P[i][i] for i in range(n)))
        if p <= n_max:
            P = [
                [sum(P[i][k] * R[k][j] for k in range(n)) for j in range(n)]
                for i in range(n)
            ]
    b = [(x0 ** nn) * M[nn] for nn in range(n_max + 1)]
    return M, b


def theta_for_coercivity(c, x0) -> Fraction:
    """theta = x0/(c + x0): claimed coercivity constant -> support bound."""
    c, x0 = Fraction(c), Fraction(x0)
    return x0 / (c + x0)


def lambda_min_lower_bound_from_theta(theta, x0) -> Fraction:
    """Invert: supp subset [0, theta]  <=>  lambda_min >= x0 (1 - theta)/theta."""
    theta, x0 = Fraction(theta), Fraction(x0)
    return x0 * (1 - theta) / theta
