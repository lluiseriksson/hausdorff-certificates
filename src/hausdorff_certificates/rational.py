"""Exact positive-semidefiniteness decision over the rationals.

For a symmetric matrix with ``fractions.Fraction`` entries this module
decides PSD *exactly* and returns an independently checkable certificate:

* ``PSD``      -- a permutation P, unit lower-triangular L and diagonal
                  D >= 0 with  P^T A P = L D L^T  (exact identity).
* ``NOT_PSD``  -- a rational witness vector v with  v^T A v < 0
                  (exact strict inequality).

The algorithm is LDL^T with symmetric (diagonal) pivoting.  For a PSD
matrix a zero diagonal entry of the running Schur complement forces its
whole row/column to vanish; if it does not, a 2x2 indefinite block gives
an explicit witness, which is lifted back through the elimination via
back-substitution and finally re-checked against the *original* matrix.

Everything is pure Python + ``fractions``; no external dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import List, Optional, Sequence

Matrix = List[List[Fraction]]


def as_fraction_matrix(rows: Sequence[Sequence]) -> Matrix:
    return [[Fraction(x) for x in row] for row in rows]


def check_symmetric(a: Matrix) -> None:
    n = len(a)
    for row in a:
        if len(row) != n:
            raise ValueError("matrix is not square")
    for i in range(n):
        for j in range(i + 1, n):
            if a[i][j] != a[j][i]:
                raise ValueError(f"matrix is not symmetric at ({i},{j})")


def quadratic_form(a: Matrix, v: Sequence[Fraction]) -> Fraction:
    n = len(a)
    return sum(Fraction(v[i]) * a[i][j] * Fraction(v[j]) for i in range(n) for j in range(n))


@dataclass
class ExactPSDResult:
    """Outcome of the exact PSD decision."""

    is_psd: bool
    n: int
    #: permutation as a list: column i of P is e_{perm[i]} (i.e. pivot order)
    perm: Optional[List[int]] = None
    #: unit lower-triangular factor (n x n, Fractions), only if is_psd
    L: Optional[Matrix] = None
    #: pivots, only if is_psd (all >= 0)
    D: Optional[List[Fraction]] = None
    #: rational witness with witness^T A witness < 0, only if not is_psd
    witness: Optional[List[Fraction]] = None
    #: exact value of the witness quadratic form (negative), if not is_psd
    witness_value: Optional[Fraction] = None
    rank: Optional[int] = None

    def verify(self, a_original: Matrix) -> bool:
        """Re-check the certificate against the original matrix, exactly."""
        n = self.n
        if self.is_psd:
            assert self.perm is not None and self.L is not None and self.D is not None
            if any(d < 0 for d in self.D):
                return False
            # reconstruct P^T A P and compare with L D L^T entrywise
            p = self.perm
            for i in range(n):
                for j in range(n):
                    ldlt = sum(self.L[i][k] * self.D[k] * self.L[j][k] for k in range(min(i, j) + 1))
                    if ldlt != a_original[p[i]][p[j]]:
                        return False
            return True
        assert self.witness is not None
        return quadratic_form(a_original, self.witness) < 0


def _lift_witness(L: Matrix, perm: List[int], step: int, w: List[Fraction], n: int) -> List[Fraction]:
    """Lift a Schur-complement witness ``w`` (living on permuted coordinates
    ``step..n-1``) to a witness for the original matrix.

    If  P^T A P = [[A11, A12], [A21, A22]]  and  S = A22 - A21 A11^{-1} A12
    with  w^T S w < 0, then  v_perm = (-A11^{-1} A12 w, w)  satisfies
    v_perm^T (P^T A P) v_perm = w^T S w.  Using the stored factors,
    -A11^{-1} A12 w = -L11^{-T} (L21^T w),  a single back-substitution.
    """
    # t = L21^T w  (length = step)
    t = [Fraction(0)] * step
    for k in range(step):
        t[k] = sum(L[i][k] * w[i - step] for i in range(step, n))
    # solve L11^T x = -t  (back substitution; L11 unit lower => L11^T unit upper)
    x = [Fraction(0)] * step
    for k in range(step - 1, -1, -1):
        s = -t[k] - sum(L[j][k] * x[j] for j in range(k + 1, step))
        x[k] = s
    v_perm = x + list(w)
    # un-permute: v_perm lives on coordinates perm[0..n-1]
    v = [Fraction(0)] * n
    for i in range(n):
        v[perm[i]] = v_perm[i]
    return v


def exact_psd(a_in: Sequence[Sequence]) -> ExactPSDResult:
    """Decide PSD for a symmetric rational matrix, with certificate."""
    a_original = as_fraction_matrix(a_in)
    check_symmetric(a_original)
    n = len(a_original)
    # working copy (will hold successive Schur complements on rows/cols step..n-1)
    s = [row[:] for row in a_original]
    perm = list(range(n))
    L: Matrix = [[Fraction(1) if i == j else Fraction(0) for j in range(n)] for i in range(n)]
    D: List[Fraction] = [Fraction(0)] * n
    rank = 0

    def swap(i: int, j: int) -> None:
        if i == j:
            return
        perm[i], perm[j] = perm[j], perm[i]
        s[i], s[j] = s[j], s[i]
        for r in range(n):
            s[r][i], s[r][j] = s[r][j], s[r][i]
        # swap already-built rows of L (columns < current step only)
        L[i], L[j] = L[j], L[i]
        L[i][i], L[i][j] = Fraction(1), Fraction(0)
        L[j][j], L[j][i] = Fraction(1), Fraction(0)

    step = 0
    while step < n:
        # 1) any negative diagonal entry in the Schur complement refutes PSD
        neg = next((j for j in range(step, n) if s[j][j] < 0), None)
        if neg is not None:
            w = [Fraction(0)] * (n - step)
            w[neg - step] = Fraction(1)
            v = _lift_witness(L, perm, step, w, n)
            val = quadratic_form(a_original, v)
            assert val < 0, "internal error: lifted witness not negative"
            return ExactPSDResult(False, n, witness=v, witness_value=val)
        # 2) choose the largest positive diagonal pivot, if any
        pos = [j for j in range(step, n) if s[j][j] > 0]
        if pos:
            piv = max(pos, key=lambda j: s[j][j])
            swap(step, piv)
            d = s[step][step]
            D[step] = d
            rank += 1
            for i in range(step + 1, n):
                L[i][step] = s[i][step] / d
            for i in range(step + 1, n):
                li = L[i][step]
                if li == 0:
                    continue
                for j in range(step + 1, n):
                    s[i][j] -= li * d * L[j][step]
            step += 1
            continue
        # 3) all remaining diagonal entries are zero: PSD forces the block to vanish
        off = None
        for i in range(step, n):
            for j in range(i + 1, n):
                if s[i][j] != 0:
                    off = (i, j)
                    break
            if off:
                break
        if off is None:
            break  # remaining Schur complement is zero -> PSD (rank deficient)
        i, j = off
        w = [Fraction(0)] * (n - step)
        w[i - step] = Fraction(1)
        w[j - step] = Fraction(-1) if s[i][j] > 0 else Fraction(1)
        v = _lift_witness(L, perm, step, w, n)
        val = quadratic_form(a_original, v)
        assert val < 0, "internal error: lifted 2x2 witness not negative"
        return ExactPSDResult(False, n, witness=v, witness_value=val)

    res = ExactPSDResult(True, n, perm=perm, L=L, D=D, rank=rank)
    assert res.verify(a_original), "internal error: LDL^T certificate failed self-check"
    return res


def solve_linear_system(a_in: Sequence[Sequence], b_cols: Sequence[Sequence[Fraction]]) -> Matrix:
    """Solve A X = B exactly over Q (Gaussian elimination with partial pivoting).

    ``b_cols`` is given column-wise; returns X column-wise.
    Raises ``ZeroDivisionError`` if A is singular.
    """
    a = as_fraction_matrix(a_in)
    n = len(a)
    cols = [list(map(Fraction, c)) for c in b_cols]
    m = len(cols)
    # augmented elimination
    for k in range(n):
        piv = next((r for r in range(k, n) if a[r][k] != 0), None)
        if piv is None:
            raise ZeroDivisionError("singular matrix")
        if piv != k:
            a[k], a[piv] = a[piv], a[k]
            for c in cols:
                c[k], c[piv] = c[piv], c[k]
        pv = a[k][k]
        for r in range(k + 1, n):
            f = a[r][k] / pv
            if f == 0:
                continue
            for c2 in range(k, n):
                a[r][c2] -= f * a[k][c2]
            for c in cols:
                c[r] -= f * c[k]
    xs: Matrix = [[Fraction(0)] * n for _ in range(m)]
    for ci in range(m):
        x = xs[ci]
        b = cols[ci]
        for k in range(n - 1, -1, -1):
            sacc = b[k] - sum(a[k][j] * x[j] for j in range(k + 1, n))
            x[k] = sacc / a[k][k]
    return xs
