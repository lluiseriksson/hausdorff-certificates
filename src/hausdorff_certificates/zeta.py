"""Riemann-side moment enclosures (demonstration tier).

Bridge Card 8 pipeline:   zeros -> b_n -> H_N, L_N -> interval certification.

For x0 > 1/4 define (single-point normalized derivatives of the resolvent
trace S_Xi(x) = sum_{gamma>0} 1/(gamma^2 + x)):

    b_n(x0) = x0^n * sum_{gamma>0} (gamma^2 + x0)^{-(n+1)}
            = int_0^1 v^n dmu_{x0}(v),     v_gamma = x0/(gamma^2 + x0).

This module produces rigorous *enclosures* of b_n from

  (a) interval evaluations over the first Z zeros (mpmath ``zetazero``
      values inflated by an explicit radius), plus
  (b) an explicit upper bound for the tail  sum_{gamma > T} f(gamma),
      f(t) = (t^2 + x0)^{-(n+1)}, via Stieltjes integration by parts against
      the explicit Riemann--von Mangoldt bound of Trudgian (2014):

        |N(t) - (t/2pi) log(t/(2pi e)) - 7/8|
              <= C1 log t + C2 log log t + C3,   t >= e,
        (C1, C2, C3) = (0.112, 0.278, 2.510).

      Writing N_up for the resulting upper bound and using that f decreases
      to 0 and N(T) >= Z (we exhibited Z zeros of height <= T):

        sum_{gamma > T} f(gamma)
            <= f(T) * max(0, N_up(T) - Z) + int_T^inf f(t) N_up'(t) dt,

      and with  f(t) <= t^{-m},  m = 2(n+1),  the integral is bounded by
      elementary closed forms (see ``tail_upper_bound``).

EPISTEMIC STATUS (see SOURCES.md): the enclosures are rigorous *modulo*
(i) the literature constants (C1, C2, C3), listed as a verification-queue
item, and (ii) the accuracy of mpmath's ``zetazero`` within the stated
inflation radius (cross-check against Odlyzko's tables is the gate).
Everything on the exact-rational side of this repository is independent of
this module.
"""

from __future__ import annotations

import json
import os
from typing import List, Tuple

import mpmath
from mpmath import iv

from .intervals import Interval, workdps_iv

TRUDGIAN_C1 = "0.112"
TRUDGIAN_C2 = "0.278"
TRUDGIAN_C3 = "2.510"


def compute_zero_ordinates(count: int, dps: int = 40) -> List[str]:
    """First ``count`` positive ordinates gamma_k of zeta zeros, as decimal
    strings at ``dps`` digits (deterministic given (count, dps))."""
    out = []
    with mpmath.workdps(dps):
        for k in range(1, count + 1):
            rho = mpmath.zetazero(k)
            out.append(mpmath.nstr(rho.imag, dps, strip_zeros=False))
    return out


def load_or_compute_zeros(path: str, count: int, dps: int = 40) -> List[str]:
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        if data.get("count") == count and data.get("dps") == dps:
            return data["gammas"]
    gammas = compute_zero_ordinates(count, dps)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"count": count, "dps": dps, "gammas": gammas}, fh, indent=1)
    return gammas


def _gamma_intervals(gammas: List[str], dps: int, inflate_exp: int):
    """Interval enclosures of the ordinates: value +- 10^{-inflate_exp}."""
    rad = mpmath.mpf(10) ** (-inflate_exp)
    out = []
    for g in gammas:
        x = mpmath.mpf(g)
        out.append(iv.mpf([x - rad, x + rad]))
    return out


def tail_upper_bound(T, x0, n: int, Z: int, dps: int = 60):
    """Certified upper bound (mpf) for  sum_{gamma > T} (gamma^2+x0)^{-(n+1)},
    given that Z zeros with ordinate <= T have been exhibited."""
    with workdps_iv(dps):
        Tiv = iv.mpf(T)
        x0iv = iv.mpf(x0)
        m = 2 * (n + 1)
        two_pi = 2 * iv.pi
        c1 = iv.mpf(TRUDGIAN_C1)
        c2 = iv.mpf(TRUDGIAN_C2)
        c3 = iv.mpf(TRUDGIAN_C3)
        # N_up(T)
        N_up_T = (Tiv / two_pi) * iv.log(Tiv / (two_pi * iv.e)) + iv.mpf("0.875") \
            + c1 * iv.log(Tiv) + c2 * iv.log(iv.log(Tiv)) + c3
        fT = (Tiv * Tiv + x0iv) ** (-(n + 1))
        boundary = fT * N_up_T - fT * iv.mpf(Z)
        if boundary.b < 0:
            boundary = iv.mpf(0)
        # integral pieces with f(t) <= t^{-m}
        T1m = Tiv ** (1 - m)      # T^{1-m}
        Tm = Tiv ** (-m)          # T^{-m}
        main = (iv.mpf(1) / two_pi) * T1m * (
            iv.log(Tiv / two_pi) / iv.mpf(m - 1) + iv.mpf(1) / iv.mpf((m - 1) ** 2)
        )
        piece_c1 = c1 * Tm / iv.mpf(m)
        piece_c2 = (c2 / iv.log(Tiv)) * Tm / iv.mpf(m)
        total = boundary + main + piece_c1 + piece_c2
        return mpmath.mpf(total.b)  # upper endpoint as mpf


def zeta_moment_enclosures(
    gammas: List[str],
    x0,
    n_max: int,
    dps: int = 60,
    inflate_exp: int = 30,
    include_tail: bool = True,
) -> List[Interval]:
    """Enclosures [lo, hi] of  b_n(x0),  n = 0..n_max.

    With ``include_tail=True`` (default): lo = interval lower bound of the
    truncated sum (valid because all terms are positive), hi = interval upper
    bound of the truncated sum + the Trudgian tail bound; the result encloses
    the moments of the FULL zero measure -- required for any refutation.

    With ``include_tail=False``: tight enclosures of the moments of the
    TRUNCATED measure mu_trunc (first Z zeros only).  By the monotone
    splitting lemma (MATH.md, Lemma S) mu = mu_trunc + mu_tail with
    mu_tail >= 0 supported in [0, s], s = x0/(gamma_Z^2 + x0), so

        H_N(mu) >= H_N(mu_trunc),   L_N(mu) >= L_N(mu_trunc),
        L_N^theta(mu) >= L_N^theta(mu_trunc)   whenever theta >= s,

    and PSD certificates for the truncated matrices transfer to the full
    measure.  Only the consistency (PSD) direction may use this mode.
    """
    Z = len(gammas)
    out: List[Interval] = []
    with workdps_iv(dps):
        gs = _gamma_intervals(gammas, dps, inflate_exp)
        x0iv = iv.mpf(x0)
        # T = certified lower bound on the largest exhibited ordinate
        T = mpmath.mpf(gs[-1].a)
        # precompute u_k = 1/(gamma_k^2 + x0) as intervals
        us = [iv.mpf(1) / (g * g + x0iv) for g in gs]
        for n in range(n_max + 1):
            acc = iv.mpf(0)
            for u in us:
                acc = acc + u ** (n + 1)
            scale = iv.mpf(x0) ** n
            full = scale * acc
            lo = mpmath.mpf(full.a)
            if include_tail:
                tail_hi = tail_upper_bound(T, x0, n, Z, dps=dps)
                hi = mpmath.mpf((scale * (acc + iv.mpf([0, tail_hi]))).b)
            else:
                hi = mpmath.mpf(full.b)
            out.append(Interval(lo, hi))
    return out


def support_bound_of_tail(gammas: List[str], x0, dps: int = 60) -> str:
    """Certified upper bound s for the support of the tail measure:
    every omitted zero has ordinate > T = (lower bound of gamma_Z), hence
    v = x0/(gamma^2 + x0) < x0/(T^2 + x0) =: s.  Decimal string (rounded up)."""
    with workdps_iv(dps):
        rad = mpmath.mpf(10) ** (-10)
        T = mpmath.mpf(gammas[-1]) - rad
        s = iv.mpf(x0) / (iv.mpf(T) * iv.mpf(T) + iv.mpf(x0))
        return mpmath.nstr(mpmath.mpf(s.b), 25)


# ----------------------------------------------------------------------
# NON-RIGOROUS cross-check: b_n from derivatives of S_Xi itself
# ----------------------------------------------------------------------

def s_xi(x, dps: int = 40):
    """S_Xi(x) = (1/(2 sqrt x)) * (xi'/xi)(1/2 + sqrt x)  via

        xi'/xi(s) = 1/s + 1/(s-1) - (1/2) log pi
                    + (1/2) psi(s/2) + zeta'(s)/zeta(s).

    Plain floating evaluation (mpmath), NOT interval arithmetic."""
    with mpmath.workdps(dps):
        x = mpmath.mpf(x) if not isinstance(x, mpmath.mpc) else x
        rt = mpmath.sqrt(x)
        sarg = mpmath.mpf("0.5") + rt
        log_der = (
            1 / sarg
            + 1 / (sarg - 1)
            - mpmath.log(mpmath.pi) / 2
            + mpmath.digamma(sarg / 2) / 2
            + mpmath.zeta(sarg, derivative=1) / mpmath.zeta(sarg)
        )
        return log_der / (2 * rt)


def s_xi_moments_via_derivatives(x0, n_max: int, dps: int = 40, radius="0.5"):
    """b_n(x0) = x0^n (-1)^n S_Xi^{(n)}(x0)/n!  computed from Cauchy-integral
    derivatives of S_Xi (mpmath ``diff`` with method='quad').

    DEMONSTRATION TIER: floating-point quadrature, no certified error bound.
    Used only to cross-check that the derivative-based values fall inside
    the rigorous zero-sum enclosures (which also validates the tail bound).
    """
    out = []
    with mpmath.workdps(dps):
        x0m = mpmath.mpf(x0)
        for n in range(n_max + 1):
            d = mpmath.diff(lambda z: s_xi(z, dps=dps), x0m, n, method="quad", radius=mpmath.mpf(radius))
            b = (x0m ** n) * ((-1) ** n) * d / mpmath.factorial(n)
            out.append(mpmath.mpf(b.real) if isinstance(b, mpmath.mpc) else b)
    return out


def enclosure_strings(encl: List[Interval], digits: int = 30) -> List[Tuple[str, str]]:
    return [(mpmath.nstr(e.lo, digits), mpmath.nstr(e.hi, digits)) for e in encl]
