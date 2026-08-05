"""High-level certification API and the on-disk certificate format.

A *certificate* is a self-contained JSON document that an independent
checker (:mod:`hausdorff_certificates.verify`) can re-validate without
trusting the generator:

* exact backend: the LDL^T factors (PSD) or a rational witness vector with
  a strictly negative quadratic form (NOT PSD) -- both re-checkable in
  exact arithmetic from the stored moments;
* interval backend: the enclosure itself plus either the successful
  interval-Cholesky parameters (re-run to check) or a rational witness
  whose interval Rayleigh quotient has negative upper bound.

Certificates are deterministic: sorted keys, no timestamps, fixed number
formatting.  ``sha256`` of the file is the reproducibility anchor.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, List, Optional, Sequence

import mpmath

from . import hankel
from .intervals import (
    Interval,
    certified_negative_rayleigh,
    interval_cholesky_pd,
    propose_rational_witness,
)
from .rational import exact_psd

FORMAT = "hausdorff-certificate/1"

VERDICT_PSD = "PSD_CERTIFIED"
VERDICT_NOT_PSD = "NOT_PSD_CERTIFIED"
VERDICT_INCONCLUSIVE = "INCONCLUSIVE"


def _fr(s) -> str:
    return str(Fraction(s))


@dataclass
class Certificate:
    name: str
    backend: str  # "exact-rational" | "interval"
    matrix_kind: str  # "hankel_H" | "hankel_L" | "hankel_L_theta" | "shifted_hankel_S"
    N: int
    verdict: str
    moments: Dict
    certificate: Dict
    theta: Optional[str] = None
    provenance: Dict = field(default_factory=dict)
    cm_check: Optional[Dict] = None

    def to_obj(self) -> Dict:
        obj = {
            "format": FORMAT,
            "name": self.name,
            "backend": self.backend,
            "matrix": {"kind": self.matrix_kind, "N": self.N},
            "moments": self.moments,
            "verdict": self.verdict,
            "certificate": self.certificate,
            "provenance": self.provenance,
        }
        if self.theta is not None:
            obj["matrix"]["theta"] = self.theta
        if self.cm_check is not None:
            obj["cm_check"] = self.cm_check
        return obj

    def to_json(self) -> str:
        return json.dumps(self.to_obj(), sort_keys=True, indent=1) + "\n"

    def save(self, path: str) -> str:
        data = self.to_json()
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(data)
        return hashlib.sha256(data.encode("utf-8")).hexdigest()


# ----------------------------------------------------------------------
# matrix construction shared by generator and verifier
# ----------------------------------------------------------------------

def build_matrix(kind: str, b: Sequence, N: int, theta=None):
    if kind == "hankel_H":
        return hankel.hankel_H(b, N)
    if kind == "hankel_L":
        return hankel.hankel_L(b, N)
    if kind == "hankel_L_theta":
        if theta is None:
            raise ValueError("theta required for hankel_L_theta")
        return hankel.hankel_L_theta(b, N, Fraction(theta))
    if kind == "shifted_hankel_S":
        return hankel.shifted_hankel_S(b, N)
    raise ValueError(f"unknown matrix kind {kind!r}")


def moments_payload_exact(b: Sequence[Fraction]) -> Dict:
    return {"type": "rational", "values": [_fr(x) for x in b]}


def moments_payload_interval(encl: Sequence[Interval], digits: int = 70) -> Dict:
    """Serialize enclosures as decimal strings, rounded OUTWARD so that the
    parsed interval always contains the original one."""
    from .intervals import workdps_iv

    vals = []
    with workdps_iv(digits + 20):
        for e in encl:
            pad_lo = abs(e.lo) * mpmath.mpf(10) ** (2 - digits) + mpmath.mpf(10) ** (-4 * digits)
            pad_hi = abs(e.hi) * mpmath.mpf(10) ** (2 - digits) + mpmath.mpf(10) ** (-4 * digits)
            vals.append(
                [mpmath.nstr(e.lo - pad_lo, digits), mpmath.nstr(e.hi + pad_hi, digits)]
            )
    return {"type": "interval", "values": vals, "digits": digits}


def parse_moments(payload: Dict):
    if payload["type"] == "rational":
        return [Fraction(s) for s in payload["values"]]
    if payload["type"] == "interval":
        from .intervals import workdps_iv

        digits = int(payload.get("digits", 30))
        out = []
        with workdps_iv(digits + 20):
            for lo, hi in payload["values"]:
                out.append(Interval(mpmath.mpf(lo), mpmath.mpf(hi)))
        return out
    raise ValueError("unknown moments payload type")


# ----------------------------------------------------------------------
# exact backend
# ----------------------------------------------------------------------

def certify_exact(
    name: str,
    b: Sequence[Fraction],
    matrix_kind: str,
    N: int,
    theta=None,
    provenance: Optional[Dict] = None,
    with_cm: bool = True,
) -> Certificate:
    b = [Fraction(x) for x in b]
    A = build_matrix(matrix_kind, b, N, theta)
    res = exact_psd(A)
    if res.is_psd:
        verdict = VERDICT_PSD
        cert = {
            "type": "ldlt",
            "perm": res.perm,
            "D": [_fr(d) for d in res.D],
            "L": [[_fr(x) for x in row] for row in res.L],
            "rank": res.rank,
        }
    else:
        verdict = VERDICT_NOT_PSD
        cert = {
            "type": "negative_witness",
            "v": [_fr(x) for x in res.witness],
            "value": _fr(res.witness_value),
        }
    cm = None
    if with_cm:
        table = hankel.complete_monotonicity_table(b)
        viol = hankel.cm_violations(table)
        cm = {"max_order": len(b) - 1, "violations": [[k, n] for (k, n) in viol]}
    return Certificate(
        name=name,
        backend="exact-rational",
        matrix_kind=matrix_kind,
        N=N,
        verdict=verdict,
        moments=moments_payload_exact(b),
        certificate=cert,
        theta=_fr(theta) if theta is not None else None,
        provenance=provenance or {},
        cm_check=cm,
    )


# ----------------------------------------------------------------------
# interval backend
# ----------------------------------------------------------------------

def certify_interval(
    name: str,
    encl: Sequence[Interval],
    matrix_kind: str,
    N: int,
    theta=None,
    dps: int = 60,
    provenance: Optional[Dict] = None,
    with_cm: bool = True,
    digits: int = 70,
) -> Certificate:
    # what you store is what you certify: serialize outward, reparse, certify
    payload = moments_payload_interval(encl, digits=digits)
    encl = parse_moments(payload)
    A = build_matrix(matrix_kind, list(encl), N, theta)
    ok, pivots = interval_cholesky_pd(A, dps=dps)
    if ok:
        verdict = VERDICT_PSD
        cert = {
            "type": "interval_cholesky",
            "dps": dps,
            "pivot_lower_bounds": [mpmath.nstr(p, 20) for p in pivots],
        }
    else:
        v = propose_rational_witness(A)
        neg, bounds = certified_negative_rayleigh(A, v, dps=dps)
        if neg:
            verdict = VERDICT_NOT_PSD
            cert = {
                "type": "interval_negative_witness",
                "dps": dps,
                "v": [_fr(x) for x in v],
                "value_interval": list(bounds),
            }
        else:
            verdict = VERDICT_INCONCLUSIVE
            cert = {
                "type": "inconclusive",
                "dps": dps,
                "note": "interval Cholesky failed and no certified negative "
                "Rayleigh witness was found; enclosure may be too wide",
                "tried_witness_value_interval": list(bounds),
            }
    cm = None
    if with_cm:
        table = hankel.complete_monotonicity_table(list(encl))
        viol = hankel.cm_violations(table)
        cm = {"max_order": len(encl) - 1, "violations": [[k, n] for (k, n) in viol]}
    return Certificate(
        name=name,
        backend="interval",
        matrix_kind=matrix_kind,
        N=N,
        verdict=verdict,
        moments=payload,
        certificate=cert,
        theta=_fr(theta) if theta is not None else None,
        provenance=provenance or {},
        cm_check=cm,
    )


# ----------------------------------------------------------------------
# convenience: necessary finite Hausdorff and coercivity screens
# ----------------------------------------------------------------------

def hausdorff_pair_exact(name: str, b, N: int, provenance=None) -> List[Certificate]:
    return [
        certify_exact(f"{name}__H", b, "hankel_H", N, provenance=provenance),
        certify_exact(f"{name}__L", b, "hankel_L", N, provenance=provenance),
    ]


def coercivity_bundle_exact(name: str, b, N: int, theta, provenance=None) -> List[Certificate]:
    return [
        certify_exact(f"{name}__H", b, "hankel_H", N, provenance=provenance),
        certify_exact(f"{name}__Ltheta", b, "hankel_L_theta", N, theta=theta, provenance=provenance),
    ]


def hausdorff_pair_interval(name: str, encl, N: int, dps=60, provenance=None) -> List[Certificate]:
    return [
        certify_interval(f"{name}__H", encl, "hankel_H", N, dps=dps, provenance=provenance),
        certify_interval(f"{name}__L", encl, "hankel_L", N, dps=dps, provenance=provenance),
    ]
