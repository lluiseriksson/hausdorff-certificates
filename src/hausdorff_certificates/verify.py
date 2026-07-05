"""Independent re-verification of stored certificates.

``python -m hausdorff_certificates.verify FILE [FILE ...]``

The verifier rebuilds the certificate matrix from the stored moments and
re-checks the stored evidence *without* re-running the decision procedure
that produced it (except for the interval-Cholesky case, where re-running
the factorisation at the stated precision *is* the check):

* ``ldlt``: exact reconstruction  P^T A P == L D L^T  and  D >= 0;
* ``negative_witness``: exact  v^T A v < 0;
* ``interval_cholesky``: re-run interval Cholesky on the enclosure;
* ``interval_negative_witness``: interval Rayleigh quotient upper bound < 0.

Exit code 0 iff every file verifies (INCONCLUSIVE verdicts verify vacuously
but are reported).
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from typing import Dict, Tuple

from .certify import FORMAT, build_matrix, parse_moments
from .intervals import certified_negative_rayleigh, interval_cholesky_pd
from .rational import ExactPSDResult, quadratic_form


def verify_obj(obj: Dict) -> Tuple[bool, str]:
    if obj.get("format") != FORMAT:
        return False, f"unknown format {obj.get('format')!r}"
    try:
        b = parse_moments(obj["moments"])
        kind = obj["matrix"]["kind"]
        N = obj["matrix"]["N"]
        theta = obj["matrix"].get("theta")
        A = build_matrix(kind, b, N, theta)
    except Exception as exc:  # noqa: BLE001 - verifier API reports malformed payloads
        return False, f"matrix rebuild failed: {exc}"
    cert = obj["certificate"]
    verdict = obj["verdict"]
    ctype = cert["type"]

    if ctype == "ldlt":
        if verdict != "PSD_CERTIFIED":
            return False, "ldlt evidence with non-PSD verdict"
        res = ExactPSDResult(
            True,
            N + 1,
            perm=list(cert["perm"]),
            L=[[Fraction(x) for x in row] for row in cert["L"]],
            D=[Fraction(x) for x in cert["D"]],
        )
        ok = res.verify(A)
        return ok, "LDL^T identity and D >= 0 re-checked exactly" if ok else "LDL^T re-check FAILED"

    if ctype == "negative_witness":
        if verdict != "NOT_PSD_CERTIFIED":
            return False, "witness evidence with non-refuting verdict"
        v = [Fraction(x) for x in cert["v"]]
        val = quadratic_form(A, v)
        ok = val < 0 and str(val) == cert["value"]
        return ok, (
            f"v^T A v = {val} < 0 re-checked exactly" if ok else "witness re-check FAILED"
        )

    if ctype == "interval_cholesky":
        if verdict != "PSD_CERTIFIED":
            return False, "interval-Cholesky evidence with non-PSD verdict"
        ok, _ = interval_cholesky_pd(A, dps=int(cert["dps"]))
        return ok, (
            "interval Cholesky re-run: every matrix in the enclosure is PD"
            if ok
            else "interval Cholesky re-run FAILED"
        )

    if ctype == "interval_negative_witness":
        if verdict != "NOT_PSD_CERTIFIED":
            return False, "interval witness with non-refuting verdict"
        v = [Fraction(x) for x in cert["v"]]
        neg, bounds = certified_negative_rayleigh(A, v, dps=int(cert["dps"]))
        return neg, (
            f"interval Rayleigh in {bounds}, upper bound < 0"
            if neg
            else "interval witness re-check FAILED"
        )

    if ctype == "inconclusive":
        return verdict == "INCONCLUSIVE", "inconclusive (nothing certified; nothing to check)"

    return False, f"unknown certificate type {ctype!r}"


def verify_file(path: str) -> Tuple[bool, str, str]:
    with open(path, "r", encoding="utf-8") as fh:
        obj = json.load(fh)
    ok, msg = verify_obj(obj)
    return ok, obj.get("verdict", "?"), msg


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(__doc__)
        return 2
    all_ok = True
    for path in argv:
        try:
            ok, verdict, msg = verify_file(path)
        except Exception as exc:  # noqa: BLE001 - report and fail
            ok, verdict, msg = False, "?", f"exception: {exc!r}"
        status = "OK " if ok else "FAIL"
        print(f"[{status}] {path}  verdict={verdict}  ({msg})")
        all_ok &= ok
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
