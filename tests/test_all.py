import json
from fractions import Fraction as F
from pathlib import Path

import pytest

from hausdorff_certificates import (
    Interval,
    certify_exact,
    certify_interval,
    cm_violations,
    coercivity_bundle_exact,
    complete_monotonicity_table,
    exact_psd,
    hankel_H,
    hankel_L,
    hankel_L_theta,
    hausdorff_pair_exact,
    lambda_min_lower_bound_from_theta,
    moments_from_atoms,
    moments_lebesgue,
    quadratic_form,
    resolvent_trace_moments,
    theta_for_coercivity,
)
from hausdorff_certificates.intervals import interval_cholesky_pd, certified_negative_rayleigh
from hausdorff_certificates.manifest import (
    EXACT_TIER,
    ZETA_TIER,
    format_manifest_digest,
    load_manifest_digest,
    main as manifest_main,
)
from hausdorff_certificates.verify import verify_obj


# ----------------------------------------------------------- exact engine

def test_exact_psd_identity_and_hilbert():
    res = exact_psd([[F(1), F(0)], [F(0), F(2)]])
    assert res.is_psd and res.rank == 2
    # Hilbert 6x6 (moments of Lebesgue) is PD
    H = hankel_H(moments_lebesgue(10), 5)
    res = exact_psd(H)
    assert res.is_psd and res.rank == 6
    assert res.verify(H)


def test_exact_not_psd_witness_is_exact():
    A = [[F(1), F(2)], [F(2), F(1)]]  # eigenvalues 3, -1
    res = exact_psd(A)
    assert not res.is_psd
    assert quadratic_form(A, res.witness) == res.witness_value < 0


def test_exact_semidefinite_rank_deficient():
    # rank-1 PSD: outer product of (1,2,3)
    v = [F(1), F(2), F(3)]
    A = [[v[i] * v[j] for j in range(3)] for i in range(3)]
    res = exact_psd(A)
    assert res.is_psd and res.rank == 1


def test_exact_zero_diagonal_indefinite():
    A = [[F(0), F(1)], [F(1), F(0)]]
    res = exact_psd(A)
    assert not res.is_psd and res.witness_value < 0


def test_exact_zero_diag_positive_elsewhere():
    # diag has a zero but matrix is PSD only if that row is zero
    A = [[F(0), F(0), F(0)], [F(0), F(2), F(1)], [F(0), F(1), F(2)]]
    assert exact_psd(A).is_psd
    B = [[F(0), F(1), F(0)], [F(1), F(2), F(1)], [F(0), F(1), F(2)]]
    res = exact_psd(B)
    assert not res.is_psd and quadratic_form(B, res.witness) < 0


# ----------------------------------------------------- moments & builders

def test_atoms_and_support_test():
    # mu = 1/2 d_{1/3} + 1/2 d_{2/3}: Hausdorff on [0,1], support max = 2/3
    b = moments_from_atoms([(F(1, 2), F(1, 3)), (F(1, 2), F(2, 3))], 12)
    assert exact_psd(hankel_H(b, 5)).is_psd
    assert exact_psd(hankel_L(b, 5)).is_psd
    # true statement theta = 3/4 >= 2/3
    assert exact_psd(hankel_L_theta(b, 5, F(3, 4))).is_psd
    # false statement theta = 1/2 < 2/3 -> refuted with exact witness
    res = exact_psd(hankel_L_theta(b, 5, F(1, 2)))
    assert not res.is_psd


def test_cm_table_detects_corruption():
    b = moments_lebesgue(12)
    assert cm_violations(complete_monotonicity_table(b)) == []
    bad = list(b)
    bad[5] += F(1, 100)  # break complete monotonicity
    viols = cm_violations(complete_monotonicity_table(bad))
    assert viols, "corruption not detected by CM table"


def test_resolvent_trace_moments_match_eigen_decomposition():
    # K = diag(1, 2, 5): closed-form check of M_n and b_n
    K = [[F(1), F(0), F(0)], [F(0), F(2), F(0)], [F(0), F(0), F(5)]]
    x0 = F(1)
    M, b = resolvent_trace_moments(K, x0, 8)
    for n in range(9):
        expected = sum(F(1, (lam + 1) ** (n + 1)) for lam in (1, 2, 5))
        assert M[n] == expected
        assert b[n] == x0**n * expected
    # coercivity: lambda_min = 1 -> theta = 1/2 passes, claiming c = 2 fails
    th_true = theta_for_coercivity(F(1), x0)
    assert exact_psd(hankel_L_theta(b, 3, th_true)).is_psd
    th_false = theta_for_coercivity(F(2), x0)
    assert not exact_psd(hankel_L_theta(b, 3, th_false)).is_psd
    assert lambda_min_lower_bound_from_theta(th_true, x0) == F(1)


def test_laplacian_coercivity_bundle():
    # path-graph Laplacian + identity: tridiag(-1, 3, -1), 5x5; lambda_min = 3 - 2cos(pi/6) = 3 - sqrt(3) ~ 1.268
    n = 5
    K = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        K[i][i] = F(3)
        if i + 1 < n:
            K[i][i + 1] = K[i + 1][i] = F(-1)
    _, b = resolvent_trace_moments(K, F(1), 12)
    certs = coercivity_bundle_exact("lap", b, 5, theta_for_coercivity(F(1), F(1)))
    assert all(c.verdict == "PSD_CERTIFIED" for c in certs)  # claim c=1 true
    bad = coercivity_bundle_exact("lap_bad", b, 5, theta_for_coercivity(F(2), F(1)))
    assert bad[1].verdict == "NOT_PSD_CERTIFIED"  # claim c=2 false


# ------------------------------------------------------------- intervals

def _encl(b, rad=F(1, 10**25)):
    return [Interval(float(x - rad), float(x + rad)) for x in b]


def test_interval_cholesky_certifies_hilbert():
    import mpmath

    b = moments_lebesgue(14)
    encl = []
    for x in b:
        with mpmath.workdps(50):
            m = mpmath.mpf(x.numerator) / x.denominator
        encl.append(Interval(m - mpmath.mpf("1e-30"), m + mpmath.mpf("1e-30")))
    ok, pivots = interval_cholesky_pd(hankel_H(encl, 6), dps=60)
    assert ok and all(p > 0 for p in pivots)


def test_interval_refutation_with_witness():
    import mpmath

    # enclosure around an indefinite matrix
    A = [
        [Interval(mpmath.mpf("1"), mpmath.mpf("1")), Interval(mpmath.mpf("2"), mpmath.mpf("2"))],
        [Interval(mpmath.mpf("2"), mpmath.mpf("2")), Interval(mpmath.mpf("1"), mpmath.mpf("1"))],
    ]
    ok, _ = interval_cholesky_pd(A)
    assert not ok
    neg, bounds = certified_negative_rayleigh(A, [F(1), F(-1)])
    assert neg and float(bounds[1]) < 0


# ------------------------------------------------- certificates roundtrip

def test_certificate_roundtrip_exact(tmp_path):
    b = moments_from_atoms([(F(1, 2), F(1, 3)), (F(1, 2), F(2, 3))], 12)
    for cert in hausdorff_pair_exact("two_atoms", b, 5):
        obj = json.loads(cert.to_json())
        ok, msg = verify_obj(obj)
        assert ok, msg
    bad = certify_exact("bad_theta", b, "hankel_L_theta", 5, theta=F(1, 2))
    obj = json.loads(bad.to_json())
    ok, msg = verify_obj(obj)
    assert ok and obj["verdict"] == "NOT_PSD_CERTIFIED", msg


def test_certificate_tamper_detected():
    b = moments_lebesgue(12)
    cert = certify_exact("hilbert", b, "hankel_H", 5)
    obj = json.loads(cert.to_json())
    obj["moments"]["values"][3] = "1/3"  # tamper with a moment
    ok, _ = verify_obj(obj)
    assert not ok


def test_certificate_roundtrip_interval():
    import mpmath

    b = moments_lebesgue(14)
    encl = []
    for x in b:
        with mpmath.workdps(50):
            m = mpmath.mpf(x.numerator) / x.denominator
        encl.append(Interval(m - mpmath.mpf("1e-28"), m + mpmath.mpf("1e-28")))
    cert = certify_interval("hilbert_iv", encl, "hankel_H", 6, dps=60)
    assert cert.verdict == "PSD_CERTIFIED"
    obj = json.loads(cert.to_json())
    ok, msg = verify_obj(obj)
    assert ok, msg


def test_determinism_same_bytes():
    b = moments_lebesgue(12)
    c1 = certify_exact("h", b, "hankel_H", 5).to_json()
    c2 = certify_exact("h", b, "hankel_H", 5).to_json()
    assert c1 == c2


def test_manifest_digest_covers_artifact_roles(capsys):
    rows = load_manifest_digest(Path("artifacts/manifest.json"))
    by_file = {row.file: row for row in rows}

    assert len(rows) == 14
    assert by_file["hilbert_lebesgue__H.cert.json"].backend == "exact-rational"
    assert by_file["hilbert_lebesgue__H.cert.json"].verdict == "PSD_CERTIFIED"
    assert by_file["hilbert_lebesgue__H.cert.json"].tier == EXACT_TIER
    assert by_file["zeta_x0_1_trunc__H_N14.cert.json"].backend == "interval"
    assert by_file["zeta_x0_1_trunc__H_N14.cert.json"].tier == ZETA_TIER
    assert by_file["zeta_x0_1_derivative_crosscheck.json"].backend == "crosscheck"
    assert by_file["zeta_x0_1_derivative_crosscheck.json"].verdict == "all_inside=True"
    assert all(row.sha256_ok for row in rows)

    text = format_manifest_digest(rows)
    assert "laplacian5_coercivity_fail.cert.json" in text
    assert "demonstration-tier; gates V-A,V-B" in text

    assert manifest_main(["artifacts/manifest.json"]) == 0
    out = capsys.readouterr().out
    assert "file" in out and "sha256" in out


def test_manifest_digest_reports_hash_failures(tmp_path, capsys):
    src = Path("artifacts")
    copied = tmp_path / "artifacts"
    copied.mkdir()
    for artifact in src.iterdir():
        copied.joinpath(artifact.name).write_bytes(artifact.read_bytes())

    target = copied / "hilbert_lebesgue__H.cert.json"
    obj = json.loads(target.read_text(encoding="utf-8"))
    obj["name"] = "tampered_hilbert_lebesgue__H"
    target.write_text(json.dumps(obj, sort_keys=True, indent=1) + "\n", encoding="utf-8")

    assert manifest_main([str(copied / "manifest.json")]) == 1
    out = capsys.readouterr().out
    assert "hilbert_lebesgue__H.cert.json" in out
    assert "FAIL" in out


# ------------------------------------------------------------- zeta smoke

@pytest.mark.slow
def test_zeta_pipeline_smoke():
    from hausdorff_certificates.zeta import (
        compute_zero_ordinates,
        zeta_moment_enclosures,
    )

    gammas = compute_zero_ordinates(12, dps=30)
    assert abs(float(gammas[0]) - 14.134725) < 1e-5
    encl = zeta_moment_enclosures(gammas, x0=1, n_max=8, dps=50, inflate_exp=20)
    # enclosures are ordered and positive
    for e in encl:
        assert 0 < e.lo <= e.hi
    # b_0 = sum 1/(gamma^2+1): with 12 zeros ~ 0.023; tail keeps hi modest
    assert float(encl[0].lo) > 0.012 and float(encl[0].hi) < 0.06
    cert = certify_interval("zeta_smoke", encl, "hankel_H", 3, dps=60)
    assert cert.verdict in ("PSD_CERTIFIED", "INCONCLUSIVE")
