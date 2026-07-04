"""Digest and hash-check helpers for ``artifacts/manifest.json``.

``python -m hausdorff_certificates.manifest artifacts/manifest.json``
prints a deterministic, mother-facing summary of the committed artifact set:
file name, backend, verdict, and epistemic tier.  It also checks every listed
SHA-256 digest against the artifact content.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


MANIFEST_FORMAT = "hausdorff-certificates-manifest/1"
ZETA_TIER = "demonstration-tier; gates V-A,V-B"
EXACT_TIER = "exact-unconditional"


@dataclass(frozen=True)
class ManifestRow:
    file: str
    backend: str
    verdict: str
    tier: str
    sha256_ok: bool


def _digest(data: bytes) -> str:
    h = hashlib.sha256(data)
    return h.hexdigest()


def _sha256_candidates(path: Path) -> set[str]:
    with path.open("rb") as fh:
        data = fh.read()
    # The committed manifest is generated from LF JSON.  On Windows checkouts
    # with autocrlf enabled, the working tree may contain CRLF text while the
    # committed content is unchanged.
    return {_digest(data), _digest(data.replace(b"\r\n", b"\n"))}


def _tier_for(filename: str, obj: dict) -> str:
    provenance = obj.get("provenance", {})
    status = " ".join(
        str(provenance.get(key, "")) for key in ("epistemic_status", "mode", "pipeline")
    )
    if filename.startswith("zeta_") or "SOURCES.md gate" in status or "gate V-" in status:
        return ZETA_TIER
    return EXACT_TIER


def _backend_and_verdict(obj: dict) -> tuple[str, str]:
    if obj.get("format") == "hausdorff-certificate/1":
        return str(obj.get("backend", "?")), str(obj.get("verdict", "?"))
    if "all_inside" in obj:
        return "crosscheck", f"all_inside={obj['all_inside']}"
    return str(obj.get("backend", "?")), str(obj.get("verdict", "?"))


def load_manifest_digest(manifest_path: Path) -> list[ManifestRow]:
    with manifest_path.open("r", encoding="utf-8") as fh:
        manifest = json.load(fh)
    if manifest.get("format") != MANIFEST_FORMAT:
        raise ValueError(f"unknown manifest format {manifest.get('format')!r}")

    base = manifest_path.parent
    rows: list[ManifestRow] = []
    for filename, expected_hash in sorted(manifest["files"].items()):
        artifact_path = base / filename
        with artifact_path.open("r", encoding="utf-8") as fh:
            obj = json.load(fh)
        backend, verdict = _backend_and_verdict(obj)
        rows.append(
            ManifestRow(
                file=filename,
                backend=backend,
                verdict=verdict,
                tier=_tier_for(filename, obj),
                sha256_ok=expected_hash in _sha256_candidates(artifact_path),
            )
        )
    return rows


def format_manifest_digest(rows: Iterable[ManifestRow]) -> str:
    ordered = list(rows)
    widths = {
        "file": max(len("file"), *(len(row.file) for row in ordered)),
        "backend": max(len("backend"), *(len(row.backend) for row in ordered)),
        "verdict": max(len("verdict"), *(len(row.verdict) for row in ordered)),
        "tier": max(len("tier"), *(len(row.tier) for row in ordered)),
    }
    lines = [
        f"{'file':<{widths['file']}}  {'backend':<{widths['backend']}}  "
        f"{'verdict':<{widths['verdict']}}  {'tier':<{widths['tier']}}  sha256",
        f"{'-' * widths['file']}  {'-' * widths['backend']}  "
        f"{'-' * widths['verdict']}  {'-' * widths['tier']}  ------",
    ]
    for row in ordered:
        sha = "OK" if row.sha256_ok else "FAIL"
        lines.append(
            f"{row.file:<{widths['file']}}  {row.backend:<{widths['backend']}}  "
            f"{row.verdict:<{widths['verdict']}}  {row.tier:<{widths['tier']}}  {sha}"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", nargs="?", default="artifacts/manifest.json")
    args = parser.parse_args(argv)

    try:
        rows = load_manifest_digest(Path(args.manifest))
    except Exception as exc:  # noqa: BLE001 - CLI should report any manifest defect
        print(f"manifest digest failed: {exc}", file=sys.stderr)
        return 1
    print(format_manifest_digest(rows))
    return 0 if all(row.sha256_ok for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
