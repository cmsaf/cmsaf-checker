#!/usr/bin/env python3
"""
download_gcmd_keywords.py

Downloads GCMD keyword scheme CSV files from the NASA KMS REST API
for offline use in a keyword checker program.

Schemes downloaded:
  - Science Keywords
  - Platforms
  - Instruments
  - Providers

NOTE: The KMS API migrated to cmr.earthdata.nasa.gov in 2025.
The old gcmd.earthdata.nasa.gov base URL is no longer reliable.

The keyword version is read from the first line of each downloaded CSV
and embedded in the output filename, e.g. gcmd_platforms_v9.1.5.csv.
Passing a `version` query parameter to the scheme endpoint is currently NOT
supported and results in an "Invalid version parameter" error.
"""

import re
import requests
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Updated base URL as of 2025 migration
KMS_BASE_URL = "https://cmr.earthdata.nasa.gov/kms"

SCHEMES = {
    "sciencekeywords": "sciencekeywords",
    "platforms":       "platforms",
    "instruments":     "instruments",
    "providers":       "providers",
}

OUTPUT_DIR = Path("share")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_version_from_csv(content: bytes) -> str | None:
    """
    Extract the keyword version from the first line of a KMS CSV response.

    The KMS CSV header looks like:
        "Keyword Version: 23.5","Revision: 2026-03-13T11:13:40.099Z",...
    """
    try:
        first_line = content.split(b"\n")[0].decode("utf-8", errors="replace")
        match = re.search(r'Keyword Version:\s*([\d.]+)', first_line)
        if match:
            return match.group(1)   # e.g. "23.5"
    except Exception:
        pass
    return None


def download_scheme(
    session: requests.Session,
    scheme: str,
    base_name: str,
) -> tuple[str, str | None]:
    """Download one keyword scheme CSV.

    Returns a tuple of (status, version) where status is one of:
      'saved'   -- file was written to disk
      'skipped' -- file already existed and user chose not to overwrite
    Raises requests.exceptions.RequestException on network/HTTP failure.
    """
    url    = f"{KMS_BASE_URL}/concepts/concept_scheme/{scheme}"
    params = {"format": "csv"}

    print(f"  Fetching {scheme!r} ...")

    try:
        response = session.get(url, params=params, timeout=60)
        response.raise_for_status()
    except requests.exceptions.HTTPError as exc:
        print(
            f"  [ERROR] HTTP {exc.response.status_code} for scheme {scheme!r}: {exc}",
            file=sys.stderr,
        )
        raise
    except requests.exceptions.RequestException as exc:
        print(f"  [ERROR] Request failed for scheme {scheme!r}: {exc}", file=sys.stderr)
        raise

    version  = parse_version_from_csv(response.content)
    suffix   = f"_v{version}" if version else ""
    filename = f"{base_name}{suffix}.csv"
    out_path = OUTPUT_DIR / filename

    if out_path.exists():
        size_existing_kb = out_path.stat().st_size / 1024
        print(f"  [WARN] File already exists: {out_path}  ({size_existing_kb:.1f} KB)")
        try:
            answer = input(f"  Overwrite '{filename}'? [y/N] ").strip().lower()
        except EOFError:
            answer = "n"   # non-interactive environment: default to skip
        if answer != "y":
            print(f"  Skipped -> {out_path}")
            return "skipped", version

    out_path.write_bytes(response.content)
    size_kb = len(response.content) / 1024
    version_label = version or "unknown (upstream bug - check CSV header)"
    print(f"  Saved  -> {out_path}  ({size_kb:.1f} KB)  [version: {version_label}]")

    return "saved", version

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory : {OUTPUT_DIR.resolve()}\n")
    print(f"KMS base URL     : {KMS_BASE_URL}\n")

    res      = {"saved": 0, "skipped": 0, "failed": 0}
    versions = set()

    with requests.Session() as session:
        session.headers.update({"Accept": "text/csv"})

        for scheme, base_name in SCHEMES.items():
            try:
                status, v = download_scheme(session, scheme, base_name)
                res[status] += 1
                if v:
                    versions.add(v)
            except requests.exceptions.RequestException:
                res["failed"] += 1

    print()
    if versions:
        print(f"Keyword version(s) detected : {', '.join(sorted(versions))}")

    total = len(SCHEMES)
    print(f"\nSummary: {res['saved']} saved, {res['skipped']} skipped, "
          f"{res['failed']} failed  (out of {total} schemes)")

    if res["failed"] > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
