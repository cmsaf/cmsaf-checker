#!/usr/bin/env python3
"""
validate_xml.py — validate CM SAF metadata standard XML files against the DTD.

Usage:
    python3 validate_xml.py [options] file [file ...]

Each XML file is validated against the DTD.  The DTD is located by searching,
in order:
  1. The path given by --dtd (if supplied)
  2. The same directory as each XML file being validated
  3. The directory containing this script

Exit code: 0 if all files pass, 1 if any file fails or the DTD cannot be found.

Requires: lxml  (pip install lxml)
"""

import argparse
import os
import sys

try:
    from lxml import etree
except ImportError:
    print("ERROR: lxml is required.  Install it with:  pip install lxml", file=sys.stderr)
    sys.exit(1)


DTD_FILENAME = "cmsaf_metadata_standard.dtd"


def _find_dtd(xml_path: str, explicit: str | None) -> str | None:
    """Return the path to the DTD, or None if it cannot be found.

    Search order:
      1. explicit path supplied via --dtd
      2. directory of the XML file being validated
      3. directory of this script
    """
    candidates = []
    if explicit:
        candidates.append(explicit)
    candidates.append(os.path.join(os.path.dirname(os.path.abspath(xml_path)), DTD_FILENAME))
    candidates.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), DTD_FILENAME))

    for path in candidates:
        if os.path.isfile(path):
            return path
    return None


def validate_file(xml_path: str, dtd_path: str) -> list[str]:
    """Validate *xml_path* against the DTD at *dtd_path*.

    Returns a list of error message strings (empty list means the file is valid).
    """
    try:
        dtd = etree.DTD(dtd_path)
    except etree.DTDParseError as exc:
        return [f"DTD parse error: {exc}"]

    try:
        tree = etree.parse(xml_path)
    except etree.XMLSyntaxError as exc:
        return [f"XML syntax error: {exc}"]

    if not dtd.validate(tree):
        return [str(e) for e in dtd.error_log.filter_from_errors()]
    return []


def main() -> int:
    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(
        prog="validate_xml.py",
        description="Validate CM SAF metadata standard XML files against the DTD.",
    )
    parser.add_argument(
        "--dtd",
        metavar="PATH",
        default=None,
        help=(
            f"path to '{DTD_FILENAME}'. "
            "If omitted the DTD is searched in each XML file's own directory, "
            f"then in '{script_dir}'."
        ),
    )
    parser.add_argument("files", nargs="+", metavar="file", help="XML file(s) to validate")
    args = parser.parse_args()

    overall_rc = 0

    for xml_path in args.files:
        if not os.path.isfile(xml_path):
            print(f"MISSING  {xml_path}  (file not found)")
            overall_rc = 1
            continue

        dtd_path = _find_dtd(xml_path, args.dtd)
        if dtd_path is None:
            print(
                f"ERROR    {xml_path}  "
                f"(cannot find '{DTD_FILENAME}'; use --dtd to specify its location)"
            )
            overall_rc = 1
            continue

        errors = validate_file(xml_path, dtd_path)

        if errors:
            print(f"INVALID  {xml_path}")
            for msg in errors:
                print(f"         {msg}")
            overall_rc = 1
        else:
            print(f"OK       {xml_path}")
            print(f"         (validated against '{dtd_path}')")

    return overall_rc


if __name__ == "__main__":
    sys.exit(main())
