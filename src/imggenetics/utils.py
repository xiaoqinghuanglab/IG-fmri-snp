"""Small shared utilities.

These helpers are used by both the CLI and the modeling code.
They intentionally avoid any project-specific I/O so they remain
safe to import from anywhere.
"""

from __future__ import annotations

import logging
import re


def setup_logging(verbose: bool = False) -> None:
    """Configure simple console logging.

    The pipeline is often run in batch jobs, so we default to WARNING
    to reduce clutter and only enable INFO messages with --verbose.
    """

    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(level=level, format="%(asctime)s | %(levelname)s | %(message)s")


def clean_column_name(name: str) -> str:
    """Make a column name safe for statsmodels formula strings.

    Statsmodels formulas treat many characters (spaces, '-', ':', etc.)
    as operators. The original notebook used a simple sanitizer; we keep
    the same logic here for backward compatibility.

    Steps:
      1) replace non-alphanumeric/_ with '_'
      2) collapse multiple underscores
      3) strip leading/trailing underscores
      4) if the name starts with a digit, prefix '_' (valid python identifier)
    """

    cleaned = re.sub(r"[^a-zA-Z0-9_]", "_", str(name))
    cleaned = re.sub(r"_{2,}", "_", cleaned)
    cleaned = cleaned.strip("_")
    if cleaned and cleaned[0].isdigit():
        cleaned = "_" + cleaned
    return cleaned


def normalize_subtype(x: str) -> str:
    """Normalize subtype labels to a consistent short form.

    Your datasets may contain multiple spellings for the same subtype.
    We map common values into the three default analysis labels:
      - Control
      - TypAD
      - AsymAD

    We also support a LowNFT label to make it easy to exclude if desired.
    Unknown labels are returned unchanged.
    """

    if x is None:
        return ""

    s = str(x).strip()
    # Normalize separators and case for matching.
    s_low = s.lower().replace(" ", "").replace("-", "").replace("_", "")

    if s_low in {"cn", "control", "controls", "cognitivenormal", "normal"}:
        return "Control"
    if s_low in {"typad", "typicalad", "typical"}:
        return "TypAD"
    if s_low in {"asymad", "asymptomaticad", "asym"}:
        return "AsymAD"

    # Some label-transfer pipelines output a separate LowNFT category.
    if s_low in {"lownft", "lownftad", "lownftadni", "lownftadcontrol", "lownftadtypicalad"} or "lownft" in s_low:
        return "LowNFT"

    return s
