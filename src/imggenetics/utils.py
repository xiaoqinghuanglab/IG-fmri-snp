from __future__ import annotations

import logging
import re

def setup_logging(verbose: bool = False) -> None:
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

def clean_column_name(name: str) -> str:
    """Clean a column name to be compatible with statsmodels formulas.

    Mirrors the original notebook logic:
    - replace non-alphanumeric/_ with '_'
    - collapse multiple underscores
    - strip leading/trailing underscores
    - prefix '_' if the name begins with a digit
    """
    cleaned = re.sub(r"[^a-zA-Z0-9_]", "_", str(name))
    cleaned = re.sub(r"_{2,}", "_", cleaned)
    cleaned = cleaned.strip("_")
    if cleaned and cleaned[0].isdigit():
        cleaned = "_" + cleaned
    return cleaned

def normalize_subtype(x: str) -> str:
    """Normalize common subtype spellings to short labels used in the pipeline.

    Output labels used by default:
      - Control
      - TypAD
      - AsymAD
      - LowNFT
    """
    if x is None:
        return ""
    s = str(x).strip()
    s_low = s.lower().replace(" ", "").replace("-", "").replace("_", "")
    if s_low in {"cn", "control", "controls", "cognitivenormal", "normal"}:
        return "Control"
    if s_low in {"typad", "typicalad", "typical"}:
        return "TypAD"
    if s_low in {"asymad", "asymptomaticad", "asym"}:
        return "AsymAD"
    if s_low in {"lownft", "lownftad", "lownftadni", "lownftadcontrol", "lownftadtypicalad"} or "lownft" in s_low:
        return "LowNFT"
    # If it's already short, keep it
    return s
