from __future__ import annotations

import re
import logging
from typing import Optional


def setup_logging(verbose: bool = False) -> None:
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def clean_column_name(name: str) -> str:
    """
    Clean a column name to be compatible with statsmodels formulas.

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
