"""Module entry point.

Enables:
    python -m imggenetics ...
"""

from .cli import main


if __name__ == "__main__":
    raise SystemExit(main())
