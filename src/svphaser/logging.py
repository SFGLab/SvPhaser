"""
svphaser.logging
================
One-liner that gives us colour-free, concise log messages on stderr.

Importing this module *once* anywhere in the program is enough.
"""

from __future__ import annotations

import logging
import sys
from typing import Literal

_Level = Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]


def init(level: _Level | int = "INFO") -> None:
    """Install a basic stderr handler – safe to call multiple times."""
    root = logging.getLogger()
    if root.handlers:  # already initialised
        return

    logging.basicConfig(
        level=level,
        stream=sys.stderr,
        format="%(levelname).1s | %(message)s",
        datefmt="%H:%M:%S",
    )


def get_logger(name: str) -> logging.Logger:
    """Convenience helper for module-level loggers."""
    return logging.getLogger(name)
