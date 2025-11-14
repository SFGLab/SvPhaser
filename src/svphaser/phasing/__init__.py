# src/svphaser/phasing/__init__.py
"""Public API for svphaser.phasing."""

from __future__ import annotations

import logging

from .algorithms import classify_haplotype, phasing_gq
from .io import phase_vcf
from .types import WorkerOpts

__all__ = [
    "phase_vcf",
    "classify_haplotype",
    "phasing_gq",
    "WorkerOpts",
]

# Library logging: don't emit anything unless the app configures it.
logging.getLogger(__name__).addHandler(logging.NullHandler())
