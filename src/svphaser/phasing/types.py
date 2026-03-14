from dataclasses import dataclass
from typing import NamedTuple

SVKeyLegacy = tuple[str, int, str]
SVKey = tuple[str, int, str, int, str]
GQBin = tuple[int, str]


@dataclass(slots=True, frozen=True)
class WorkerOpts:
    min_support: int
    min_tagged_support: int
    major_delta: float
    equal_delta: float
    tie_to_hom_alt: bool
    support_mode: str
    bp_window: int
    dynamic_window: bool

    # NEW: strict DEL/INS size consistency filtering
    size_match_required: bool
    size_tol_abs: int
    size_tol_frac: float

    gq_bins: list[GQBin]


class CallTuple(NamedTuple):
    gt: str
    gq: int
    gq_label: str | None
