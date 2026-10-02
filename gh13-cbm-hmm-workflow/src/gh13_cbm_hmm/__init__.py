"""Standalone profile-HMM screening and cutoff evaluation."""

from .scan import DomainHit, load_hmms, scan_fasta, scan_records
from .validation import (
    confusion_at_threshold,
    roc_auc,
    roc_curve,
    select_zero_observed_false_positive_cutoff,
    wilson_interval,
)

__all__ = [
    "DomainHit",
    "confusion_at_threshold",
    "load_hmms",
    "roc_auc",
    "roc_curve",
    "scan_fasta",
    "scan_records",
    "select_zero_observed_false_positive_cutoff",
    "wilson_interval",
]

__version__ = "0.1.0"
