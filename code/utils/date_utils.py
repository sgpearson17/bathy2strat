"""Date utilities for interval-based analyses."""

from __future__ import annotations

import numpy as np
import pandas as pd


def interval_highlight_mask(start_dates, end_dates, highlight_dates):
    """Return a boolean mask for intervals containing any highlight date.

    Args:
        start_dates: Iterable of interval start dates.
        end_dates: Iterable of interval end dates.
        highlight_dates: Iterable of dates to flag (strings or datetime-like).
    """
    if not highlight_dates:
        return np.zeros(len(start_dates), dtype=bool)

    start = pd.to_datetime(start_dates, errors="coerce").dt.normalize()
    end = pd.to_datetime(end_dates, errors="coerce").dt.normalize()
    highlights = pd.to_datetime(list(highlight_dates), errors="coerce")
    highlights = [d.normalize() for d in highlights if pd.notna(d)]

    if not highlights:
        return np.zeros(len(start), dtype=bool)

    mask = np.zeros(len(start), dtype=bool)
    for date in highlights:
        mask |= (start <= date) & (end >= date)

    return mask
