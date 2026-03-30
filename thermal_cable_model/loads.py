"""Dynamic electrical load (current) profiles.

Supports constant, cyclic, arbitrary time-series, and CSV import.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


class LoadProfile:
    """Piece-wise-linear current profile I(t).

    Parameters
    ----------
    times : array-like
        Monotonically increasing time stamps [s].
    currents : array-like
        RMS current values [A] at each time stamp.
    """

    def __init__(self, times: np.ndarray, currents: np.ndarray):
        self.times = np.asarray(times, dtype=float)
        self.currents = np.asarray(currents, dtype=float)
        if self.times.shape != self.currents.shape:
            raise ValueError("times and currents must have the same length")

    def current_at(self, t: float) -> float:
        """Linearly interpolated current at time *t* [s]."""
        return float(np.interp(t, self.times, self.currents))

    def current_array(self, times: np.ndarray) -> np.ndarray:
        """Vectorised interpolation."""
        return np.interp(times, self.times, self.currents)

    @property
    def duration(self) -> float:
        return float(self.times[-1] - self.times[0])

    # ── factory methods ──────────────────────────────────────────────

    @classmethod
    def constant(cls, current_a: float, duration_s: float) -> LoadProfile:
        """Constant current for a given duration."""
        return cls(
            times=np.array([0.0, duration_s]),
            currents=np.array([current_a, current_a]),
        )

    @classmethod
    def cyclic(
        cls,
        peak_current: float,
        base_current: float,
        period_s: float,
        duty_cycle: float,
        n_cycles: int,
    ) -> LoadProfile:
        """Rectangular on/off cyclic load.

        Parameters
        ----------
        peak_current, base_current : float
            Current during on and off phases [A].
        period_s : float
            Duration of one full cycle [s].
        duty_cycle : float
            Fraction of period at peak current (0–1).
        n_cycles : int
            Number of complete cycles.
        """
        t_on = period_s * duty_cycle
        times, currents = [0.0], [peak_current]
        for i in range(n_cycles):
            t_start = i * period_s
            times += [t_start + t_on - 1e-3, t_start + t_on, t_start + period_s - 1e-3]
            currents += [peak_current, base_current, base_current]
            if i < n_cycles - 1:
                times.append(t_start + period_s)
                currents.append(peak_current)
        return cls(times=np.array(times), currents=np.array(currents))

    @classmethod
    def daily_pattern(
        cls,
        hourly_currents: list[float],
        n_days: int = 1,
    ) -> LoadProfile:
        """Create a profile from 24 hourly current values repeated for *n_days*.

        Parameters
        ----------
        hourly_currents : list of 24 floats
            Current [A] for each hour of the day (00:00–23:00).
        n_days : int
            Number of days to repeat.
        """
        if len(hourly_currents) != 24:
            raise ValueError("Need exactly 24 hourly values")
        t_list, i_list = [], []
        for day in range(n_days):
            for h, amp in enumerate(hourly_currents):
                t_list.append((day * 24 + h) * 3600.0)
                i_list.append(amp)
        t_list.append(n_days * 24 * 3600.0)
        i_list.append(hourly_currents[0])
        return cls(times=np.array(t_list), currents=np.array(i_list))

    @classmethod
    def from_csv(
        cls,
        filepath: str | Path,
        time_col: int = 0,
        current_col: int = 1,
        delimiter: str = ",",
        skip_header: int = 1,
    ) -> LoadProfile:
        """Load a profile from a two-column CSV (time [s], current [A])."""
        data = np.loadtxt(filepath, delimiter=delimiter, skiprows=skip_header)
        return cls(times=data[:, time_col], currents=data[:, current_col])
