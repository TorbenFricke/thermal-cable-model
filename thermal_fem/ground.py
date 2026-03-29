"""Seasonal ground and ambient temperature models.

The Kasuda equation models the undisturbed ground temperature as a function
of depth and day-of-year, accounting for seasonal thermal wave penetration.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class GroundTemperatureModel(ABC):
    """Base class for ground / ambient temperature models."""

    @abstractmethod
    def temperature(self, depth: float, time_s: float) -> float:
        """Return ground temperature [°C] at *depth* [m] and *time_s* [s]."""

    def temperature_array(
        self, depth: float, times_s: np.ndarray
    ) -> np.ndarray:
        """Vectorised convenience wrapper."""
        return np.array([self.temperature(depth, t) for t in times_s])


class ConstantGroundTemperature(GroundTemperatureModel):
    """Spatially and temporally uniform ground temperature."""

    def __init__(self, temperature_c: float = 15.0):
        self.temp = temperature_c

    def temperature(self, depth: float, time_s: float) -> float:
        return self.temp


class KasudaModel(GroundTemperatureModel):
    r"""Kasuda & Archenbach (1965) undisturbed ground temperature.

    .. math::

        T(z,t) = T_\text{mean} - T_\text{amp} \exp\!\Bigl(
            -z \sqrt{\frac{\pi}{P\,\alpha}}\Bigr)
            \cos\!\Bigl(\frac{2\pi}{P}\bigl(t - t_0
            - \tfrac{z}{2}\sqrt{\frac{P}{\pi\,\alpha}}\bigr)\Bigr)

    Parameters
    ----------
    mean_surface_temp : float
        Annual mean surface temperature [°C].
    annual_amplitude : float
        Half the peak-to-peak annual surface temperature swing [°C].
    day_of_min_surface_temp : float
        Day of year at which the surface temperature is at its minimum (e.g. ~35
        for central Europe ≈ early February).
    soil_diffusivity : float
        Thermal diffusivity α of the soil [m²/s].  Typical range 0.3–0.8 × 10⁻⁶.
    """

    SECONDS_PER_DAY = 86_400.0
    DAYS_PER_YEAR = 365.25

    def __init__(
        self,
        mean_surface_temp: float = 12.0,
        annual_amplitude: float = 10.0,
        day_of_min_surface_temp: float = 35.0,
        soil_diffusivity: float = 0.5e-6,
    ):
        self.T_mean = mean_surface_temp
        self.T_amp = annual_amplitude
        self.t0_days = day_of_min_surface_temp
        self.alpha = soil_diffusivity
        self._period_s = self.DAYS_PER_YEAR * self.SECONDS_PER_DAY

    def temperature(self, depth: float, time_s: float) -> float:
        P = self._period_s
        z = abs(depth)
        omega = 2.0 * np.pi / P
        decay = np.exp(-z * np.sqrt(np.pi / (P * self.alpha)))
        phase_lag = z / 2.0 * np.sqrt(P / (np.pi * self.alpha))
        t0_s = self.t0_days * self.SECONDS_PER_DAY
        return float(
            self.T_mean
            - self.T_amp * decay * np.cos(omega * (time_s - t0_s - phase_lag))
        )

    def temperature_array(
        self, depth: float, times_s: np.ndarray
    ) -> np.ndarray:
        P = self._period_s
        z = abs(depth)
        omega = 2.0 * np.pi / P
        decay = np.exp(-z * np.sqrt(np.pi / (P * self.alpha)))
        phase_lag = z / 2.0 * np.sqrt(P / (np.pi * self.alpha))
        t0_s = self.t0_days * self.SECONDS_PER_DAY
        return self.T_mean - self.T_amp * decay * np.cos(
            omega * (times_s - t0_s - phase_lag)
        )
