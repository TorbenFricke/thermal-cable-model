"""Thermal analysis of cables crossing at an angle.

When two cable circuits cross in the ground at different depths, the deeper
(or hotter) cable raises the temperature of the shallower one and vice-versa.
This module computes the mutual temperature rise at the crossing point
using the analytical line-source / image-source superposition method.

Reference
---------
CIGRÉ Technical Brochure 640 (2015) – *A Guide for Rating Calculations of
Insulated Cables*, §5.4 "Crossing cables".
"""

from __future__ import annotations

import math

import numpy as np
from scipy import integrate

from thermal_cable_model.cable import Cable
from thermal_cable_model.materials import ThermalMaterial


class CableCrossing:
    """Model for two cable systems crossing at an arbitrary angle.

    Parameters
    ----------
    cable_upper, cable_lower : Cable
        The two cables involved in the crossing.
    depth_upper, depth_lower : float
        Burial depth (to centre) of each cable [m].
    crossing_angle_deg : float
        Angle between the cable axes projected onto the horizontal plane [°].
        90° = perpendicular crossing, 0° = parallel (degenerate case).
    soil : ThermalMaterial
        Surrounding soil.
    """

    def __init__(
        self,
        cable_upper: Cable,
        cable_lower: Cable,
        depth_upper: float,
        depth_lower: float,
        crossing_angle_deg: float,
        soil: ThermalMaterial,
    ):
        self.cable_upper = cable_upper
        self.cable_lower = cable_lower
        self.depth_upper = depth_upper
        self.depth_lower = depth_lower
        self.angle_rad = math.radians(crossing_angle_deg)
        self.soil = soil

        if abs(self.angle_rad) < 1e-6:
            raise ValueError(
                "Crossing angle ≈ 0° (parallel) is not a crossing — "
                "use parallel cable modelling instead."
            )

    def temperature_rise_at_upper(
        self,
        heat_rate_lower: float,
        integration_half_length: float = 50.0,
    ) -> float:
        """Steady-state temperature rise at the upper cable centre-line due
        to the lower cable's heat emission [°C].

        Parameters
        ----------
        heat_rate_lower : float
            Total heat output of the lower cable per unit length [W/m].
        integration_half_length : float
            Half-length of the integration domain along the lower cable [m].
            50 m is sufficient for most practical cases.
        """
        return self._mutual_rise(
            heat_rate_lower,
            self.depth_upper,
            self.depth_lower,
            integration_half_length,
        )

    def temperature_rise_at_lower(
        self,
        heat_rate_upper: float,
        integration_half_length: float = 50.0,
    ) -> float:
        """Steady-state temperature rise at the lower cable due to the upper."""
        return self._mutual_rise(
            heat_rate_upper,
            self.depth_lower,
            self.depth_upper,
            integration_half_length,
        )

    def _mutual_rise(
        self,
        W: float,
        depth_target: float,
        depth_source: float,
        L: float,
    ) -> float:
        r"""Compute ΔT at target cable axis due to a line heat source.

        The source cable runs along direction :math:`\hat{s}` at depth
        *depth_source*.  At the crossing point (taken as the origin) the
        target cable is directly above/below.  The vertical separation is
        Δh = |depth_source − depth_target|.

        .. math::

            \Delta T = \frac{W}{4\pi\lambda}
                \int_{-L}^{L}
                \left(
                    \frac{1}{\sqrt{s^2 \sin^2\alpha + \Delta h^2}}
                  - \frac{1}{\sqrt{s^2 \sin^2\alpha + (h_1 + h_2)^2}}
                \right) ds

        The second term is the image correction for the ground surface.
        """
        lam = self.soil.thermal_conductivity
        sin_a = math.sin(self.angle_rad)
        dh = abs(depth_source - depth_target)
        h_sum = depth_source + depth_target

        def integrand(s: float) -> float:
            r_sq = (s * sin_a) ** 2
            term_real = 1.0 / math.sqrt(r_sq + dh ** 2)
            term_image = 1.0 / math.sqrt(r_sq + h_sum ** 2)
            return term_real - term_image

        result, _ = integrate.quad(integrand, -L, L, points=[0.0], limit=100)
        return W / (4.0 * math.pi * lam) * result

    # ── transient mutual temperature rise ────────────────────────────

    def transient_temperature_rise_at_upper(
        self,
        heat_rate_lower: float,
        time_s: float,
        integration_half_length: float = 50.0,
    ) -> float:
        """Transient temperature rise at the upper cable [°C].

        Uses the transient line-source Green's function for a semi-infinite
        medium (exponential integral formulation).
        """
        return self._transient_mutual_rise(
            heat_rate_lower,
            self.depth_upper,
            self.depth_lower,
            time_s,
            integration_half_length,
        )

    def transient_temperature_rise_at_lower(
        self,
        heat_rate_upper: float,
        time_s: float,
        integration_half_length: float = 50.0,
    ) -> float:
        return self._transient_mutual_rise(
            heat_rate_upper,
            self.depth_lower,
            self.depth_upper,
            time_s,
            integration_half_length,
        )

    def _transient_mutual_rise(
        self,
        W: float,
        depth_target: float,
        depth_source: float,
        t: float,
        L: float,
    ) -> float:
        r"""Transient mutual temperature rise using the exponential integral.

        .. math::

            \Delta T(t) = \frac{W}{4\pi\lambda}
                \int_{-L}^{L}
                \left[
                    \text{E}_1\!\left(\frac{r^2}{4\alpha t}\right)
                  - \text{E}_1\!\left(\frac{{r'}^2}{4\alpha t}\right)
                \right] ds

        where E₁ is the exponential integral, r and r' are distances
        to the source and its image respectively.
        """
        from scipy.special import exp1
        import warnings

        if t <= 0:
            return 0.0

        lam = self.soil.thermal_conductivity
        alpha = self.soil.thermal_diffusivity
        sin_a = math.sin(self.angle_rad)
        dh = abs(depth_source - depth_target)
        h_sum = depth_source + depth_target
        four_alpha_t = 4.0 * alpha * t

        # Adaptive integration limit: beyond this distance the integrand
        # contribution is negligible (thermal diffusion length scale).
        L_eff = min(L, max(5.0 * math.sqrt(four_alpha_t) / sin_a, 2.0))

        def integrand(s: float) -> float:
            r_sq = (s * sin_a) ** 2 + dh ** 2
            ri_sq = (s * sin_a) ** 2 + h_sum ** 2
            return float(exp1(r_sq / four_alpha_t) - exp1(ri_sq / four_alpha_t))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", integrate.IntegrationWarning)
            result, _ = integrate.quad(
                integrand, -L_eff, L_eff, points=[0.0], limit=200,
            )
        return W / (4.0 * math.pi * lam) * result


def crossing_derating_factor(
    cable: Cable,
    depth: float,
    delta_T_crossing: float,
    soil: ThermalMaterial,
) -> float:
    """Compute the derating factor due to a crossing-induced temperature rise.

    The cable's permissible current is reduced so that the total conductor
    temperature (self-heating + crossing contribution) does not exceed the
    maximum allowed temperature.

    Returns a factor in (0, 1] to multiply with the isolated cable rating.
    """
    from thermal_cable_model.thermal_network import (
        InternalThermalResistances,
        external_thermal_resistance,
    )

    tr = InternalThermalResistances.from_cable(cable)
    T4 = external_thermal_resistance(depth, cable.outer_radius, soil)
    n = cable.n_conductors
    lam1 = cable.loss_factor_sheath
    lam2 = cable.loss_factor_armour

    T_total = (
        tr.T1 / n
        + (1 + lam1) * tr.T2 / n
        + (1 + lam1 + lam2) * (tr.T3 / n + T4)
    )
    delta_T_max = cable.max_conductor_temp - 20.0  # assume 20 °C ambient

    if delta_T_crossing >= delta_T_max:
        return 0.0

    return math.sqrt((delta_T_max - delta_T_crossing) / delta_T_max)
