"""High-level simulation interface.

``CableInstallation`` defines the physical arrangement (cables, positions,
soil, ground temperature).  ``ThermalSimulation`` runs the analysis and
returns results that can be queried or plotted.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from thermal_cable_model.cable import Cable
from thermal_cable_model.crossing import CableCrossing
from thermal_cable_model.ground import (
    ConstantGroundTemperature,
    GroundTemperatureModel,
)
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import SOIL_STANDARD, ThermalMaterial
from thermal_cable_model.solver import TransientResult, TransientSolver
from thermal_cable_model.thermal_network import CableThermalNetwork


@dataclass
class _CableEntry:
    cable: Cable
    x: float
    depth: float
    load: LoadProfile


class CableInstallation:
    """Describes a buried cable arrangement.

    Example
    -------
    >>> inst = CableInstallation(soil=SOIL_STANDARD)
    >>> inst.add_cable(cable_a, x=0.0, depth=1.2, load=load_a)
    >>> inst.add_cable(cable_b, x=0.3, depth=1.2, load=load_b)
    """

    def __init__(
        self,
        soil: ThermalMaterial = SOIL_STANDARD,
        ground_temp_model: GroundTemperatureModel | None = None,
    ):
        self.soil = soil
        self.ground_model = ground_temp_model or ConstantGroundTemperature(15.0)
        self._entries: list[_CableEntry] = []
        self._crossings: list[CableCrossing] = []

    def add_cable(
        self,
        cable: Cable,
        x: float,
        depth: float,
        load: LoadProfile,
    ) -> int:
        """Add a cable to the installation.

        Parameters
        ----------
        cable : Cable
        x : float
            Horizontal position [m] (arbitrary reference).
        depth : float
            Burial depth to cable centre [m] (positive downward).
        load : LoadProfile
            Time-varying current for this cable.

        Returns
        -------
        int
            Index of the added cable.
        """
        idx = len(self._entries)
        self._entries.append(_CableEntry(cable, x, depth, load))
        return idx

    def add_crossing(self, crossing: CableCrossing) -> None:
        """Register a cable crossing for derating analysis."""
        self._crossings.append(crossing)

    @property
    def n_cables(self) -> int:
        return len(self._entries)

    @property
    def cables(self) -> list[Cable]:
        return [e.cable for e in self._entries]

    @property
    def positions_x(self) -> list[float]:
        return [e.x for e in self._entries]

    @property
    def depths(self) -> list[float]:
        return [e.depth for e in self._entries]

    @property
    def loads(self) -> list[LoadProfile]:
        return [e.load for e in self._entries]

    @property
    def crossings(self) -> list[CableCrossing]:
        return list(self._crossings)


class ThermalSimulation:
    """Run steady-state or transient thermal analysis on an installation.

    Parameters
    ----------
    installation : CableInstallation
    """

    def __init__(self, installation: CableInstallation):
        self.inst = installation
        self.network = CableThermalNetwork(
            cables=installation.cables,
            positions_x=installation.positions_x,
            depths=installation.depths,
            soil=installation.soil,
        )

    def run_steady_state(
        self,
        time_s: float = 0.0,
    ) -> TransientResult:
        """Compute the steady-state temperature at a given time instant.

        The currents and ambient temperatures are evaluated at *time_s*.
        """
        currents = [lp.current_at(time_s) for lp in self.inst.loads]
        amb = [
            self.inst.ground_model.temperature(d, time_s)
            for d in self.inst.depths
        ]
        theta = self.network.steady_state(currents, amb)

        n = self.inst.n_cables
        result = TransientResult(
            times=np.array([time_s]),
            conductor_temps=np.array([self.network.get_conductor_temperatures(theta)]),
            insulation_temps=np.array([self.network.get_insulation_temperatures(theta)]),
            sheath_temps=np.array([self.network.get_sheath_temperatures(theta)]),
            armour_temps=np.array([self.network.get_armour_temperatures(theta)]),
            surface_temps=np.array([self.network.get_surface_temperatures(theta)]),
            soil_temps=np.array([self.network.get_soil_temperatures(theta)]),
            currents=np.array([currents]),
            ambient_temps=np.array([amb]),
            cable_names=[c.name for c in self.inst.cables],
        )

        # Apply crossing temperature rises if any
        if self.inst.crossings:
            for crossing in self.inst.crossings:
                W_upper = crossing.cable_upper.total_heat_per_length(
                    currents[0], result.conductor_temps[0, 0]
                )
                W_lower = crossing.cable_lower.total_heat_per_length(
                    currents[-1] if n > 1 else currents[0],
                    result.conductor_temps[0, -1] if n > 1 else result.conductor_temps[0, 0],
                )
                dT_upper = crossing.temperature_rise_at_upper(W_lower)
                dT_lower = crossing.temperature_rise_at_lower(W_upper)
                result.conductor_temps[0, 0] += dT_upper
                if n > 1:
                    result.conductor_temps[0, -1] += dT_lower

        return result

    def run_transient(
        self,
        dt: float = 60.0,
        duration: float | None = None,
        store_full_state: bool = False,
    ) -> TransientResult:
        """Run a transient simulation.

        Parameters
        ----------
        dt : float
            Time step [s].  Default 60 s.
        duration : float, optional
            Total simulation time [s].  Defaults to the longest load profile.
        store_full_state : bool
            Keep the full node-temperature array at every step.
        """
        if duration is None:
            duration = max(lp.duration for lp in self.inst.loads)

        solver = TransientSolver(
            network=self.network,
            load_profiles=self.inst.loads,
            ground_model=self.inst.ground_model,
        )

        result = solver.solve(
            dt=dt,
            duration=duration,
            store_full_state=store_full_state,
        )

        # Superimpose crossing temperature rises (added to conductor temps)
        if self.inst.crossings:
            n = self.inst.n_cables
            for step_idx in range(result.n_steps):
                t = result.times[step_idx]
                for crossing in self.inst.crossings:
                    W_lower = crossing.cable_lower.total_heat_per_length(
                        result.currents[step_idx, -1] if n > 1 else result.currents[step_idx, 0],
                        result.conductor_temps[step_idx, -1] if n > 1 else result.conductor_temps[step_idx, 0],
                    )
                    W_upper = crossing.cable_upper.total_heat_per_length(
                        result.currents[step_idx, 0],
                        result.conductor_temps[step_idx, 0],
                    )
                    dT_upper = crossing.transient_temperature_rise_at_upper(
                        W_lower, t if t > 0 else 1.0
                    )
                    dT_lower = crossing.transient_temperature_rise_at_lower(
                        W_upper, t if t > 0 else 1.0
                    )
                    result.conductor_temps[step_idx, 0] += dT_upper
                    if n > 1:
                        result.conductor_temps[step_idx, -1] += dT_lower

        return result
