"""Transient thermal solver using implicit Euler time integration.

Solves the ODE system  C · dθ/dt + G · θ = P(t)  arising from the
lumped-parameter thermal network, with temperature-dependent conductor
resistance updated at each time step.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from thermal_cable_model.thermal_network import CableThermalNetwork
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.ground import GroundTemperatureModel


@dataclass
class TransientResult:
    """Container for time-domain simulation results."""

    times: np.ndarray                      # [s]
    conductor_temps: np.ndarray            # (n_steps, n_cables)  [°C]
    insulation_temps: np.ndarray           # (n_steps, n_cables)
    sheath_temps: np.ndarray               # (n_steps, n_cables)
    armour_temps: np.ndarray               # (n_steps, n_cables)
    surface_temps: np.ndarray              # (n_steps, n_cables)
    soil_temps: np.ndarray                 # (n_steps, n_cables)
    currents: np.ndarray                   # (n_steps, n_cables)
    ambient_temps: np.ndarray              # (n_steps, n_cables)
    full_state: np.ndarray | None = None   # (n_steps, n_nodes)  — optional
    cable_names: list[str] = field(default_factory=list)

    @property
    def n_steps(self) -> int:
        return len(self.times)

    @property
    def n_cables(self) -> int:
        return self.conductor_temps.shape[1]

    def max_conductor_temp(self, cable_index: int = 0) -> float:
        return float(np.max(self.conductor_temps[:, cable_index]))

    def max_insulation_temp(self, cable_index: int = 0) -> float:
        return float(np.max(self.insulation_temps[:, cable_index]))

    def max_sheath_temp(self, cable_index: int = 0) -> float:
        return float(np.max(self.sheath_temps[:, cable_index]))

    def max_armour_temp(self, cable_index: int = 0) -> float:
        return float(np.max(self.armour_temps[:, cable_index]))


class TransientSolver:
    """Implicit-Euler solver for the cable thermal network.

    Parameters
    ----------
    network : CableThermalNetwork
        The assembled thermal network.
    load_profiles : list[LoadProfile]
        One per cable.
    ground_model : GroundTemperatureModel
        Provides ambient temperature at each cable's burial depth.
    """

    def __init__(
        self,
        network: CableThermalNetwork,
        load_profiles: list[LoadProfile],
        ground_model: GroundTemperatureModel,
    ):
        self.net = network
        self.loads = load_profiles
        self.ground = ground_model
        if len(load_profiles) != network.n_cables:
            raise ValueError(
                f"Need {network.n_cables} load profiles, got {len(load_profiles)}"
            )

    def _ambient_temps(self, time_s: float) -> list[float]:
        return [
            self.ground.temperature(self.net.depths[k], time_s)
            for k in range(self.net.n_cables)
        ]

    def _currents(self, time_s: float) -> list[float]:
        return [lp.current_at(time_s) for lp in self.loads]

    def solve(
        self,
        dt: float,
        duration: float,
        initial_state: np.ndarray | None = None,
        store_full_state: bool = False,
        nonlinear_iterations: int = 3,
    ) -> TransientResult:
        """Run transient simulation.

        Parameters
        ----------
        dt : float
            Time step [s].
        duration : float
            Total simulation time [s].
        initial_state : ndarray, optional
            Starting temperatures.  Defaults to steady-state at t = 0.
        store_full_state : bool
            If True, store the entire state vector at each step.
        nonlinear_iterations : int
            Number of Picard iterations per step for the resistance
            non-linearity.
        """
        net = self.net
        N = net.n_nodes
        n_steps = int(duration / dt) + 1
        times = np.linspace(0, duration, n_steps)

        nc = net.n_cables
        cond_T  = np.zeros((n_steps, nc))
        ins_T   = np.zeros((n_steps, nc))
        sh_T    = np.zeros((n_steps, nc))
        arm_T   = np.zeros((n_steps, nc))
        surf_T  = np.zeros((n_steps, nc))
        soil_T  = np.zeros((n_steps, nc))
        cur_arr = np.zeros((n_steps, nc))
        amb_arr = np.zeros((n_steps, nc))
        full = np.zeros((n_steps, N)) if store_full_state else None

        C_diag = net.C.copy()
        C_diag[C_diag < 1e-12] = 1e-12

        if initial_state is not None:
            theta = initial_state.copy()
        else:
            amb0 = self._ambient_temps(0.0)
            cur0 = self._currents(0.0)
            theta = net.steady_state(cur0, amb0)

        self._record(0, theta, times[0],
                     cond_T, ins_T, sh_T, arm_T, surf_T, soil_T,
                     cur_arr, amb_arr, full)

        A = np.diag(C_diag / dt) + net.G

        for step in range(1, n_steps):
            t = times[step]
            currents = self._currents(t)
            amb_temps = self._ambient_temps(t)

            theta_pred = theta.copy()
            for _ in range(nonlinear_iterations):
                cond_temps = net.get_conductor_temperatures(theta_pred)
                P = net.forcing_vector(currents, cond_temps, amb_temps)
                rhs = C_diag / dt * theta + P
                theta_pred = np.linalg.solve(A, rhs)

            theta = theta_pred
            self._record(step, theta, t,
                         cond_T, ins_T, sh_T, arm_T, surf_T, soil_T,
                         cur_arr, amb_arr, full)

        return TransientResult(
            times=times,
            conductor_temps=cond_T,
            insulation_temps=ins_T,
            sheath_temps=sh_T,
            armour_temps=arm_T,
            surface_temps=surf_T,
            soil_temps=soil_T,
            currents=cur_arr,
            ambient_temps=amb_arr,
            full_state=full,
            cable_names=[c.name for c in net.cables],
        )

    def _record(self, step, theta, t,
                cond_T, ins_T, sh_T, arm_T, surf_T, soil_T,
                cur_arr, amb_arr, full):
        net = self.net
        cond_T[step] = net.get_conductor_temperatures(theta)
        ins_T[step]  = net.get_insulation_temperatures(theta)
        sh_T[step]   = net.get_sheath_temperatures(theta)
        arm_T[step]  = net.get_armour_temperatures(theta)
        surf_T[step] = net.get_surface_temperatures(theta)
        soil_T[step] = net.get_soil_temperatures(theta)
        cur_arr[step] = self._currents(t)
        amb_arr[step] = self._ambient_temps(t)
        if full is not None:
            full[step] = theta
