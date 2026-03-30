"""IEC 60287 / IEC 60853 thermal resistance and capacitance network.

Builds a lumped-parameter thermal circuit for each cable, including the
cable-internal resistances T1–T3 and the external (soil) resistance T4.
Mutual heating between parallel cables is computed via the image method.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from thermal_cable_model.cable import Cable
from thermal_cable_model.materials import ThermalMaterial


# ═══════════════════════════════════════════════════════════════════════
#  Cable-internal thermal resistances (T1, T2, T3)
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class InternalThermalResistances:
    """Per-unit-length thermal resistances of the cable layers [K·m/W]."""
    T1: float  # insulation
    T2: float  # bedding / metallic screen-to-armour
    T3: float  # outer serving / jacket

    @classmethod
    def from_cable(cls, cable: Cable) -> InternalThermalResistances:
        """Compute T1–T3 by summing layer thermal resistances.

        The first insulation layer is T1, intermediate layers form T2,
        and the outermost jacket is T3.  This mapping follows IEC 60287-2-1
        conventions for typical cable constructions.
        """
        if not cable.layers:
            return cls(T1=0.0, T2=0.0, T3=0.0)

        resistances = [l.thermal_resistance_per_length for l in cable.layers]

        if len(resistances) == 1:
            return cls(T1=resistances[0], T2=0.0, T3=0.0)
        if len(resistances) == 2:
            return cls(T1=resistances[0], T2=0.0, T3=resistances[1])

        return cls(
            T1=resistances[0],
            T2=sum(resistances[1:-1]),
            T3=resistances[-1],
        )


@dataclass
class InternalThermalCapacitances:
    """Per-unit-length thermal capacitances [J/(m·K)]."""
    Q_conductor: float
    Q_insulation: float  # includes screen
    Q_bedding: float
    Q_armour: float
    Q_jacket: float

    @classmethod
    def from_cable(cls, cable: Cable) -> InternalThermalCapacitances:
        caps = [l.thermal_capacitance_per_length for l in cable.layers]
        q_c = cable.conductor_capacitance_per_length
        return cls(
            Q_conductor=q_c,
            Q_insulation=caps[0] if len(caps) > 0 else 0.0,
            Q_bedding=sum(caps[1:-1]) if len(caps) > 2 else 0.0,
            Q_armour=0.0,
            Q_jacket=caps[-1] if len(caps) > 1 else 0.0,
        )


# ═══════════════════════════════════════════════════════════════════════
#  External (soil) thermal resistance T4
# ═══════════════════════════════════════════════════════════════════════

def external_thermal_resistance(
    depth: float,
    cable_outer_radius: float,
    soil: ThermalMaterial,
) -> float:
    """T4 for a single isolated cable buried at *depth* [m].

    IEC 60287-2-1 §2.2.7:
        T4 = (ρ_soil / 2π) · ln(2L / D_e)

    where L = burial depth to cable centre, D_e = external diameter.
    """
    D_e = 2.0 * cable_outer_radius
    u = 2.0 * depth / D_e
    if u <= 1.0:
        raise ValueError(
            f"Cable outer diameter ({D_e*1e3:.1f} mm) exceeds twice the "
            f"burial depth ({depth*1e3:.1f} mm) — geometry invalid."
        )
    return soil.thermal_resistivity / (2.0 * math.pi) * math.log(
        u + math.sqrt(u ** 2 - 1.0)
    )


def mutual_heating_resistance(
    xi: float, yi: float,
    xj: float, yj: float,
    soil: ThermalMaterial,
) -> float:
    """Temperature rise at cable *i* due to unit heat from cable *j*.

    Uses the image method for a semi-infinite conducting half-space with
    isothermal surface (ground level at y = 0, cables at y < 0 i.e.
    depth is positive downward so yi, yj > 0).

    Returns ΔT4_ij = ρ_soil/(2π) · ln(d'_ij / d_ij)

    d_ij   = real distance between cables i and j
    d'_ij  = distance from cable i to the *image* of cable j
    """
    dx = xi - xj
    dy_real = yi - yj
    dy_image = yi + yj  # image is at −yj (mirrored above ground)

    d_real = math.sqrt(dx ** 2 + dy_real ** 2)
    d_image = math.sqrt(dx ** 2 + dy_image ** 2)

    if d_real < 1e-9:
        return 0.0  # same cable

    return soil.thermal_resistivity / (2.0 * math.pi) * math.log(
        d_image / d_real
    )


# ═══════════════════════════════════════════════════════════════════════
#  Complete thermal network for a cable group
# ═══════════════════════════════════════════════════════════════════════

class CableThermalNetwork:
    """State-space thermal circuit for one or more parallel cables.

    For *N* cables the state vector has *N × M* entries where *M* is the
    number of thermal nodes per cable (conductor, insulation-midpoint,
    screen/sheath, jacket-surface).  The network is represented as

        C · dθ/dt + G · θ = P(t)

    where C is the (diagonal) capacitance matrix, G the conductance matrix,
    and P the forcing vector (heat sources + boundary coupling).

    Nodes per cable (4-node model)
    ──────────────────────────────
    0 : conductor            θ_c
    1 : insulation mid-point θ_i
    2 : cable surface        θ_s  (outer jacket)
    3 : soil node            θ_soil (near cable)

    Coupling to ambient (ground temperature) is through T4 from node 3.
    """

    NODES_PER_CABLE = 4

    def __init__(
        self,
        cables: list[Cable],
        positions_x: list[float],
        depths: list[float],
        soil: ThermalMaterial,
    ):
        self.cables = cables
        self.positions_x = list(positions_x)
        self.depths = list(depths)
        self.soil = soil
        self.n_cables = len(cables)
        self.n_nodes = self.n_cables * self.NODES_PER_CABLE
        self._build_network()

    def _build_network(self) -> None:
        N = self.n_nodes
        self.C = np.zeros(N)       # diagonal capacitance
        self.G = np.zeros((N, N))  # conductance matrix
        # which node indices are "conductor" / "surface" / "soil"
        self._conductor_idx = []
        self._insulation_idx = []
        self._surface_idx = []
        self._soil_idx = []

        for k, cable in enumerate(self.cables):
            i0 = k * self.NODES_PER_CABLE
            ic, ii, is_, ig = i0, i0 + 1, i0 + 2, i0 + 3

            self._conductor_idx.append(ic)
            self._insulation_idx.append(ii)
            self._surface_idx.append(is_)
            self._soil_idx.append(ig)

            tr = InternalThermalResistances.from_cable(cable)
            tc = InternalThermalCapacitances.from_cable(cable)

            # Van-Wormer coefficients for splitting resistances
            p1 = _van_wormer(cable.layers[0]) if cable.layers else 0.5
            T1 = tr.T1 / cable.n_conductors if cable.n_conductors > 0 else tr.T1

            T_inner = p1 * T1          # conductor → insulation mid
            T_outer = (1.0 - p1) * T1  # insulation mid → screen
            T23 = (tr.T2 + tr.T3) / cable.n_conductors if cable.n_conductors > 0 else (tr.T2 + tr.T3)
            T4 = external_thermal_resistance(
                self.depths[k], cable.outer_radius, self.soil
            )

            # Conductances (G_ij = 1 / R_ij)
            g_c_i = 1.0 / max(T_inner, 1e-12)
            g_i_s = 1.0 / max(T_outer + T23, 1e-12)
            g_s_g = 1.0 / max(T4 * 0.5, 1e-12)

            self.G[ic, ic] += g_c_i
            self.G[ic, ii] -= g_c_i
            self.G[ii, ic] -= g_c_i
            self.G[ii, ii] += g_c_i + g_i_s
            self.G[ii, is_] -= g_i_s
            self.G[is_, ii] -= g_i_s
            self.G[is_, is_] += g_i_s + g_s_g
            self.G[is_, ig] -= g_s_g
            self.G[ig, is_] -= g_s_g
            self.G[ig, ig] += g_s_g

            # Capacitances
            self.C[ic] = tc.Q_conductor * cable.n_conductors
            self.C[ii] = tc.Q_insulation * cable.n_conductors
            self.C[is_] = tc.Q_jacket + tc.Q_bedding
            # Soil node capacitance — effective soil annulus
            r_out_soil = min(self.depths[k], 0.5)  # up to 0.5 m annulus
            r_in_soil = cable.outer_radius
            self.C[ig] = (
                self.soil.volumetric_heat_capacity
                * math.pi * (r_out_soil ** 2 - r_in_soil ** 2)
            )

        # Mutual heating resistance matrix (IEC 60287 image method).
        # Rm[i][j] = temperature rise at cable i per unit heat from cable j.
        # Used in the forcing vector to raise the effective ambient
        # temperature seen by each cable's soil node.
        self._Rm = np.zeros((self.n_cables, self.n_cables))
        for i in range(self.n_cables):
            for j in range(self.n_cables):
                if i == j:
                    continue
                self._Rm[i, j] = mutual_heating_resistance(
                    self.positions_x[i], self.depths[i],
                    self.positions_x[j], self.depths[j],
                    self.soil,
                )

        # External coupling: soil-node to ambient (other half of T4)
        self._g_ambient = np.zeros(N)
        for k, cable in enumerate(self.cables):
            T4 = external_thermal_resistance(
                self.depths[k], cable.outer_radius, self.soil
            )
            ig = self._soil_idx[k]
            g_ext = 1.0 / max(T4 * 0.5, 1e-12)
            self.G[ig, ig] += g_ext
            self._g_ambient[ig] = g_ext

    def forcing_vector(
        self,
        currents: list[float],
        conductor_temps: list[float],
        ambient_temps: list[float],
    ) -> np.ndarray:
        """Build the right-hand-side vector P.

        Parameters
        ----------
        currents : per-cable RMS current [A]
        conductor_temps : current conductor temperatures for R(T) [°C]
        ambient_temps : per-cable ambient ground temperature [°C]
        """
        P = np.zeros(self.n_nodes)

        # Total heat per cable (needed for mutual heating calculation)
        W_total = np.array([
            cable.total_heat_per_length(currents[k], conductor_temps[k])
            for k, cable in enumerate(self.cables)
        ])

        for k, cable in enumerate(self.cables):
            ic = self._conductor_idx[k]
            Wc = cable.n_conductors * cable.conductor_loss(
                currents[k], conductor_temps[k]
            )
            Wd = cable.n_conductors * cable.dielectric_loss
            P[ic] += Wc
            ii = self._insulation_idx[k]
            P[ii] += Wd * 0.5
            P[ic] += Wd * 0.5

            # Sheath + armour losses attributed to surface node
            Ws = Wc * cable.loss_factor_sheath * cable.n_conductors
            Wa = Wc * cable.loss_factor_armour * cable.n_conductors
            is_ = self._surface_idx[k]
            P[is_] += Ws + Wa

            # Mutual heating: raise effective ambient by ΔT_mutual_k
            # per IEC 60287: ΔT_k = Σ_{j≠k} W_j × Rm[k,j]
            delta_T_mutual = float(self._Rm[k, :] @ W_total)
            ig = self._soil_idx[k]
            P[ig] += self._g_ambient[ig] * (ambient_temps[k] + delta_T_mutual)

        return P

    def steady_state(
        self,
        currents: list[float],
        ambient_temps: list[float],
        tol: float = 0.01,
        max_iter: int = 50,
    ) -> np.ndarray:
        """Iterative steady-state solve (resistance is temperature-dependent).

        Returns the full state vector θ [°C].
        """
        theta = np.full(self.n_nodes, np.mean(ambient_temps))

        for _ in range(max_iter):
            cond_temps = [theta[ic] for ic in self._conductor_idx]
            P = self.forcing_vector(currents, cond_temps, ambient_temps)
            theta_new = np.linalg.solve(self.G, P)
            if np.max(np.abs(theta_new - theta)) < tol:
                return theta_new
            theta = theta_new

        return theta

    def get_conductor_temperatures(self, theta: np.ndarray) -> list[float]:
        return [float(theta[i]) for i in self._conductor_idx]

    def get_insulation_temperatures(self, theta: np.ndarray) -> list[float]:
        return [float(theta[i]) for i in self._insulation_idx]

    def get_surface_temperatures(self, theta: np.ndarray) -> list[float]:
        return [float(theta[i]) for i in self._surface_idx]

    def get_soil_temperatures(self, theta: np.ndarray) -> list[float]:
        return [float(theta[i]) for i in self._soil_idx]


def _van_wormer(layer: "CableLayer") -> float:
    """Van Wormer coefficient for splitting a cylindrical thermal resistance.

    p = (1 / (2·ln(r2/r1))) - 1 / ((r2/r1)² - 1)
    Returns value in [0, 0.5]; defaults to 0.5 for thin layers.
    """
    ratio = layer.outer_radius / max(layer.inner_radius, 1e-9)
    if ratio < 1.001:
        return 0.5
    ln_r = math.log(ratio)
    return 1.0 / (2.0 * ln_r) - 1.0 / (ratio ** 2 - 1.0)
