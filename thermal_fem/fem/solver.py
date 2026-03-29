"""2-D finite-element solver for the soil temperature field.

Solves the transient heat equation

    ρ·c_p · ∂T/∂t = λ · ∇²T + Q

on a rectangular domain with:
  - prescribed temperature (Dirichlet) on the top surface and bottom boundary,
  - zero-flux (Neumann) on the left and right sides,
  - volumetric heat sources at cable locations.

Uses implicit Euler time integration with a direct sparse solver.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from thermal_fem.cable import Cable
from thermal_fem.fem.elements import assemble_system
from thermal_fem.fem.mesh import RectangularMesh
from thermal_fem.ground import GroundTemperatureModel
from thermal_fem.materials import ThermalMaterial


@dataclass
class FEMResult:
    """Container for the 2-D temperature field results."""

    x: np.ndarray            # (nx,)
    y: np.ndarray            # (ny,)
    times: np.ndarray        # (n_steps,)
    fields: np.ndarray       # (n_steps, ny, nx)  temperature [°C]

    def field_at(self, step: int) -> np.ndarray:
        return self.fields[step]


class FEMSolver:
    """2-D soil thermal FEM solver.

    Parameters
    ----------
    domain_x : (float, float)
        Horizontal extent [m].
    domain_y : (float, float)
        Depth extent [m] — (0, max_depth).
    soil : ThermalMaterial
    ground_model : GroundTemperatureModel
    base_nx, base_ny : int
        Baseline mesh divisions.
    """

    def __init__(
        self,
        domain_x: tuple[float, float] = (-2.0, 2.0),
        domain_y: tuple[float, float] = (0.0, 4.0),
        soil: ThermalMaterial | None = None,
        ground_model: GroundTemperatureModel | None = None,
        base_nx: int = 60,
        base_ny: int = 60,
    ):
        from thermal_fem.materials import SOIL_STANDARD
        from thermal_fem.ground import ConstantGroundTemperature

        self.domain_x = domain_x
        self.domain_y = domain_y
        self.soil = soil or SOIL_STANDARD
        self.ground = ground_model or ConstantGroundTemperature(15.0)
        self.base_nx = base_nx
        self.base_ny = base_ny

        self._cables: list[Cable] = []
        self._cable_positions: list[tuple[float, float]] = []
        self._cable_heat_rates: list[float] = []

    def add_cable(
        self,
        cable: Cable,
        x: float,
        depth: float,
        heat_rate: float,
    ) -> None:
        """Register a cable as a heat source in the domain.

        Parameters
        ----------
        cable : Cable
        x, depth : float
            Position [m].
        heat_rate : float
            Total heat per unit length [W/m] to inject.
        """
        self._cables.append(cable)
        self._cable_positions.append((x, depth))
        self._cable_heat_rates.append(heat_rate)

    def solve_steady_state(self, time_s: float = 0.0) -> FEMResult:
        """Solve the steady-state temperature field.

        Parameters
        ----------
        time_s : float
            Time [s] for evaluating the ground temperature boundary.
        """
        mesh = RectangularMesh.create(
            self.domain_x, self.domain_y,
            self.base_nx, self.base_ny,
            cable_positions=self._cable_positions,
        )
        K, _ = assemble_system(
            mesh, self.soil.thermal_conductivity,
            self.soil.volumetric_heat_capacity,
        )

        n = mesh.n_nodes
        rhs = np.zeros(n)

        # Heat sources — distribute cable heat to nearest nodes
        coords = mesh.node_coords
        for (cx, cy), W in zip(self._cable_positions, self._cable_heat_rates):
            dist = np.sqrt((coords[:, 0] - cx) ** 2 + (coords[:, 1] - cy) ** 2)
            near = np.where(dist < 0.05)[0]
            if len(near) == 0:
                near = np.array([np.argmin(dist)])
            rhs[near] += W / len(near)

        # Boundary conditions
        K_bc, rhs_bc, free = self._apply_bc(mesh, K, rhs, time_s)
        T = np.zeros(n)
        T_top = self.ground.temperature(0.0, time_s)
        T_bot = self.ground.temperature(self.domain_y[1], time_s)

        # Set Dirichlet values
        for ix in range(mesh.nx):
            T[mesh.node_index(ix, 0)] = T_top
            T[mesh.node_index(ix, mesh.ny - 1)] = T_bot

        T[free] = spsolve(K_bc, rhs_bc)

        field = T.reshape(mesh.ny, mesh.nx)
        return FEMResult(
            x=mesh.x, y=mesh.y,
            times=np.array([time_s]),
            fields=field[np.newaxis, :, :],
        )

    def solve_transient(
        self,
        dt: float,
        duration: float,
        heat_rate_func: callable | None = None,
        store_every: int = 1,
    ) -> FEMResult:
        """Solve the transient temperature field.

        Parameters
        ----------
        dt : float
            Time step [s].
        duration : float
            Simulation duration [s].
        heat_rate_func : callable(time_s) → list[float], optional
            Time-dependent heat rates [W/m] per cable.  Defaults to the
            constant values set via ``add_cable``.
        store_every : int
            Store results every N steps (to limit memory).
        """
        mesh = RectangularMesh.create(
            self.domain_x, self.domain_y,
            self.base_nx, self.base_ny,
            cable_positions=self._cable_positions,
        )
        K, M = assemble_system(
            mesh, self.soil.thermal_conductivity,
            self.soil.volumetric_heat_capacity,
        )
        n = mesh.n_nodes
        coords = mesh.node_coords
        n_steps = int(duration / dt) + 1

        # Pre-compute cable node mapping
        cable_nodes = []
        for cx, cy in self._cable_positions:
            dist = np.sqrt((coords[:, 0] - cx) ** 2 + (coords[:, 1] - cy) ** 2)
            near = np.where(dist < 0.05)[0]
            if len(near) == 0:
                near = np.array([np.argmin(dist)])
            cable_nodes.append(near)

        # Initial condition
        T = np.full(n, self.ground.temperature(self.domain_y[1] / 2, 0.0))
        for ix in range(mesh.nx):
            T[mesh.node_index(ix, 0)] = self.ground.temperature(0.0, 0.0)
            T[mesh.node_index(ix, mesh.ny - 1)] = self.ground.temperature(
                self.domain_y[1], 0.0
            )

        # Storage
        store_times = []
        store_fields = []

        A_system = M / dt + K

        for step in range(n_steps):
            t = step * dt

            if step % store_every == 0:
                store_times.append(t)
                store_fields.append(T.reshape(mesh.ny, mesh.nx).copy())

            # Build RHS
            rhs = (M / dt).dot(T)

            # Add cable heat sources
            if heat_rate_func is not None:
                heat_rates = heat_rate_func(t)
            else:
                heat_rates = self._cable_heat_rates
            for nodes, W in zip(cable_nodes, heat_rates):
                rhs[nodes] += W / len(nodes)

            # Apply BCs and solve
            A_bc, rhs_bc, free = self._apply_bc(
                mesh, A_system, rhs, t + dt
            )
            T_new = T.copy()
            T_top = self.ground.temperature(0.0, t + dt)
            T_bot = self.ground.temperature(self.domain_y[1], t + dt)
            for ix in range(mesh.nx):
                T_new[mesh.node_index(ix, 0)] = T_top
                T_new[mesh.node_index(ix, mesh.ny - 1)] = T_bot

            T_new[free] = spsolve(A_bc, rhs_bc)
            T = T_new

        return FEMResult(
            x=mesh.x, y=mesh.y,
            times=np.array(store_times),
            fields=np.array(store_fields),
        )

    def _apply_bc(
        self,
        mesh: RectangularMesh,
        A: sparse.csr_matrix,
        rhs: np.ndarray,
        time_s: float,
    ) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
        """Apply Dirichlet BCs on top and bottom rows.

        Returns the reduced system matrix, modified RHS, and the free DOF indices.
        """
        n = mesh.n_nodes
        bc_nodes = set()
        bc_values = {}

        T_top = self.ground.temperature(0.0, time_s)
        T_bot = self.ground.temperature(self.domain_y[1], time_s)

        for ix in range(mesh.nx):
            idx_top = mesh.node_index(ix, 0)
            idx_bot = mesh.node_index(ix, mesh.ny - 1)
            bc_nodes.add(idx_top)
            bc_nodes.add(idx_bot)
            bc_values[idx_top] = T_top
            bc_values[idx_bot] = T_bot

        free = np.array(sorted(set(range(n)) - bc_nodes))
        bc_arr = np.array(sorted(bc_nodes))

        A_dense = A if isinstance(A, np.ndarray) else A.toarray()
        rhs_mod = rhs.copy()

        for bc_idx in bc_arr:
            rhs_mod -= A_dense[:, bc_idx] * bc_values[bc_idx]

        A_free = sparse.csr_matrix(A_dense[np.ix_(free, free)])
        rhs_free = rhs_mod[free]
        return A_free, rhs_free, free
