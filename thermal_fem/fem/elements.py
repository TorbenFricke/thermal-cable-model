"""Bilinear quadrilateral (Q4) element matrices and global assembly.

Provides the element-level stiffness and mass matrices for the 2-D
steady-state and transient heat equation, plus the global assembly routine
using sparse matrices.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

from thermal_fem.fem.mesh import RectangularMesh


# ── Gauss quadrature (2×2) in natural coordinates ───────────────────
_GP = 1.0 / np.sqrt(3.0)
_GAUSS_PTS = np.array([[-_GP, -_GP], [_GP, -_GP], [_GP, _GP], [-_GP, _GP]])
_GAUSS_W = np.ones(4)


def _shape_functions(xi: float, eta: float) -> np.ndarray:
    """Evaluate the four Q4 shape functions at (ξ, η)."""
    return 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ])


def _shape_derivatives(xi: float, eta: float) -> np.ndarray:
    """dN/dξ and dN/dη, shape (2, 4)."""
    return 0.25 * np.array([
        [-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],
        [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)],
    ])


def _element_matrices(
    coords: np.ndarray,
    conductivity: float,
    vol_heat_cap: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute K_e (4×4) and M_e (4×4) for a single Q4 element.

    Parameters
    ----------
    coords : (4, 2) array of nodal (x, y) coordinates.
    conductivity : thermal conductivity λ [W/(m·K)].
    vol_heat_cap : volumetric heat capacity ρ·c_p [J/(m³·K)].
    """
    Ke = np.zeros((4, 4))
    Me = np.zeros((4, 4))

    for gp, w in zip(_GAUSS_PTS, _GAUSS_W):
        xi, eta = gp
        dN_dxi = _shape_derivatives(xi, eta)
        J = dN_dxi @ coords  # (2, 2) Jacobian
        detJ = np.linalg.det(J)
        dN_dx = np.linalg.solve(J, dN_dxi)  # (2, 4)

        N = _shape_functions(xi, eta)

        Ke += conductivity * (dN_dx.T @ dN_dx) * detJ * w
        Me += vol_heat_cap * np.outer(N, N) * detJ * w

    return Ke, Me


def assemble_system(
    mesh: RectangularMesh,
    conductivity: float | np.ndarray,
    vol_heat_cap: float | np.ndarray,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Assemble global stiffness K and mass M in sparse CSR format.

    Parameters
    ----------
    mesh : RectangularMesh
    conductivity : float or (n_elements,) array
        Element-wise thermal conductivity.
    vol_heat_cap : float or (n_elements,) array
        Element-wise volumetric heat capacity.

    Returns
    -------
    K, M : sparse CSR matrices, shape (n_nodes, n_nodes).
    """
    n = mesh.n_nodes
    elems = mesh.elements
    coords_all = mesh.node_coords

    scalar_k = np.isscalar(conductivity)
    scalar_c = np.isscalar(vol_heat_cap)

    rows, cols, k_vals, m_vals = [], [], [], []

    for e_idx in range(mesh.n_elements):
        nodes = elems[e_idx]
        e_coords = coords_all[nodes]
        k = conductivity if scalar_k else conductivity[e_idx]
        c = vol_heat_cap if scalar_c else vol_heat_cap[e_idx]
        Ke, Me = _element_matrices(e_coords, k, c)

        for i_loc in range(4):
            for j_loc in range(4):
                rows.append(nodes[i_loc])
                cols.append(nodes[j_loc])
                k_vals.append(Ke[i_loc, j_loc])
                m_vals.append(Me[i_loc, j_loc])

    K = sparse.csr_matrix((k_vals, (rows, cols)), shape=(n, n))
    M = sparse.csr_matrix((m_vals, (rows, cols)), shape=(n, n))
    return K, M
