"""Structured rectangular mesh for 2-D soil temperature field analysis.

The mesh uses bilinear quadrilateral elements on a graded Cartesian grid
with finer spacing near the cable locations.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class RectangularMesh:
    """Axis-aligned rectangular FEM mesh.

    Attributes
    ----------
    x : ndarray, shape (nx,)
        Horizontal node coordinates [m].
    y : ndarray, shape (ny,)
        Vertical (depth) coordinates [m], increasing downward.
    nx, ny : int
        Number of nodes in each direction.
    node_coords : ndarray, shape (n_nodes, 2)
    elements : ndarray, shape (n_elements, 4)
        Node indices (counter-clockwise) for each Q4 element.
    """

    x: np.ndarray
    y: np.ndarray

    @property
    def nx(self) -> int:
        return len(self.x)

    @property
    def ny(self) -> int:
        return len(self.y)

    @property
    def n_nodes(self) -> int:
        return self.nx * self.ny

    @property
    def n_elements(self) -> int:
        return (self.nx - 1) * (self.ny - 1)

    def node_index(self, ix: int, iy: int) -> int:
        return iy * self.nx + ix

    @property
    def node_coords(self) -> np.ndarray:
        """(n_nodes, 2) array of (x, y) for every node."""
        X, Y = np.meshgrid(self.x, self.y)
        return np.column_stack([X.ravel(), Y.ravel()])

    @property
    def elements(self) -> np.ndarray:
        """(n_elements, 4) connectivity — CCW node indices per Q4 element."""
        elems = []
        for iy in range(self.ny - 1):
            for ix in range(self.nx - 1):
                n0 = self.node_index(ix, iy)
                n1 = self.node_index(ix + 1, iy)
                n2 = self.node_index(ix + 1, iy + 1)
                n3 = self.node_index(ix, iy + 1)
                elems.append([n0, n1, n2, n3])
        return np.array(elems, dtype=int)

    @classmethod
    def create(
        cls,
        x_range: tuple[float, float],
        y_range: tuple[float, float],
        base_nx: int = 60,
        base_ny: int = 60,
        cable_positions: list[tuple[float, float]] | None = None,
        refinement_radius: float = 0.3,
        refinement_factor: float = 3.0,
    ) -> RectangularMesh:
        """Build a graded mesh with optional refinement near cables.

        Parameters
        ----------
        x_range : (x_min, x_max)
        y_range : (y_min, y_max)   y_min is ground surface (≥ 0), y_max is
            maximum depth.
        base_nx, base_ny : int
            Baseline number of divisions.
        cable_positions : list of (x, depth)
            Locations around which the mesh is refined.
        refinement_radius : float
            Radius [m] of the refinement zone around each cable.
        refinement_factor : float
            Ratio of coarse to fine element size.
        """
        x = _graded_coords(
            x_range[0], x_range[1], base_nx,
            [p[0] for p in cable_positions] if cable_positions else [],
            refinement_radius, refinement_factor,
        )
        y = _graded_coords(
            y_range[0], y_range[1], base_ny,
            [p[1] for p in cable_positions] if cable_positions else [],
            refinement_radius, refinement_factor,
        )
        return cls(x=x, y=y)


def _graded_coords(
    lo: float,
    hi: float,
    n_base: int,
    focal_points: list[float],
    r_refine: float,
    factor: float,
) -> np.ndarray:
    """Create a 1-D graded coordinate array.

    The spacing is uniform at *h_fine* near each focal point and gradually
    coarsens to *h_coarse = factor × h_fine* elsewhere.
    """
    h_coarse = (hi - lo) / n_base
    h_fine = h_coarse / factor

    coords = set()
    coords.add(lo)
    coords.add(hi)

    for fp in focal_points:
        zone_lo = max(lo, fp - r_refine)
        zone_hi = min(hi, fp + r_refine)
        n_fine = max(int((zone_hi - zone_lo) / h_fine), 4)
        coords.update(np.linspace(zone_lo, zone_hi, n_fine).tolist())

    # Fill the rest with coarse spacing
    pos = lo
    while pos < hi - 1e-12:
        coords.add(pos)
        # Find distance to nearest focal point
        min_dist = min(
            (abs(pos - fp) for fp in focal_points), default=r_refine + 1.0
        )
        if min_dist < r_refine:
            pos += h_fine
        else:
            blend = min((min_dist - r_refine) / r_refine, 1.0)
            h = h_fine + blend * (h_coarse - h_fine)
            pos += h
    coords.add(hi)

    return np.array(sorted(coords))
