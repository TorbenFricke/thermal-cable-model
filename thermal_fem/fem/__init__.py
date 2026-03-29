"""2-D finite-element sub-package for soil temperature field computation."""

from thermal_fem.fem.mesh import RectangularMesh
from thermal_fem.fem.elements import assemble_system
from thermal_fem.fem.solver import FEMSolver
