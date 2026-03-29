2-D Finite Element Solver
=========================

In addition to the lumped-parameter thermal network, Thermal FEM includes a
2-D finite-element solver for computing the full temperature distribution in
the soil surrounding the cables.  This is useful for:

- Visualising the temperature field and isotherms
- Validating the lumped-parameter results
- Assessing the thermal influence radius of a cable installation

The FEM solver operates on a rectangular domain with bilinear quadrilateral
(Q4) elements on a graded Cartesian mesh.

Setting up the FEM solver
--------------------------

.. code-block:: python

   from thermal_fem.cable import Cable
   from thermal_fem.ground import KasudaModel
   from thermal_fem.materials import SOIL_STANDARD
   from thermal_fem.fem.solver import FEMSolver

   cable = Cable.single_core_xlpe_cu(240, voltage_class="MV")
   ground = KasudaModel(mean_surface_temp=10.0, annual_amplitude=12.0)

   fem = FEMSolver(
       domain_x=(-1.5, 1.5),   # horizontal extent [m]
       domain_y=(0.0, 3.0),    # depth extent [m]
       soil=SOIL_STANDARD,
       ground_model=ground,
       base_nx=50,              # baseline horizontal divisions
       base_ny=50,              # baseline vertical divisions
   )

Adding cable heat sources
--------------------------

Register cables as volumetric heat sources.  Provide the total heat per
unit length (compute this from the cable and load conditions):

.. code-block:: python

   I = 400.0                # current [A]
   T_est = 70.0             # estimated conductor temperature [°C]
   W = cable.total_heat_per_length(I, T_est)

   fem.add_cable(cable, x=-0.2, depth=1.2, heat_rate=W)
   fem.add_cable(cable, x=0.2,  depth=1.2, heat_rate=W)

Steady-state solution
---------------------

.. code-block:: python

   t_summer = 200 * 86400  # mid-summer for boundary temperatures
   result = fem.solve_steady_state(time_s=t_summer)

   field = result.field_at(0)
   print(f"Temperature range: {field.min():.1f} – {field.max():.1f} °C")

Transient solution
------------------

.. code-block:: python

   result = fem.solve_transient(dt=3600, duration=24 * 3600)

For time-dependent heat sources, provide a callable:

.. code-block:: python

   def heat_func(time_s):
       """Return a list of heat rates [W/m] per cable."""
       return [W * (1 if time_s < 12*3600 else 0.5)] * 2

   result = fem.solve_transient(dt=3600, duration=24*3600, heat_rate_func=heat_func)

Plotting the temperature field
-------------------------------

.. code-block:: python

   from thermal_fem import plot_soil_temperature_field

   fig = plot_soil_temperature_field(
       result.field_at(0), result.x, result.y,
       cable_positions=[(-0.2, 1.2), (0.2, 1.2)],
   )
   fig.savefig("temperature_field.png", dpi=150)

Boundary conditions
-------------------

The FEM solver applies:

- **Dirichlet (prescribed temperature)** on the top and bottom boundaries,
  evaluated from the ground temperature model at the appropriate depth and
  time.
- **Neumann (zero heat flux)** on the left and right boundaries, representing
  undisturbed far-field conditions.

Mesh refinement
---------------

The mesh is automatically refined near cable locations.  The
:class:`~thermal_fem.fem.mesh.RectangularMesh` uses a grading function that
transitions from fine elements (near cables) to coarse elements (far field).

Parameters controlling the mesh density:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Description
   * - ``base_nx``, ``base_ny``
     - 60
     - Baseline number of divisions in each direction
   * - ``refinement_radius``
     - 0.3 m
     - Radius of the fine-mesh zone around each cable
   * - ``refinement_factor``
     - 3.0
     - Ratio of coarse to fine element size
