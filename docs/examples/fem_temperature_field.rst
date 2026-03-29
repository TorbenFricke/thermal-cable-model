FEM Temperature Field
=====================

**Script:** ``examples/fem_temperature_field.py``

Computes and plots the steady-state 2-D soil temperature field around two
parallel MV cables using the finite element method.

Setup
-----

- **Cable:** 1 × 240 mm\ :sup:`2` Cu XLPE, 20 kV
- **Load:** 400 A (heat rate computed at estimated 70 °C conductor)
- **Domain:** −1.5 to 1.5 m horizontal, 0 to 3 m depth
- **Cable positions:** x = ±0.2 m at 1.2 m depth
- **Ground temperature:** Kasuda model, evaluated at mid-summer (day 200)
- **Mesh:** 50 × 50 baseline divisions with automatic refinement near cables

What it demonstrates
--------------------

1. Setting up the :class:`~thermal_fem.fem.solver.FEMSolver`
2. Adding cable heat sources
3. Solving the steady-state temperature field
4. Plotting with :func:`~thermal_fem.visualization.plot_soil_temperature_field`

Output files
------------

- ``fem_temperature_field.png`` — contour plot of the soil temperature field
  with cable positions marked
