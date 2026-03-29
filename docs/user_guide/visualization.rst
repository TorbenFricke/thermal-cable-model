Visualisation
=============

Thermal FEM provides three plotting functions that produce publication-ready
Matplotlib figures.

Temperature history
-------------------

:func:`~thermal_fem.visualization.plot_temperature_history` plots conductor,
insulation, surface, and soil node temperatures over time, with an optional
current panel:

.. code-block:: python

   from thermal_fem import plot_temperature_history

   fig = plot_temperature_history(
       result,
       cable_indices=[0, 1],  # which cables to show (default: all)
       show_ambient=True,     # overlay ambient ground temperature
       show_current=True,     # add current panel below
       time_unit="hours",     # "seconds", "hours", or "days"
       figsize=(14, 8),
   )
   fig.savefig("temperatures.png", dpi=150)

Cross-section
-------------

:func:`~thermal_fem.visualization.plot_cross_section` draws a 2-D view of
the cable arrangement in the soil with optional temperature annotations:

.. code-block:: python

   from thermal_fem.visualization import plot_cross_section

   fig = plot_cross_section(
       positions_x=[0.0],
       depths=[1.2],
       outer_radii=[cable.outer_radius],
       temperatures=[result.max_conductor_temp(0)],
       cable_names=[cable.name],
       soil_extent=(-1.0, 1.0, 2.5),
   )
   fig.savefig("cross_section.png", dpi=150)

Soil temperature field
----------------------

:func:`~thermal_fem.visualization.plot_soil_temperature_field` creates a
filled contour plot from a 2-D FEM result:

.. code-block:: python

   from thermal_fem import plot_soil_temperature_field

   fig = plot_soil_temperature_field(
       fem_result.field_at(0),
       fem_result.x,
       fem_result.y,
       cable_positions=[(-0.2, 1.2), (0.2, 1.2)],
       cmap="hot_r",
   )
   fig.savefig("field.png", dpi=150)

All functions return a :class:`matplotlib.figure.Figure` object that can be
further customised before saving.
