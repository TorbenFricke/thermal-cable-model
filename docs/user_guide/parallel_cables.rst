Parallel Cable Installations
============================

When multiple cables are buried near each other, each cable raises the soil
temperature in its vicinity, increasing the operating temperature of its
neighbours.  Thermal Cable Model accounts for this *mutual heating* automatically
when multiple cables are added to a :class:`~thermal_cable_model.CableInstallation`.

Setting up parallel cables
--------------------------

Add cables at different horizontal positions (``x``) and burial depths:

.. code-block:: python

   from thermal_cable_model import Cable, CableInstallation, ThermalSimulation, LoadProfile
   from thermal_cable_model.materials import SOIL_STANDARD

   cable = Cable.single_core_xlpe_cu(240, voltage_class="MV")

   inst = CableInstallation(soil=SOIL_STANDARD)
   inst.add_cable(cable, x=-0.3, depth=1.0, load=LoadProfile.constant(450, 48*3600))
   inst.add_cable(cable, x=0.0,  depth=1.0, load=LoadProfile.constant(400, 48*3600))
   inst.add_cable(cable, x=0.3,  depth=1.0, load=LoadProfile.constant(350, 48*3600))

   sim = ThermalSimulation(inst)
   result = sim.run_transient(dt=600, duration=48 * 3600)

   for i in range(3):
       print(f"Cable {i}: max conductor = {result.max_conductor_temp(i):.1f} °C")

Image method for mutual heating
--------------------------------

The mutual thermal resistance between cables *i* and *j* is computed using
the image method for a semi-infinite conducting half-space with an
isothermal ground surface:

.. math::

   \Delta T_{4,ij} = \frac{\rho_\text{soil}}{2\pi}
   \ln\!\left(\frac{d'_{ij}}{d_{ij}}\right)

where:

- :math:`d_{ij}` is the real distance between cables *i* and *j*
- :math:`d'_{ij}` is the distance from cable *i* to the **image** of cable *j*
  (reflected about the ground surface)

The temperature rise at cable *i* due to heat from all other cables is:

.. math::

   \Delta T_{\text{mutual},i} = \sum_{j \ne i} W_j \cdot \Delta T_{4,ij}

This coupling is included in the forcing vector of the thermal network at
every time step, so the mutual heating evolves dynamically with the changing
cable currents.

Different load profiles
-----------------------

Each cable can carry a different load profile — this is useful for modelling
cables from different circuits (industrial, residential, constant base load)
sharing the same trench:

.. code-block:: python

   load_industrial = LoadProfile.cyclic(500, 100, 8*3600, 0.6, 6)
   load_base       = LoadProfile.constant(450, 48*3600)
   load_residential = LoadProfile.daily_pattern(hourly_currents, n_days=2)

   inst.add_cable(cable, x=-0.3, depth=1.0, load=load_industrial)
   inst.add_cable(cable, x=0.0,  depth=1.0, load=load_base)
   inst.add_cable(cable, x=0.3,  depth=1.0, load=load_residential)

Cross-section visualisation
---------------------------

Use :func:`~thermal_cable_model.visualization.plot_cross_section` to visualise the
cable arrangement and annotate temperatures:

.. code-block:: python

   from thermal_cable_model.visualization import plot_cross_section

   fig = plot_cross_section(
       positions_x=[-0.3, 0.0, 0.3],
       depths=[1.0, 1.0, 1.0],
       outer_radii=[cable.outer_radius] * 3,
       temperatures=[result.max_conductor_temp(i) for i in range(3)],
       cable_names=["Left", "Centre", "Right"],
   )
   fig.savefig("parallel_cross_section.png", dpi=150)
