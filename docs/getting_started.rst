Getting Started
===============

This tutorial walks through a complete thermal analysis of a single buried
medium-voltage cable — from defining the cable and load, through running a
transient simulation, to inspecting and plotting the results.

Step 1: Define a cable
----------------------

Thermal FEM includes factory methods for common cable constructions.  A
single-core 240 mm\ :sup:`2` copper XLPE cable rated at 20 kV is created with:

.. code-block:: python

   from thermal_fem import Cable

   cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
   print(cable.name)                      # "1×240 mm² Cu XLPE 20 kV"
   print(f"{cable.outer_diameter*1e3:.1f} mm")  # outer diameter in mm

The factory method builds the insulation, screen, and jacket layers
automatically.  See :doc:`user_guide/cables` for how to define custom cables
layer-by-layer.

Step 2: Define a load profile
-----------------------------

A load profile describes the cable current as a function of time.  The
simplest case is a constant current:

.. code-block:: python

   from thermal_fem import LoadProfile

   duration = 48 * 3600       # 48 hours in seconds
   load = LoadProfile.constant(400.0, duration)  # 400 A for 48 h

Other load shapes are available — see :doc:`user_guide/loads`.

Step 3: Set up the installation
-------------------------------

A :class:`~thermal_fem.CableInstallation` describes the physical arrangement:
cable positions, burial depths, soil properties, and the ground temperature
model.

.. code-block:: python

   from thermal_fem import CableInstallation, KasudaModel
   from thermal_fem.materials import SOIL_STANDARD

   ground = KasudaModel(
       mean_surface_temp=10.0,
       annual_amplitude=12.0,
       day_of_min_surface_temp=35.0,
       soil_diffusivity=0.5e-6,
   )

   inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
   inst.add_cable(cable, x=0.0, depth=1.2, load=load)

Step 4: Run the simulation
--------------------------

Wrap the installation in a :class:`~thermal_fem.ThermalSimulation` and run
either a steady-state or transient analysis:

.. code-block:: python

   from thermal_fem import ThermalSimulation

   sim = ThermalSimulation(inst)

   # Steady-state at a particular moment in time
   ss = sim.run_steady_state(time_s=200 * 86400)  # mid-summer
   print(f"Conductor: {ss.conductor_temps[0, 0]:.1f} °C")

   # Transient — 48 hours at 1-minute steps
   result = sim.run_transient(dt=60, duration=48 * 3600)
   print(f"Max conductor temp: {result.max_conductor_temp(0):.1f} °C")

Step 5: Visualise the results
-----------------------------

Built-in plotting functions create publication-ready figures:

.. code-block:: python

   from thermal_fem import plot_temperature_history

   fig = plot_temperature_history(result, time_unit="hours")
   fig.savefig("my_cable_temperatures.png", dpi=150)

The :func:`~thermal_fem.plot_temperature_history` function shows conductor,
insulation, surface, and soil node temperatures on one axis with an optional
current panel below.

What next?
----------

- **Multiple cables** — see :doc:`user_guide/parallel_cables`
- **Cable crossings** — see :doc:`user_guide/crossings`
- **2-D soil temperature field** — see :doc:`user_guide/fem`
- **Full API reference** — see :doc:`api/index`
- **Worked examples** — see :doc:`examples/index`
