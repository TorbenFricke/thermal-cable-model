Single Cable
============

**Script:** ``examples/single_cable.py``

Demonstrates the complete workflow for a single buried cable: defining the
cable, setting up the ground model, running steady-state and transient
analyses, and plotting the results.

Setup
-----

- **Cable:** 1 × 240 mm\ :sup:`2` Cu XLPE, 20 kV
- **Load:** constant 400 A for one year
- **Soil:** standard soil (ρ = 1.0 K·m/W)
- **Ground temperature:** Kasuda model (central European climate, mean 10 °C,
  amplitude 12 °C)
- **Burial depth:** 1.2 m

What it demonstrates
--------------------

1. Creating a cable with :meth:`Cable.single_core_xlpe_cu
   <thermal_cable_model.Cable.single_core_xlpe_cu>`
2. Setting up a :class:`~thermal_cable_model.KasudaModel` for seasonal ground
   temperature variation
3. Running a **steady-state** analysis at mid-summer (day 200)
4. Running a **year-long transient** at 1-hour time steps
5. Plotting temperature history and cross-section

Key code
--------

.. code-block:: python

   cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
   ground = KasudaModel(mean_surface_temp=10.0, annual_amplitude=12.0)
   load = LoadProfile.constant(400.0, 365.25 * 24 * 3600)

   inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
   inst.add_cable(cable, x=0.0, depth=1.2, load=load)

   sim = ThermalSimulation(inst)
   result = sim.run_transient(dt=3600, duration=365.25 * 24 * 3600)

Output files
------------

- ``single_cable_temperatures.png`` — temperature history over one year
- ``single_cable_cross_section.png`` — cable cross-section with peak temperature
