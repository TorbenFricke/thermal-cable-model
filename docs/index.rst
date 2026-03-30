Thermal FEM
===========

**Cable thermal rating analysis for buried power cables.**

Thermal FEM is a Python library for computing the temperature of buried
low-voltage and medium-voltage power cables under steady-state and transient
loading conditions.  It implements the lumped-parameter thermal network approach
from IEC 60287 / IEC 60853.

.. code-block:: python

   from thermal_cable_model import Cable, CableInstallation, ThermalSimulation, LoadProfile

   cable = Cable.single_core_xlpe_cu(240, voltage_class="MV")
   inst = CableInstallation()
   inst.add_cable(cable, x=0.0, depth=1.2,
                  load=LoadProfile.constant(400, 48 * 3600))
   sim = ThermalSimulation(inst)
   result = sim.run_transient(dt=60)
   print(f"Max conductor temperature: {result.max_conductor_temp(0):.1f} °C")

Key Features
------------

- **Multi-layer cable models** with factory methods for common LV/MV constructions
- **Dynamic load profiles** — constant, cyclic, daily patterns, or CSV import
- **Seasonal ground temperature** via the Kasuda & Archenbach equation
- **Parallel cables** with mutual heating (IEC 60287 image method)
- **Cable crossings** at arbitrary angles with derating factors (CIGRE TB 640)
- **Transient solver** — implicit Euler with Picard iteration for temperature-dependent resistance

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   getting_started
   user_guide/index
   examples/index
   theory
   api/index
   changelog
