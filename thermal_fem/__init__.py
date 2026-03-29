"""
Thermal FEM — Cable thermal rating and temperature field analysis.

Supports low-voltage and medium-voltage cables with:
  - Dynamic (time-varying) load profiles
  - Seasonally varying ground and ambient temperatures
  - Parallel cable installations with mutual heating
  - Cable crossings at arbitrary angles
  - 2-D finite-element soil temperature field visualisation

Typical workflow
----------------
>>> from thermal_fem import Cable, CableInstallation, ThermalSimulation
>>> cable = Cable.single_core_xlpe_cu(240e-6, voltage_class="MV")
>>> inst = CableInstallation()
>>> inst.add_cable(cable, x=0.0, depth=1.2, load=LoadProfile.constant(400, 8760))
>>> sim = ThermalSimulation(inst)
>>> results = sim.run_transient(dt=3600, duration=8760*3600)
>>> results.plot_temperature_history()
"""

from thermal_fem.materials import ThermalMaterial, MATERIALS
from thermal_fem.cable import Cable, CableLayer
from thermal_fem.ground import GroundTemperatureModel, KasudaModel
from thermal_fem.loads import LoadProfile
from thermal_fem.thermal_network import CableThermalNetwork
from thermal_fem.solver import TransientSolver
from thermal_fem.crossing import CableCrossing
from thermal_fem.simulation import CableInstallation, ThermalSimulation
from thermal_fem.visualization import plot_temperature_history, plot_soil_temperature_field

__version__ = "0.1.0"
