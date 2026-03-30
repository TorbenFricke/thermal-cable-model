"""
Thermal FEM — Cable thermal rating analysis.

Supports low-voltage and medium-voltage cables with:
  - Dynamic (time-varying) load profiles
  - Seasonally varying ground and ambient temperatures
  - Parallel cable installations with mutual heating
  - Cable crossings at arbitrary angles

Typical workflow
----------------
>>> from thermal_cable_model import Cable, CableInstallation, ThermalSimulation
>>> cable = Cable.single_core_xlpe_cu(240e-6, voltage_class="MV")
>>> inst = CableInstallation()
>>> inst.add_cable(cable, x=0.0, depth=1.2, load=LoadProfile.constant(400, 8760))
>>> sim = ThermalSimulation(inst)
>>> results = sim.run_transient(dt=3600, duration=8760*3600)
>>> results.plot_temperature_history()
"""

from thermal_cable_model.materials import ThermalMaterial, MATERIALS
from thermal_cable_model.cable import Cable, CableLayer
from thermal_cable_model.ground import GroundTemperatureModel, KasudaModel
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.thermal_network import CableThermalNetwork
from thermal_cable_model.solver import TransientSolver
from thermal_cable_model.crossing import CableCrossing
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation
from thermal_cable_model.visualization import plot_temperature_history

__version__ = "0.1.0"
