"""Example: Single MV cable — steady-state and transient analysis.

Demonstrates a 240 mm² Cu XLPE 20 kV cable buried at 1.2 m depth
carrying a constant 400 A load with seasonal ground temperature variation.
"""

import matplotlib
matplotlib.use("Agg")

import numpy as np

from thermal_cable_model.cable import Cable
from thermal_cable_model.ground import KasudaModel
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import SOIL_STANDARD
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation
from thermal_cable_model.visualization import plot_temperature_history, plot_cross_section


def main():
    # ── Cable definition ─────────────────────────────────────────────
    cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
    print(f"Cable: {cable.name}")
    print(f"  Outer diameter: {cable.outer_diameter * 1e3:.1f} mm")
    print(f"  R_ac(20°C): {cable.ac_resistance_20c * 1e3:.4f} mΩ/m")
    print(f"  R_ac(90°C): {cable.ac_resistance(90.0) * 1e3:.4f} mΩ/m")
    print()

    # ── Ground temperature — Central European climate ────────────────
    ground = KasudaModel(
        mean_surface_temp=10.0,
        annual_amplitude=12.0,
        day_of_min_surface_temp=35.0,    # early February
        soil_diffusivity=0.5e-6,
    )

    # ── Load profile — constant 400 A for one year ───────────────────
    one_year_s = 365.25 * 24 * 3600
    load = LoadProfile.constant(400.0, one_year_s)

    # ── Assembly ─────────────────────────────────────────────────────
    inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
    inst.add_cable(cable, x=0.0, depth=1.2, load=load)

    sim = ThermalSimulation(inst)

    # ── Steady-state at mid-summer (day 200) ─────────────────────────
    t_summer = 200 * 86400
    ss = sim.run_steady_state(time_s=t_summer)
    print("Steady-state at mid-summer (day 200):")
    print(f"  Conductor temperature: {ss.conductor_temps[0, 0]:.1f} °C")
    print(f"  Insulation temperature: {ss.insulation_temps[0, 0]:.1f} °C")
    print(f"  Cable surface temperature: {ss.surface_temps[0, 0]:.1f} °C")
    print(f"  Soil near cable: {ss.soil_temps[0, 0]:.1f} °C")
    print(f"  Ambient at depth: {ss.ambient_temps[0, 0]:.1f} °C")
    print()

    # ── Transient — one full year at 1-hour steps ────────────────────
    print("Running transient simulation (1 year, hourly steps)...")
    result = sim.run_transient(dt=3600, duration=one_year_s)
    print(f"  Max conductor temperature: {result.max_conductor_temp(0):.1f} °C")
    print(f"  Max insulation temperature: {result.max_insulation_temp(0):.1f} °C")

    # ── Plots ────────────────────────────────────────────────────────
    fig1 = plot_temperature_history(result, time_unit="days")
    fig1.savefig("single_cable_temperatures.png", dpi=150)
    print("\nSaved: single_cable_temperatures.png")

    fig2 = plot_cross_section(
        positions_x=[0.0],
        depths=[1.2],
        outer_radii=[cable.outer_radius],
        temperatures=[result.max_conductor_temp(0)],
        cable_names=[cable.name],
    )
    fig2.savefig("single_cable_cross_section.png", dpi=150)
    print("Saved: single_cable_cross_section.png")


if __name__ == "__main__":
    main()
