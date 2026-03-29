"""Example: Parallel cable installation — three MV circuits side-by-side.

Three single-core 20 kV 240 mm² Cu XLPE cables in flat formation at 1.0 m
depth with 0.3 m spacing.  Each carries a different load profile.
"""

import matplotlib
matplotlib.use("Agg")

import numpy as np

from thermal_fem.cable import Cable
from thermal_fem.ground import KasudaModel
from thermal_fem.loads import LoadProfile
from thermal_fem.materials import SOIL_WET
from thermal_fem.simulation import CableInstallation, ThermalSimulation
from thermal_fem.visualization import plot_temperature_history, plot_cross_section


def main():
    # ── Cables ───────────────────────────────────────────────────────
    cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)

    # ── Load profiles (48 h, 10-minute steps) ────────────────────────
    duration_s = 48 * 3600
    load_center = LoadProfile.constant(450.0, duration_s)

    # Left cable: cyclic industrial load
    load_left = LoadProfile.cyclic(
        peak_current=500.0,
        base_current=100.0,
        period_s=8 * 3600,
        duty_cycle=0.6,
        n_cycles=6,
    )

    # Right cable: daily residential pattern
    hourly = [
        120, 100, 90, 85, 80, 90, 150, 280,
        350, 320, 300, 290, 280, 270, 260, 270,
        310, 380, 420, 400, 350, 280, 200, 150,
    ]
    load_right = LoadProfile.daily_pattern(hourly, n_days=2)

    # ── Ground ───────────────────────────────────────────────────────
    ground = KasudaModel(
        mean_surface_temp=11.0,
        annual_amplitude=11.0,
        day_of_min_surface_temp=35.0,
        soil_diffusivity=0.6e-6,
    )
    # Start in July (day 200)
    t_offset = 200 * 86400

    # ── Installation ─────────────────────────────────────────────────
    inst = CableInstallation(soil=SOIL_WET, ground_temp_model=ground)
    inst.add_cable(cable, x=-0.3, depth=1.0, load=load_left)
    inst.add_cable(cable, x=0.0,  depth=1.0, load=load_center)
    inst.add_cable(cable, x=0.3,  depth=1.0, load=load_right)

    sim = ThermalSimulation(inst)

    # ── Transient ────────────────────────────────────────────────────
    print("Running 48-hour transient with 3 parallel cables...")
    result = sim.run_transient(dt=600, duration=duration_s)

    for i in range(3):
        names = ["Left (cyclic)", "Centre (constant)", "Right (residential)"]
        print(f"  {names[i]}:")
        print(f"    Max conductor temp: {result.max_conductor_temp(i):.1f} °C")
        print(f"    Max insulation temp: {result.max_insulation_temp(i):.1f} °C")

    # ── Plots ────────────────────────────────────────────────────────
    fig1 = plot_temperature_history(result, time_unit="hours")
    fig1.savefig("parallel_cables_temperatures.png", dpi=150)
    print("\nSaved: parallel_cables_temperatures.png")

    fig2 = plot_cross_section(
        positions_x=[-0.3, 0.0, 0.3],
        depths=[1.0, 1.0, 1.0],
        outer_radii=[cable.outer_radius] * 3,
        temperatures=[
            result.max_conductor_temp(0),
            result.max_conductor_temp(1),
            result.max_conductor_temp(2),
        ],
        cable_names=["Left", "Centre", "Right"],
    )
    fig2.savefig("parallel_cables_cross_section.png", dpi=150)
    print("Saved: parallel_cables_cross_section.png")


if __name__ == "__main__":
    main()
