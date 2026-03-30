"""Example: Seasonal simulation with dynamic daily load pattern.

Runs a full year with a repeating daily load pattern superimposed on
seasonal ground temperature variation.  Demonstrates how conductor
temperature varies with both load and ambient conditions.
"""

import matplotlib
matplotlib.use("Agg")

import numpy as np

from thermal_cable_model.cable import Cable
from thermal_cable_model.ground import KasudaModel
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import SOIL_STANDARD
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation
from thermal_cable_model.visualization import plot_temperature_history


def main():
    cable = Cable.single_core_xlpe_cu(300, voltage_class="MV", voltage_kv=20.0)
    print(f"Cable: {cable.name}")

    # ── Seasonal ground temperature ──────────────────────────────────
    ground = KasudaModel(
        mean_surface_temp=10.5,
        annual_amplitude=13.0,
        day_of_min_surface_temp=30.0,
        soil_diffusivity=0.5e-6,
    )

    # ── Daily load pattern (repeated for 365 days) ───────────────────
    hourly_current = [
        180, 160, 140, 130, 130, 150, 220, 350,
        420, 440, 460, 450, 430, 420, 400, 410,
        440, 480, 500, 480, 420, 350, 280, 220,
    ]
    load = LoadProfile.daily_pattern(hourly_current, n_days=365)

    print(f"Peak daily current: {max(hourly_current)} A")
    print(f"Min daily current:  {min(hourly_current)} A")
    print()

    # ── Installation ─────────────────────────────────────────────────
    inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
    inst.add_cable(cable, x=0.0, depth=1.2, load=load)

    sim = ThermalSimulation(inst)

    # ── Run year-long transient ──────────────────────────────────────
    one_year = 365 * 24 * 3600
    print("Running 1-year transient simulation (hourly steps)...")
    result = sim.run_transient(dt=3600, duration=one_year)

    print(f"  Max conductor temperature: {result.max_conductor_temp(0):.1f} °C")
    print(f"  Max insulation temperature: {result.max_insulation_temp(0):.1f} °C")

    # Find when peak occurs
    peak_step = np.argmax(result.conductor_temps[:, 0])
    peak_day = result.times[peak_step] / 86400
    print(f"  Peak occurs at day {peak_day:.0f} "
          f"({_day_to_month(peak_day)})")
    print(f"  Ambient at peak: {result.ambient_temps[peak_step, 0]:.1f} °C")
    print()

    # ── Seasonal summary ─────────────────────────────────────────────
    days = result.times / 86400
    for season, (d_start, d_end) in {
        "Winter (Dec–Feb)": (335, 365 + 59),
        "Spring (Mar–May)": (59, 151),
        "Summer (Jun–Aug)": (151, 243),
        "Autumn (Sep–Nov)": (243, 335),
    }.items():
        mask = (days >= d_start) & (days < d_end) if d_end <= 365 else (
            (days >= d_start) | (days < d_end - 365)
        )
        if mask.any():
            tc = result.conductor_temps[mask, 0]
            print(f"  {season}: "
                  f"conductor {np.mean(tc):.1f} avg / {np.max(tc):.1f} max °C")

    fig = plot_temperature_history(result, time_unit="days")
    fig.savefig("seasonal_dynamic_temperatures.png", dpi=150)
    print("\nSaved: seasonal_dynamic_temperatures.png")


def _day_to_month(day: float) -> str:
    months = [
        (31, "Jan"), (59, "Feb"), (90, "Mar"), (120, "Apr"),
        (151, "May"), (181, "Jun"), (212, "Jul"), (243, "Aug"),
        (273, "Sep"), (304, "Oct"), (334, "Nov"), (365, "Dec"),
    ]
    d = int(day) % 365
    for end, name in months:
        if d < end:
            return name
    return "Dec"


if __name__ == "__main__":
    main()
