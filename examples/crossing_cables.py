"""Example: Cable crossing — MV cable over an LV cable at 60°.

Analyses the mutual temperature rise when a 20 kV MV cable crosses over
a 0.6 kV LV cable at an angle of 60° with a vertical separation of 0.3 m.

Includes a comparison with each cable simulated individually (no neighbour)
at the same depth and current to quantify the mutual heating from the crossing.
"""

import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt

from thermal_cable_model.cable import Cable
from thermal_cable_model.crossing import CableCrossing, crossing_derating_factor
from thermal_cable_model.ground import KasudaModel
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import SOIL_STANDARD
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation


def main():
    # ── Cable definitions ────────────────────────────────────────────
    mv_cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
    lv_cable = Cable.three_core_xlpe_cu(150, voltage_class="LV", voltage_kv=0.6)

    print(f"Upper cable: {mv_cable.name}")
    print(f"Lower cable: {lv_cable.name}")
    print()

    depth_upper = 0.9   # MV cable at 0.9 m
    depth_lower = 1.2   # LV cable at 1.2 m
    crossing_angle = 60  # degrees

    I_mv = 350.0  # A
    I_lv = 280.0  # A

    # ── Crossing analysis ────────────────────────────────────────────
    crossing = CableCrossing(
        cable_upper=mv_cable,
        cable_lower=lv_cable,
        depth_upper=depth_upper,
        depth_lower=depth_lower,
        crossing_angle_deg=crossing_angle,
        soil=SOIL_STANDARD,
    )

    # Steady-state mutual heating
    W_mv = mv_cable.total_heat_per_length(I_mv, 70.0)
    W_lv = lv_cable.total_heat_per_length(I_lv, 55.0)

    dT_at_mv = crossing.temperature_rise_at_upper(W_lv)
    dT_at_lv = crossing.temperature_rise_at_lower(W_mv)

    print(f"Crossing angle: {crossing_angle}°")
    print(f"Vertical separation: {abs(depth_lower - depth_upper)*1e3:.0f} mm")
    print(f"MV heat rate: {W_mv:.2f} W/m at {I_mv:.0f} A")
    print(f"LV heat rate: {W_lv:.2f} W/m at {I_lv:.0f} A")
    print()
    print(f"Steady-state ΔT at MV cable (from LV): {dT_at_mv:.2f} °C")
    print(f"Steady-state ΔT at LV cable (from MV): {dT_at_lv:.2f} °C")
    print()

    # ── Derating factor ──────────────────────────────────────────────
    df_mv = crossing_derating_factor(mv_cable, depth_upper, dT_at_mv, SOIL_STANDARD)
    df_lv = crossing_derating_factor(lv_cable, depth_lower, dT_at_lv, SOIL_STANDARD)
    print(f"Derating factor for MV cable: {df_mv:.3f}")
    print(f"Derating factor for LV cable: {df_lv:.3f}")
    print()

    # ── Transient crossing temperature rise ──────────────────────────
    times_h = np.array([0.5, 1, 2, 4, 8, 12, 24, 48, 100, 200])
    times_s = times_h * 3600
    print("Transient ΔT at MV cable from LV heat source:")
    for th, ts in zip(times_h, times_s):
        dT = crossing.transient_temperature_rise_at_upper(W_lv, ts)
        print(f"  t = {th:6.1f} h  →  ΔT = {dT:.2f} °C")

    # ── Common simulation parameters ─────────────────────────────────
    ground = KasudaModel(mean_surface_temp=12.0, annual_amplitude=10.0)
    duration_s = 24 * 3600
    load_mv = LoadProfile.constant(I_mv, duration_s)
    load_lv = LoadProfile.constant(I_lv, duration_s)

    # ── Isolated cable simulations ────────────────────────────────────
    # Each cable is simulated alone so the thermal network does not
    # introduce spurious parallel mutual heating between them.
    print("\nRunning individual cable simulations (24 hours)...")

    inst_mv = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
    inst_mv.add_cable(mv_cable, x=0.0, depth=depth_upper, load=load_mv)
    result_mv = ThermalSimulation(inst_mv).run_transient(dt=300, duration=duration_s)

    inst_lv = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
    inst_lv.add_cable(lv_cable, x=0.0, depth=depth_lower, load=load_lv)
    result_lv = ThermalSimulation(inst_lv).run_transient(dt=300, duration=duration_s)

    print(f"  MV conductor temp (isolated): {result_mv.max_conductor_temp(0):.1f} °C")
    print(f"  LV conductor temp (isolated): {result_lv.max_conductor_temp(0):.1f} °C")

    # ── Superimpose crossing temperature rise ─────────────────────────
    # The crossing model is applied analytically on top of the isolated
    # results — this correctly captures only the localised interaction
    # at the crossing point.
    print("\nSuperimposing transient crossing effect...")

    mv_cond_isolated = result_mv.conductor_temps[:, 0]
    lv_cond_isolated = result_lv.conductor_temps[:, 0]
    mv_cond_crossing = mv_cond_isolated.copy()
    lv_cond_crossing = lv_cond_isolated.copy()

    for i, t in enumerate(result_mv.times):
        t_eff = max(t, 1.0)
        W_lv_t = lv_cable.total_heat_per_length(I_lv, lv_cond_isolated[i])
        W_mv_t = mv_cable.total_heat_per_length(I_mv, mv_cond_isolated[i])
        mv_cond_crossing[i] += crossing.transient_temperature_rise_at_upper(W_lv_t, t_eff)
        lv_cond_crossing[i] += crossing.transient_temperature_rise_at_lower(W_mv_t, t_eff)

    max_mv_cross = float(np.max(mv_cond_crossing))
    max_lv_cross = float(np.max(lv_cond_crossing))

    print(f"  MV conductor temp (with crossing): {max_mv_cross:.1f} °C")
    print(f"  LV conductor temp (with crossing): {max_lv_cross:.1f} °C")

    # ── Comparison summary ────────────────────────────────────────────
    diff_mv = max_mv_cross - result_mv.max_conductor_temp(0)
    diff_lv = max_lv_cross - result_lv.max_conductor_temp(0)
    print(f"\n── Crossing vs Isolated ────────────────────────────────────")
    print(f"  MV cable: crossing adds {diff_mv:+.1f} °C")
    print(f"  LV cable: crossing adds {diff_lv:+.1f} °C")

    # ── Comparison plot — conductor temperatures ──────────────────────
    t_hours = result_mv.times / 3600.0

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    axes[0].plot(t_hours, mv_cond_crossing,
                 linewidth=1.8, label="With crossing (60°)")
    axes[0].plot(t_hours, mv_cond_isolated,
                 linewidth=1.8, linestyle="--", label="Isolated")
    axes[0].set_ylabel("Temperature [°C]")
    axes[0].set_title(f"MV Cable — {mv_cable.name} @ {depth_upper} m, {I_mv:.0f} A")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t_hours, lv_cond_crossing,
                 linewidth=1.8, label="With crossing (60°)")
    axes[1].plot(t_hours, lv_cond_isolated,
                 linewidth=1.8, linestyle="--", label="Isolated")
    axes[1].set_ylabel("Temperature [°C]")
    axes[1].set_xlabel("Time [hours]")
    axes[1].set_title(f"LV Cable — {lv_cable.name} @ {depth_lower} m, {I_lv:.0f} A")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    fig.suptitle("Conductor Temperature: Crossing vs Isolated", fontsize=14)
    fig.tight_layout()
    fig.savefig("crossing_cables_temperatures.png", dpi=150)
    print("\nSaved: crossing_cables_temperatures.png")


if __name__ == "__main__":
    main()
