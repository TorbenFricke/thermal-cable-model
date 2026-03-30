"""Example: Flat vs trefoil formation — 3-phase MV single-core cables.

Compares two common installation geometries for a symmetric 3-phase
20 kV 240 mm² Cu XLPE cable system at 400 A per phase:

  1. **Flat formation** — three cables side-by-side (touching) at 1.0 m depth.
  2. **Trefoil formation** — three cables in an equilateral triangle
     (touching), centroid at 1.0 m depth.

Both carry the same balanced load.  The comparison shows how the
geometric arrangement affects mutual heating, peak conductor temperature,
and the thermal asymmetry between phases.
"""

import math

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from thermal_cable_model.cable import Cable
from thermal_cable_model.ground import KasudaModel
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import SOIL_STANDARD
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation
from thermal_cable_model.visualization import plot_cross_section, plot_temperature_history


def main():
    # ── Cable & load ──────────────────────────────────────────────────
    cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
    I_phase = 400.0  # A, balanced three-phase
    duration_s = 48 * 3600
    load = LoadProfile.constant(I_phase, duration_s)

    ground = KasudaModel(mean_surface_temp=12.0, annual_amplitude=10.0)
    depth = 1.0  # m — nominal burial depth

    d = 2.0 * cable.outer_radius  # centre-to-centre when touching
    OD_mm = cable.outer_radius * 2e3

    print(f"Cable: {cable.name}")
    print(f"  Outer diameter: {OD_mm:.1f} mm")
    print(f"  Centre-to-centre (touching): {d*1e3:.1f} mm")
    print(f"  Phase current: {I_phase:.0f} A (symmetric)")
    print()

    # ── Flat formation ────────────────────────────────────────────────
    #   L      C      R
    #   ●──d──●──d──●       all at depth = 1.0 m
    flat_x = [-d, 0.0, d]
    flat_d = [depth, depth, depth]

    inst_flat = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
    for x in flat_x:
        inst_flat.add_cable(cable, x=x, depth=depth, load=load)

    print("Running flat formation (48 h)...")
    result_flat = ThermalSimulation(inst_flat).run_transient(dt=300, duration=duration_s)

    flat_names = ["L (flat)", "C (flat)", "R (flat)"]
    for i, name in enumerate(flat_names):
        print(f"  {name}: max conductor = {result_flat.max_conductor_temp(i):.1f} °C, "
              f"max insulation = {result_flat.max_insulation_temp(i):.1f} °C")

    # ── Trefoil formation ─────────────────────────────────────────────
    #        ●  (top)           depth_top  = depth - d·√3/3
    #       / \
    #      ●───●  (bottom)     depth_bot  = depth + d·√3/6
    h_tri = d * math.sqrt(3) / 2.0
    depth_top = depth - h_tri * 2.0 / 3.0   # centroid at 'depth'
    depth_bot = depth + h_tri * 1.0 / 3.0

    tri_x = [-d / 2.0, d / 2.0, 0.0]
    tri_d = [depth_bot, depth_bot, depth_top]

    inst_tri = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
    for x, dp in zip(tri_x, tri_d):
        inst_tri.add_cable(cable, x=x, depth=dp, load=load)

    print("\nRunning trefoil formation (48 h)...")
    result_tri = ThermalSimulation(inst_tri).run_transient(dt=300, duration=duration_s)

    tri_names = ["BL (trefoil)", "BR (trefoil)", "Top (trefoil)"]
    for i, name in enumerate(tri_names):
        print(f"  {name}: max conductor = {result_tri.max_conductor_temp(i):.1f} °C, "
              f"max insulation = {result_tri.max_insulation_temp(i):.1f} °C")

    # ── Comparison summary ────────────────────────────────────────────
    flat_max = max(result_flat.max_conductor_temp(i) for i in range(3))
    flat_min = min(result_flat.max_conductor_temp(i) for i in range(3))
    tri_max = max(result_tri.max_conductor_temp(i) for i in range(3))
    tri_min = min(result_tri.max_conductor_temp(i) for i in range(3))

    print("\n── Flat vs Trefoil Comparison ──────────────────────────────")
    print(f"  Flat:    hottest phase = {flat_max:.1f} °C,  spread = {flat_max - flat_min:.1f} °C")
    print(f"  Trefoil: hottest phase = {tri_max:.1f} °C,  spread = {tri_max - tri_min:.1f} °C")
    diff = flat_max - tri_max
    better = "trefoil" if diff > 0 else "flat"
    print(f"  Trefoil vs flat hottest: {abs(diff):.1f} °C {'cooler' if diff > 0 else 'warmer'} "
          f"→ {better} is better for this layout")

    # ── Plots ─────────────────────────────────────────────────────────

    # 1. Cross-section comparison
    fig_cs, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    for ax, xs, ds, names, title in [
        (ax1, flat_x, flat_d, ["L", "C", "R"], "Flat Formation"),
        (ax2, tri_x, tri_d, ["BL", "BR", "Top"], "Trefoil Formation"),
    ]:
        ax.axhspan(0, 2.0, color="#d2b48c", alpha=0.3)
        ax.axhline(0, color="green", linewidth=2)
        res = result_flat if "Flat" in title else result_tri
        for i, (x, dp) in enumerate(zip(xs, ds)):
            circle = plt.Circle((x, dp), cable.outer_radius,
                                color="gray", ec="black", linewidth=1.5)
            ax.add_patch(circle)
            conductor = plt.Circle((x, dp), cable.outer_radius * 0.4,
                                   color="#b87333", ec="black", linewidth=0.8)
            ax.add_patch(conductor)
            T = res.max_conductor_temp(i)
            ax.annotate(
                f"{names[i]}\n{T:.1f} °C",
                xy=(x, dp), xytext=(x, dp - cable.outer_radius - 0.07),
                fontsize=9, ha="center", va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray"),
                arrowprops=dict(arrowstyle="->", color="gray"),
            )
        ax.set_xlim(-0.15, 0.15)
        ax.set_ylim(1.15, 0.85)
        ax.set_xlabel("Horizontal position [m]")
        ax.set_ylabel("Depth [m]")
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)

    fig_cs.tight_layout()
    fig_cs.savefig("flat_vs_trefoil_cross_section.png", dpi=150)
    print("\nSaved: flat_vs_trefoil_cross_section.png")

    # 2. Conductor temperature comparison
    t_hours = result_flat.times / 3600.0

    fig_t, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    colors_flat = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    for i, (name, c) in enumerate(zip(["L", "C", "R"], colors_flat)):
        axes[0].plot(t_hours, result_flat.conductor_temps[:, i],
                     linewidth=1.8, color=c, label=name)
    axes[0].set_ylabel("Temperature [°C]")
    axes[0].set_title("Flat Formation — Conductor Temperatures")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    colors_tri = ["#d62728", "#9467bd", "#8c564b"]
    for i, (name, c) in enumerate(zip(["BL", "BR", "Top"], colors_tri)):
        axes[1].plot(t_hours, result_tri.conductor_temps[:, i],
                     linewidth=1.8, color=c, label=name)
    axes[1].set_ylabel("Temperature [°C]")
    axes[1].set_xlabel("Time [hours]")
    axes[1].set_title("Trefoil Formation — Conductor Temperatures")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    fig_t.suptitle(f"3×{cable.name} @ {I_phase:.0f} A symmetric", fontsize=13, y=1.01)
    fig_t.tight_layout()
    fig_t.savefig("flat_vs_trefoil_temperatures.png", dpi=150, bbox_inches="tight")
    print("Saved: flat_vs_trefoil_temperatures.png")


if __name__ == "__main__":
    main()
