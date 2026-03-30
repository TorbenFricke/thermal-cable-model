"""Plotting utilities for thermal simulation results."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from thermal_cable_model.solver import TransientResult


def plot_temperature_history(
    result: TransientResult,
    cable_indices: list[int] | None = None,
    show_ambient: bool = True,
    show_current: bool = True,
    time_unit: str = "hours",
    figsize: tuple[float, float] = (14, 8),
) -> Figure:
    """Plot conductor, insulation, surface, and soil temperatures vs time.

    Parameters
    ----------
    result : TransientResult
    cable_indices : list[int], optional
        Which cables to plot.  Defaults to all.
    show_ambient : bool
        Overlay the ambient ground temperature.
    show_current : bool
        Add a secondary axis with the load current.
    time_unit : {"seconds", "hours", "days"}
    figsize : tuple
    """
    divisor = {"seconds": 1.0, "hours": 3600.0, "days": 86400.0}[time_unit]
    t = result.times / divisor

    if cable_indices is None:
        cable_indices = list(range(result.n_cables))

    n_panels = 1 + int(show_current)
    fig, axes = plt.subplots(
        n_panels, 1, figsize=figsize, sharex=True,
        gridspec_kw={"height_ratios": [3, 1] if n_panels == 2 else [1]},
    )
    if n_panels == 1:
        axes = [axes]
    ax_temp = axes[0]

    colors = plt.cm.tab10.colors  # type: ignore[attr-defined]

    for ci in cable_indices:
        c = colors[ci % len(colors)]
        label = result.cable_names[ci] if result.cable_names else f"Cable {ci}"
        ax_temp.plot(t, result.conductor_temps[:, ci], color=c,
                     linewidth=1.8, label=f"{label} — conductor")
        ax_temp.plot(t, result.insulation_temps[:, ci], color=c,
                     linewidth=1.4, linestyle="--", label=f"{label} — insulation")
        ax_temp.plot(t, result.sheath_temps[:, ci], color=c,
                     linewidth=1.1, linestyle=(0, (5, 2)),
                     label=f"{label} — sheath")
        ax_temp.plot(t, result.armour_temps[:, ci], color=c,
                     linewidth=0.9, linestyle=(0, (3, 2, 1, 2)),
                     label=f"{label} — armour")
        ax_temp.plot(t, result.surface_temps[:, ci], color=c,
                     linewidth=0.8, linestyle=":", label=f"{label} — surface")
        ax_temp.plot(t, result.soil_temps[:, ci], color=c,
                     linewidth=0.7, linestyle="-.", label=f"{label} — soil")
        if show_ambient:
            ax_temp.plot(t, result.ambient_temps[:, ci], color=c,
                         linewidth=0.6, linestyle=(0, (5, 10)),
                         alpha=0.5, label=f"{label} — ambient")

    ax_temp.set_ylabel("Temperature [°C]")
    ax_temp.legend(fontsize=7, ncol=2, loc="upper right")
    ax_temp.grid(True, alpha=0.3)
    ax_temp.set_title("Cable Temperature History")

    if show_current and n_panels == 2:
        ax_cur = axes[1]
        for ci in cable_indices:
            c = colors[ci % len(colors)]
            label = result.cable_names[ci] if result.cable_names else f"Cable {ci}"
            ax_cur.plot(t, result.currents[:, ci], color=c, linewidth=1.2,
                        label=label)
        ax_cur.set_ylabel("Current [A]")
        ax_cur.set_xlabel(f"Time [{time_unit}]")
        ax_cur.legend(fontsize=7)
        ax_cur.grid(True, alpha=0.3)
    else:
        ax_temp.set_xlabel(f"Time [{time_unit}]")

    fig.tight_layout()
    return fig


def plot_cross_section(
    positions_x: list[float],
    depths: list[float],
    outer_radii: list[float],
    temperatures: list[float] | None = None,
    cable_names: list[str] | None = None,
    soil_extent: tuple[float, float, float] = (-1.0, 1.0, 2.5),
    figsize: tuple[float, float] = (10, 8),
) -> Figure:
    """Plot a 2-D cross-section of the cable arrangement.

    Parameters
    ----------
    positions_x : list[float]
        Horizontal positions [m].
    depths : list[float]
        Burial depths [m] (positive downward).
    outer_radii : list[float]
        Cable outer radii [m].
    temperatures : list[float], optional
        Conductor temperatures to annotate.
    cable_names : list[str], optional
    soil_extent : (x_min, x_max, y_max_depth)
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    x_min, x_max, y_max = soil_extent

    # Soil background
    ax.axhspan(0, y_max, color="#d2b48c", alpha=0.3, label="Soil")
    ax.axhline(0, color="green", linewidth=2, label="Ground surface")

    for i, (x, d, r) in enumerate(zip(positions_x, depths, outer_radii)):
        circle = plt.Circle((x, d), r, color="gray", ec="black", linewidth=1.5)
        ax.add_patch(circle)
        # conductor
        r_cond = r * 0.4
        conductor = plt.Circle((x, d), r_cond, color="#b87333", ec="black",
                               linewidth=0.8)
        ax.add_patch(conductor)

        name = cable_names[i] if cable_names else f"Cable {i}"
        label_text = name
        if temperatures is not None:
            label_text += f"\n{temperatures[i]:.1f} °C"
        ax.annotate(
            label_text,
            xy=(x, d),
            xytext=(x, d - r - 0.08),
            fontsize=8,
            ha="center",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray"),
            arrowprops=dict(arrowstyle="->", color="gray"),
        )

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_max, -0.2)
    ax.set_xlabel("Horizontal position [m]")
    ax.set_ylabel("Depth [m]")
    ax.set_title("Cable Cross-Section")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


