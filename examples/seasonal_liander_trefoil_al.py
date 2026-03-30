"""Example: Full-year trefoil 240 mm² Al XLPE with Liander 15-minute load shape.

Uses open *historical 15-minute operational measurements* (e.g. OS Apeldoorn,
station installation series) from Liander.  The CSV ``load`` column is **not**
in amperes; we take its magnitude as a **dimensionless proxy** for how demand
varies in time and map it linearly to a plausible phase-current band for a
buried MV circuit.

Bundled sample data (relative to this file via ``pathlib``):

  ``Path(__file__).resolve().parent / "data" / "OS_Apeldoorn.csv"``

Original source (same dataset; check Liander’s disclaimer before redistribution):

  https://www.liander.nl/over-ons/open-data#historische-15-minuten-bedrijfsmetingen

Some files omit ``load`` for a few rows at DST boundaries; those values are
forward-filled from the previous sample.

Tuning
------
* ``I_PEAK_A`` / ``I_MIN_A`` — scale the mapped current (default peak **360 A**
  is a reasonable order of magnitude for continuous loading of 240 mm² Al XLPE
  MV in soil; adjust for your rating case).
* ``KasudaModel`` parameters — Netherlands-oriented mean temperature, annual
  amplitude, and day of minimum surface temperature.

Soil at burial depth follows the Kasuda attenuation/phase shift from these
surface parameters.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np

from thermal_cable_model.cable import Cable
from thermal_cable_model.ground import KasudaModel
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import SOIL_STANDARD
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation
from thermal_cable_model.visualization import plot_temperature_history

# This script’s directory and bundled Liander sample CSV (pathlib, portable).
_EXAMPLES_DIR = Path(__file__).resolve().parent
LIANDER_CSV_PATH = _EXAMPLES_DIR / "data" / "OS_Apeldoorn.csv"

# Phase current scaling (240 mm² Al XLPE MV — adjust for your study).
I_PEAK_A = 360.0
I_MIN_A = 0.22 * I_PEAK_A  # ~22% of peak at troughs

# Burial: trefoil centroid depth [m] (matches flat_vs_trefoil-style layout).
DEPTH_CENTROID_M = 1.2


def load_liander_station_csv(filepath: Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse Liander ``datetime,load`` CSV into t [s] from first row and |load|."""
    times_s: list[float] = []
    magnitudes: list[float] = []
    t0: datetime | None = None

    last_mag: float | None = None

    with filepath.open(newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError(f"Expected header with at least 2 columns in {filepath}")

        for row in reader:
            if len(row) < 2:
                continue
            raw_dt = row[0].strip()
            if not raw_dt:
                continue
            # Liander uses "2024-01-01 00:00:00+00:00"
            dt = datetime.fromisoformat(raw_dt.replace("Z", "+00:00"))
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
            else:
                dt = dt.astimezone(timezone.utc)

            raw_load = row[1].strip()
            if raw_load == "":
                if last_mag is None:
                    continue
                mag_v = last_mag
            else:
                mag_v = abs(float(raw_load))
                last_mag = mag_v

            if t0 is None:
                t0 = dt
            times_s.append((dt - t0).total_seconds())
            magnitudes.append(mag_v)

    if not times_s:
        raise ValueError(f"No data rows in {filepath}")

    t = np.asarray(times_s, dtype=float)
    mag = np.asarray(magnitudes, dtype=float)

    if not np.all(np.diff(t) > 0):
        raise ValueError("Timestamps must be strictly increasing")

    dt_expected = 900.0
    gaps = np.diff(t)
    if not np.allclose(gaps, dt_expected, rtol=0.0, atol=1.0):
        bad = np.where(np.abs(gaps - dt_expected) > 1.0)[0]
        print(
            f"Warning: {len(bad)} interval(s) differ from {dt_expected:.0f} s "
            f"(first at row index {bad[0] if len(bad) else 'n/a'})",
            file=sys.stderr,
        )

    return t, mag


def magnitudes_to_load_profile(
    times_s: np.ndarray,
    mag: np.ndarray,
    i_peak: float,
    i_min: float,
) -> LoadProfile:
    """Linear map mag to [i_min, i_peak] (shape only; not physical substation amps)."""
    lo, hi = float(mag.min()), float(mag.max())
    if hi <= lo:
        currents = np.full_like(mag, (i_peak + i_min) / 2.0)
    else:
        currents = i_min + (i_peak - i_min) * (mag - lo) / (hi - lo)
    return LoadProfile(times_s, currents)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Seasonal trefoil Al thermal run driven by Liander 15-min CSV.",
    )
    parser.add_argument(
        "csv_path",
        nargs="?",
        type=Path,
        default=None,
        help=(
            "Path to Liander CSV (optional; default: examples/data/OS_Apeldoorn.csv "
            "next to this script)"
        ),
    )
    args = parser.parse_args()
    if args.csv_path is None:
        csv_path = LIANDER_CSV_PATH
    else:
        p = args.csv_path.expanduser()
        csv_path = p if p.is_absolute() else (Path.cwd() / p).resolve()

    if not csv_path.is_file():
        print(
            f"CSV not found: {csv_path}\n"
            f"Expected bundled file at:\n  {LIANDER_CSV_PATH}\n"
            f"Or download from Liander open data (historische 15-minuten bedrijfsmetingen) "
            f"and pass the path as an argument.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Loading: {csv_path}")
    t_raw, mag = load_liander_station_csv(csv_path)
    load = magnitudes_to_load_profile(t_raw, mag, I_PEAK_A, I_MIN_A)
    currents_scaled = load.current_array(t_raw)
    print(
        f"Scaled phase current [A]: min={currents_scaled.min():.1f}, "
        f"median={float(np.median(currents_scaled)):.1f}, max={currents_scaled.max():.1f}"
    )
    print(f"Profile span: {t_raw[0]:.0f} … {t_raw[-1]:.0f} s ({len(t_raw)} samples)\n")

    base = Cable.single_core_xlpe_al(240, voltage_class="MV", voltage_kv=20.0)
    cable_bl = replace(base, name=f"{base.name} — BL")
    cable_br = replace(base, name=f"{base.name} — BR")
    cable_top = replace(base, name=f"{base.name} — top")
    print(f"Cable: {base.name} (trefoil, balanced load)")

    # Netherlands-oriented Kasuda surface parameters (soil temp at depth follows model).
    ground = KasudaModel(
        mean_surface_temp=10.7,
        annual_amplitude=10.5,
        day_of_min_surface_temp=40.0,
        soil_diffusivity=0.5e-6,
    )

    d = 2.0 * base.outer_radius
    h_tri = d * math.sqrt(3) / 2.0
    depth_top = DEPTH_CENTROID_M - h_tri * 2.0 / 3.0
    depth_bot = DEPTH_CENTROID_M + h_tri * 1.0 / 3.0
    tri_x = [-d / 2.0, d / 2.0, 0.0]
    tri_d = [depth_bot, depth_bot, depth_top]

    inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
    for cab, x, dp in zip((cable_bl, cable_br, cable_top), tri_x, tri_d):
        inst.add_cable(cab, x=x, depth=dp, load=load)

    sim = ThermalSimulation(inst)

    duration = float(t_raw[-1] + 900.0)
    dt = 900.0
    print(
        f"Running transient: duration={duration/86400:.2f} d, dt={dt:.0f} s "
        f"(~{int(duration / dt)} steps)..."
    )

    result = sim.run_transient(dt=dt, duration=duration)

    n_c = inst.n_cables
    for i in range(n_c):
        print(
            f"  Cable {i} ({inst.cables[i].name}): "
            f"max conductor {result.max_conductor_temp(i):.1f} °C, "
            f"max insulation {result.max_insulation_temp(i):.1f} °C"
        )

    peak_step = int(np.argmax(result.conductor_temps[:, 0]))
    peak_day = result.times[peak_step] / 86400
    print(f"  Peak (phase 0) at day {peak_day:.1f} ({_day_to_month(peak_day)})")
    print(f"  Ambient at peak (cable 0): {result.ambient_temps[peak_step, 0]:.1f} °C\n")

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
            print(f"  {season}: conductor (phase 0) "
                  f"{np.mean(tc):.1f} avg / {np.max(tc):.1f} max °C")

    out_png = _EXAMPLES_DIR / "seasonal_liander_trefoil_al_temperatures.png"
    fig = plot_temperature_history(result, time_unit="days")
    fig.savefig(out_png, dpi=150)
    print(f"\nSaved: {out_png}")


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
