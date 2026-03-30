"""Reference-oriented tests for parallel burial and cable crossings.

Inventory: CIGRE TB 880 parallel cases vs this package
=======================================================
**Case 1-1** (132 kV single-core, trefoil touching, solid bonding) and **Case 1-2**
(flat spaced, cross-bonding, mutual-heating variants) in *CIGRE TB 880* quote full
circuit **ampacities** and depend on iterative sheath / bonding loss treatment.

Reproducible *today* with `thermal_cable_model`:

- Pairwise **mutual heating resistance** ``R_m`` (IEC 60287 image method), via
  :func:`thermal_cable_model.thermal_network.mutual_heating_resistance` and the
  network’s ``_Rm`` matrix in :class:`thermal_cable_model.thermal_network.CableThermalNetwork`.
- **Per-cable external resistance** ``T4`` from
  :func:`thermal_cable_model.thermal_network.external_thermal_resistance` applied
  independently to each cable (not the group-equivalent ``T4`` for touching trefoil
  as a single thermal object in IEC 60287-2-1 §4.2.4).
- **Steady-state conductor temperatures** for multi-cable installations via
  :meth:`thermal_cable_model.simulation.ThermalSimulation.run_steady_state`, using
  the package’s mutual-heating forcing at the soil node.

Not reproduced (without further physics / APIs):

- Published **rated current** for Case 1-1 / 1-2 (e.g. values tabulated on
  third-party validation sites or in TB 880 reports).
- **Touching trefoil** modelled as one equivalent external thermal resistance for
  the group vs three separate ``T4`` + image mutual terms.
- **Cross-bonding**, eddy-current, and **λ′** variants in Case 1-2.
- **2D horizontal** trefoil offset (only one horizontal axis ``x`` is available;
  the examples use a **vertical-plane** trefoil layout with different ``x`` and
  ``depth``, consistent with :file:`examples/flat_vs_trefoil.py`).

Crossing reference (CIGRE TB 640 methodology)
==============================================
:mod:`thermal_cable_model.crossing` follows the line-source / image formulation
described in *CIGRE Technical Brochure 640* (2015), §5.4.  The golden tests below
lock **numerical regression** for a fixed geometry and heat rate (not a transcribed
table from the brochure, which is not redistributed here).
"""

from __future__ import annotations

import math

import pytest

from thermal_cable_model.cable import Cable
from thermal_cable_model.crossing import CableCrossing
from thermal_cable_model.ground import ConstantGroundTemperature
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import SOIL_STANDARD, ThermalMaterial
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation
from thermal_cable_model.thermal_network import mutual_heating_resistance

# ── Parallel burial: IEC image formula for Rm (two cables, same depth) ─────

REF_RM_FLAT_0_20M_DEPTH_1M_SOIL_RHO1 = 0.3672596216100546


def _rm_closed_form_horizontal(
    soil: ThermalMaterial,
    x1: float,
    y_depth_1: float,
    x2: float,
    y_depth_2: float,
) -> float:
    """R_m = ρ/(2π)·ln(d'/d) for parallel horizontal cylinders, same reference as
    :func:`mutual_heating_resistance`.
    """
    dx = x1 - x2
    dy_real = y_depth_1 - y_depth_2
    dy_image = y_depth_1 + y_depth_2
    d_real = math.sqrt(dx**2 + dy_real**2)
    d_image = math.sqrt(dx**2 + dy_image**2)
    rho = soil.thermal_resistivity
    return rho / (2.0 * math.pi) * math.log(d_image / d_real)


class TestParallelBurialMutualResistanceReference:
    """``mutual_heating_resistance`` vs closed IEC image form; trefoil-plane Rm."""

    def test_flat_spacing_0_2_m_matches_closed_form(self):
        soil = SOIL_STANDARD
        got = mutual_heating_resistance(0.0, 1.0, 0.2, 1.0, soil)
        expected = _rm_closed_form_horizontal(soil, 0.0, 1.0, 0.2, 1.0)
        assert got == pytest.approx(expected, rel=0, abs=1e-12)
        assert got == pytest.approx(REF_RM_FLAT_0_20M_DEPTH_1M_SOIL_RHO1, rel=0, abs=1e-12)

    def test_vertical_plane_trefoil_rm_symmetry(self):
        """2-D trefoil cross-section (see ``examples/flat_vs_trefoil.py``): phases
        at two depths; symmetry implies R_{0,2} = R_{1,2} for bottom pair symmetry."""
        cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
        d = 2.0 * cable.outer_radius
        h_tri = d * math.sqrt(3) / 2.0
        depth = 1.0
        depth_top = depth - h_tri * 2.0 / 3.0
        depth_bot = depth + h_tri * 1.0 / 3.0
        tri_x = [-d / 2.0, d / 2.0, 0.0]
        tri_d = [depth_bot, depth_bot, depth_top]

        r01 = mutual_heating_resistance(
            tri_x[0], tri_d[0], tri_x[1], tri_d[1], SOIL_STANDARD
        )
        r02 = mutual_heating_resistance(
            tri_x[0], tri_d[0], tri_x[2], tri_d[2], SOIL_STANDARD
        )
        r12 = mutual_heating_resistance(
            tri_x[1], tri_d[1], tri_x[2], tri_d[2], SOIL_STANDARD
        )
        assert r02 == pytest.approx(r12, rel=0, abs=1e-9)
        assert r01 > 0.0 and r02 > 0.0
        # Bottom–bottom spacing equals d → same Rm as flat pair at same depth
        r_flat_d = mutual_heating_resistance(
            0.0, depth_bot, d, depth_bot, SOIL_STANDARD
        )
        assert r01 == pytest.approx(r_flat_d, rel=0, abs=1e-9)


class TestParallelBurialSteadyStateSymmetry:
    """Multi-cable ``ThermalSimulation`` symmetry for identical phases (flat)."""

    def test_three_flat_touching_outer_conductors_equal_and_centre_hottest(self):
        cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
        d = 2.0 * cable.outer_radius
        ground = ConstantGroundTemperature(15.0)
        load = LoadProfile.constant(400.0, 3600.0)
        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        for x in (-d, 0.0, d):
            inst.add_cable(cable, x=x, depth=1.0, load=load)
        ss = ThermalSimulation(inst).run_steady_state()
        t_left, t_mid, t_right = (
            float(ss.conductor_temps[0, 0]),
            float(ss.conductor_temps[0, 1]),
            float(ss.conductor_temps[0, 2]),
        )
        assert t_left == pytest.approx(t_right, abs=0.02)
        assert t_mid > t_left + 0.5


# ── Crossing: golden regression (TB 640 §5.4 methodology, fixed inputs) ──────

REF_CROSSING_60DEG_DEPTHS_0p9_1p2_W50_UPPER_K = 17.875320217680393
REF_CROSSING_60DEG_DEPTHS_0p9_1p2_W100_UPPER_K = 35.750640435360786


@pytest.fixture
def _tb640_style_crossing_cables():
    mv = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
    lv = Cable.three_core_xlpe_cu(150, voltage_class="LV", voltage_kv=0.6)
    return mv, lv


class TestCrossingCigreTb640StyleGolden:
    """Numerical regression for :class:`CableCrossing` at documented geometry."""

    def test_temperature_rise_upper_W50_W_per_m(self, _tb640_style_crossing_cables):
        mv, lv = _tb640_style_crossing_cables
        cx = CableCrossing(
            cable_upper=mv,
            cable_lower=lv,
            depth_upper=0.9,
            depth_lower=1.2,
            crossing_angle_deg=60.0,
            soil=SOIL_STANDARD,
        )
        d_t = cx.temperature_rise_at_upper(50.0, integration_half_length=50.0)
        assert d_t == pytest.approx(REF_CROSSING_60DEG_DEPTHS_0p9_1p2_W50_UPPER_K, abs=1e-3)

    def test_temperature_rise_scales_linearly_with_heat_rate(
        self, _tb640_style_crossing_cables
    ):
        mv, lv = _tb640_style_crossing_cables
        cx = CableCrossing(
            cable_upper=mv,
            cable_lower=lv,
            depth_upper=0.9,
            depth_lower=1.2,
            crossing_angle_deg=60.0,
            soil=SOIL_STANDARD,
        )
        d50 = cx.temperature_rise_at_upper(50.0)
        d100 = cx.temperature_rise_at_upper(100.0)
        assert d100 == pytest.approx(2.0 * d50, rel=0, abs=1e-6)
        assert d100 == pytest.approx(
            REF_CROSSING_60DEG_DEPTHS_0p9_1p2_W100_UPPER_K, abs=1e-3
        )

    def test_upper_and_lower_rise_equal_for_same_W_symmetric_depth_swap(
        self, _tb640_style_crossing_cables
    ):
        """For fixed *W*, swapping which cable is labelled upper/lower swaps depths;
        with the same angle and soil, both ``temperature_rise_at_*`` use the same
        |Δh| structure for the reciprocal problem — here depths 0.9 and 1.2 give
        identical rise for equal *W*."""
        mv, lv = _tb640_style_crossing_cables
        cx = CableCrossing(mv, lv, 0.9, 1.2, 60.0, SOIL_STANDARD)
        W = 40.0
        assert cx.temperature_rise_at_upper(W) == pytest.approx(
            cx.temperature_rise_at_lower(W), rel=0, abs=1e-6
        )
