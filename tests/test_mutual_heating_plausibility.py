"""Plausibility tests for mutual heating between parallel cables.

Verifies that the thermal network correctly models the IEC 60287
superposition-based mutual heating between cables in a group.
Each test encodes a physical invariant or monotonicity relationship —
not a regression against specific numbers.
"""

import math

import numpy as np
import pytest

from thermal_cable_model.cable import Cable
from thermal_cable_model.ground import ConstantGroundTemperature
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import SOIL_STANDARD, ThermalMaterial
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation
from thermal_cable_model.thermal_network import (
    CableThermalNetwork,
    external_thermal_resistance,
    mutual_heating_resistance,
)


# ── Shared constants ──────────────────────────────────────────────────

DEPTH = 1.0         # m
I_PHASE = 400.0     # A
DURATION = 24 * 3600  # 24 h
DT = 300             # 5-min steps
GROUND_TEMP = 15.0   # °C constant


# ── Shared fixtures ───────────────────────────────────────────────────


@pytest.fixture
def cable():
    return Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)


@pytest.fixture
def ground():
    return ConstantGroundTemperature(GROUND_TEMP)


def _run_transient(cable, xs, depths, ground, current=I_PHASE):
    """Helper: build installation and run a transient simulation."""
    load = LoadProfile.constant(current, DURATION)
    inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
    for x, d in zip(xs, depths):
        inst.add_cable(cable, x=x, depth=d, load=load)
    return ThermalSimulation(inst).run_transient(dt=DT, duration=DURATION)


def _run_steady(cable, xs, depths, ground, current=I_PHASE):
    """Helper: build installation and run a steady-state solve."""
    load = LoadProfile.constant(current, DURATION)
    inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
    for x, d in zip(xs, depths):
        inst.add_cable(cable, x=x, depth=d, load=load)
    return ThermalSimulation(inst).run_steady_state()


# ═══════════════════════════════════════════════════════════════════════
#  1. Centre cable must be hotter than outer cables in flat formation
# ═══════════════════════════════════════════════════════════════════════


class TestFlatFormationCentreHotter:
    """In a symmetric 3-cable flat formation the centre cable receives
    mutual heating from both neighbours and must be the hottest."""

    def test_transient_centre_hotter(self, cable, ground):
        d = 2.0 * cable.outer_radius
        result = _run_transient(cable, [-d, 0.0, d], [DEPTH]*3, ground)
        T_L = result.max_conductor_temp(0)
        T_C = result.max_conductor_temp(1)
        T_R = result.max_conductor_temp(2)
        assert T_C > T_L, (
            f"Centre ({T_C:.2f} °C) must exceed left ({T_L:.2f} °C)"
        )
        assert T_C > T_R, (
            f"Centre ({T_C:.2f} °C) must exceed right ({T_R:.2f} °C)"
        )

    def test_steady_state_centre_hotter(self, cable, ground):
        d = 2.0 * cable.outer_radius
        result = _run_steady(cable, [-d, 0.0, d], [DEPTH]*3, ground)
        T_L = result.conductor_temps[0, 0]
        T_C = result.conductor_temps[0, 1]
        T_R = result.conductor_temps[0, 2]
        assert T_C > T_L
        assert T_C > T_R

    def test_wider_spacing(self, cable, ground):
        """Still holds at 200 mm spacing."""
        d = 0.2
        result = _run_transient(cable, [-d, 0.0, d], [DEPTH]*3, ground)
        assert result.max_conductor_temp(1) > result.max_conductor_temp(0)
        assert result.max_conductor_temp(1) > result.max_conductor_temp(2)


# ═══════════════════════════════════════════════════════════════════════
#  2. Outer cables symmetric in flat formation
# ═══════════════════════════════════════════════════════════════════════


class TestFlatFormationSymmetry:
    """Left and right cables must have identical temperatures by symmetry."""

    def test_transient_symmetry(self, cable, ground):
        d = 2.0 * cable.outer_radius
        result = _run_transient(cable, [-d, 0.0, d], [DEPTH]*3, ground)
        T_L = result.max_conductor_temp(0)
        T_R = result.max_conductor_temp(2)
        assert T_L == pytest.approx(T_R, abs=0.01), (
            f"Left ({T_L:.2f}) and right ({T_R:.2f}) must be equal by symmetry"
        )

    def test_steady_state_symmetry(self, cable, ground):
        d = 2.0 * cable.outer_radius
        result = _run_steady(cable, [-d, 0.0, d], [DEPTH]*3, ground)
        T_L = result.conductor_temps[0, 0]
        T_R = result.conductor_temps[0, 2]
        assert T_L == pytest.approx(T_R, abs=0.01)


# ═══════════════════════════════════════════════════════════════════════
#  3. Multi-cable always hotter than single cable
# ═══════════════════════════════════════════════════════════════════════


class TestMultiCableHotterThanSingle:
    """Adding neighbouring cables must never cool any cable down."""

    def test_three_cables_hotter_than_one(self, cable, ground):
        d = 2.0 * cable.outer_radius
        result_single = _run_transient(cable, [0.0], [DEPTH], ground)
        result_three = _run_transient(
            cable, [-d, 0.0, d], [DEPTH]*3, ground
        )
        T_single = result_single.max_conductor_temp(0)
        T_centre = result_three.max_conductor_temp(1)
        T_outer = result_three.max_conductor_temp(0)
        assert T_centre > T_single, (
            f"Centre of 3-cable group ({T_centre:.2f}) must exceed "
            f"isolated cable ({T_single:.2f})"
        )
        assert T_outer > T_single, (
            f"Outer cable ({T_outer:.2f}) must also exceed "
            f"isolated cable ({T_single:.2f})"
        )

    def test_two_cables_hotter_than_one(self, cable, ground):
        d = 2.0 * cable.outer_radius
        result_single = _run_transient(cable, [0.0], [DEPTH], ground)
        result_two = _run_transient(cable, [0.0, d], [DEPTH]*2, ground)
        T_single = result_single.max_conductor_temp(0)
        T_near = result_two.max_conductor_temp(0)
        assert T_near > T_single


# ═══════════════════════════════════════════════════════════════════════
#  4. Closer spacing → higher temperature
# ═══════════════════════════════════════════════════════════════════════


class TestCloserSpacingHotter:
    """Cables closer together receive more mutual heating."""

    def test_touching_hotter_than_spaced(self, cable, ground):
        d_touch = 2.0 * cable.outer_radius
        d_wide = 0.3  # 300 mm

        res_touch = _run_transient(
            cable, [-d_touch, 0.0, d_touch], [DEPTH]*3, ground
        )
        res_wide = _run_transient(
            cable, [-d_wide, 0.0, d_wide], [DEPTH]*3, ground
        )
        T_touch_C = res_touch.max_conductor_temp(1)
        T_wide_C = res_wide.max_conductor_temp(1)
        assert T_touch_C > T_wide_C, (
            f"Touching ({T_touch_C:.2f}) must exceed 300 mm spacing "
            f"({T_wide_C:.2f})"
        )

    def test_spread_increases_with_closeness(self, cable, ground):
        """Centre–outer spread is larger when cables are closer."""
        d_touch = 2.0 * cable.outer_radius
        d_wide = 0.3

        res_touch = _run_transient(
            cable, [-d_touch, 0.0, d_touch], [DEPTH]*3, ground
        )
        res_wide = _run_transient(
            cable, [-d_wide, 0.0, d_wide], [DEPTH]*3, ground
        )
        spread_touch = (
            res_touch.max_conductor_temp(1) - res_touch.max_conductor_temp(0)
        )
        spread_wide = (
            res_wide.max_conductor_temp(1) - res_wide.max_conductor_temp(0)
        )
        assert spread_touch > spread_wide


# ═══════════════════════════════════════════════════════════════════════
#  5. Higher current → more mutual heating
# ═══════════════════════════════════════════════════════════════════════


class TestCurrentScaling:
    """Mutual heating effect grows with current (approximately I²)."""

    def test_double_current_more_mutual_rise(self, cable, ground):
        d = 2.0 * cable.outer_radius
        res_lo = _run_transient(
            cable, [-d, 0.0, d], [DEPTH]*3, ground, current=200.0
        )
        res_hi = _run_transient(
            cable, [-d, 0.0, d], [DEPTH]*3, ground, current=400.0
        )
        single_lo = _run_transient(
            cable, [0.0], [DEPTH], ground, current=200.0
        )
        single_hi = _run_transient(
            cable, [0.0], [DEPTH], ground, current=400.0
        )
        mutual_rise_lo = (
            res_lo.max_conductor_temp(1) - single_lo.max_conductor_temp(0)
        )
        mutual_rise_hi = (
            res_hi.max_conductor_temp(1) - single_hi.max_conductor_temp(0)
        )
        assert mutual_rise_hi > mutual_rise_lo, (
            f"Mutual rise at 400 A ({mutual_rise_hi:.2f}) must exceed "
            f"mutual rise at 200 A ({mutual_rise_lo:.2f})"
        )


# ═══════════════════════════════════════════════════════════════════════
#  6. Trefoil: bottom cables hotter than top
# ═══════════════════════════════════════════════════════════════════════


class TestTrefoilTemperatureOrder:
    """In a trefoil with centroid at depth, bottom cables are deeper and
    closer together → they should be at least as hot as the top cable."""

    def test_bottom_at_least_as_hot_as_top(self, cable, ground):
        d = 2.0 * cable.outer_radius
        h_tri = d * math.sqrt(3) / 2.0
        depth_top = DEPTH - h_tri * 2.0 / 3.0
        depth_bot = DEPTH + h_tri * 1.0 / 3.0

        result = _run_transient(
            cable,
            [-d/2, d/2, 0.0],
            [depth_bot, depth_bot, depth_top],
            ground,
        )
        T_BL = result.max_conductor_temp(0)
        T_BR = result.max_conductor_temp(1)
        T_top = result.max_conductor_temp(2)
        assert T_BL >= T_top - 0.05, (
            f"Bottom-left ({T_BL:.2f}) must be ≥ top ({T_top:.2f})"
        )
        assert T_BL == pytest.approx(T_BR, abs=0.01), (
            "Bottom cables must be symmetric"
        )


# ═══════════════════════════════════════════════════════════════════════
#  7. Trefoil spread is smaller than flat spread
# ═══════════════════════════════════════════════════════════════════════


class TestTrefoilSmallerSpread:
    """The trefoil is more symmetric → smaller temperature spread."""

    def test_spread_comparison(self, cable, ground):
        d = 2.0 * cable.outer_radius
        h_tri = d * math.sqrt(3) / 2.0
        depth_top = DEPTH - h_tri * 2.0 / 3.0
        depth_bot = DEPTH + h_tri * 1.0 / 3.0

        res_flat = _run_transient(
            cable, [-d, 0.0, d], [DEPTH]*3, ground
        )
        res_tri = _run_transient(
            cable,
            [-d/2, d/2, 0.0],
            [depth_bot, depth_bot, depth_top],
            ground,
        )
        flat_temps = [res_flat.max_conductor_temp(i) for i in range(3)]
        tri_temps = [res_tri.max_conductor_temp(i) for i in range(3)]

        spread_flat = max(flat_temps) - min(flat_temps)
        spread_tri = max(tri_temps) - min(tri_temps)
        assert spread_tri < spread_flat, (
            f"Trefoil spread ({spread_tri:.2f}) must be < "
            f"flat spread ({spread_flat:.2f})"
        )


# ═══════════════════════════════════════════════════════════════════════
#  8. Mutual heating resistance is positive
# ═══════════════════════════════════════════════════════════════════════


class TestMutualHeatingResistancePositive:
    """Rm must be > 0 for distinct cables buried in the ground."""

    @pytest.mark.parametrize("dx, dy", [
        (0.034, 0.0),
        (0.2, 0.0),
        (0.0, 0.3),
        (0.1, 0.1),
        (1.0, 0.0),
    ])
    def test_positive(self, dx, dy):
        Rm = mutual_heating_resistance(0.0, DEPTH, dx, DEPTH + dy, SOIL_STANDARD)
        assert Rm > 0, f"Rm must be positive for dx={dx}, dy={dy}"

    def test_symmetric(self):
        """Rm(i→j) == Rm(j→i) for cables at the same depth."""
        Rm_ij = mutual_heating_resistance(0.0, DEPTH, 0.2, DEPTH, SOIL_STANDARD)
        Rm_ji = mutual_heating_resistance(0.2, DEPTH, 0.0, DEPTH, SOIL_STANDARD)
        assert Rm_ij == pytest.approx(Rm_ji, rel=1e-12)

    def test_closer_cables_larger_Rm(self):
        """Closer cables have larger Rm (more mutual temperature rise)."""
        Rm_close = mutual_heating_resistance(0.0, DEPTH, 0.05, DEPTH, SOIL_STANDARD)
        Rm_far = mutual_heating_resistance(0.0, DEPTH, 0.5, DEPTH, SOIL_STANDARD)
        assert Rm_close > Rm_far


# ═══════════════════════════════════════════════════════════════════════
#  9. Steady-state matches IEC 60287 analytical formula
# ═══════════════════════════════════════════════════════════════════════


class TestSteadyStateMatchesIEC:
    """For a single cable, T_conductor = T_ambient + W × (T1+T2+T3+T4).
    For two cables, the mutual rise ΔT = W_j × Rm_ij should appear."""

    def test_single_cable_baseline(self, cable, ground):
        """Verify the single-cable temperature is consistent with the
        thermal resistance chain."""
        result = _run_steady(cable, [0.0], [DEPTH], ground)
        T_cond = result.conductor_temps[0, 0]
        assert T_cond > GROUND_TEMP, "Conductor must be above ambient"

    def test_mutual_rise_magnitude(self, cable, ground):
        """The temperature difference between a cable in a group and the
        same cable in isolation should approximately equal the IEC mutual
        temperature rise: ΔT ≈ Σ W_j × Rm_ij."""
        d = 0.2

        res_single = _run_steady(cable, [0.0], [DEPTH], ground)
        res_pair = _run_steady(cable, [0.0, d], [DEPTH]*2, ground)

        T_single = res_single.conductor_temps[0, 0]
        T_pair_0 = res_pair.conductor_temps[0, 0]
        T_pair_1 = res_pair.conductor_temps[0, 1]

        Rm = mutual_heating_resistance(0.0, DEPTH, d, DEPTH, SOIL_STANDARD)
        W_approx = cable.total_heat_per_length(I_PHASE, T_single)
        expected_rise = W_approx * Rm

        actual_rise_0 = T_pair_0 - T_single
        actual_rise_1 = T_pair_1 - T_single

        assert actual_rise_0 == pytest.approx(expected_rise, rel=0.15), (
            f"Cable 0 mutual rise ({actual_rise_0:.2f}) should ≈ "
            f"W×Rm = {expected_rise:.2f}"
        )
        assert actual_rise_1 == pytest.approx(expected_rise, rel=0.15), (
            f"Cable 1 mutual rise ({actual_rise_1:.2f}) should ≈ "
            f"W×Rm = {expected_rise:.2f}"
        )
        assert T_pair_0 == pytest.approx(T_pair_1, abs=0.01), (
            "Symmetric pair must have equal temperatures"
        )


# ═══════════════════════════════════════════════════════════════════════
#  10. Deeper burial → higher temperature (weaker cooling to surface)
# ═══════════════════════════════════════════════════════════════════════


class TestDeeperBurialHotter:
    """Deeper burial increases T4, so the cable runs hotter."""

    def test_single_cable(self, cable, ground):
        res_shallow = _run_transient(cable, [0.0], [0.7], ground)
        res_deep = _run_transient(cable, [0.0], [1.5], ground)
        assert res_deep.max_conductor_temp(0) > res_shallow.max_conductor_temp(0)

    def test_group_of_three(self, cable, ground):
        d = 0.2
        res_shallow = _run_transient(
            cable, [-d, 0.0, d], [0.7]*3, ground
        )
        res_deep = _run_transient(
            cable, [-d, 0.0, d], [1.5]*3, ground
        )
        assert res_deep.max_conductor_temp(1) > res_shallow.max_conductor_temp(1)


# ═══════════════════════════════════════════════════════════════════════
#  11. Soil resistivity effect on mutual heating
# ═══════════════════════════════════════════════════════════════════════


class TestSoilResistivityEffect:
    """Higher soil thermal resistivity → more mutual heating."""

    def test_higher_resistivity_hotter(self, cable):
        soil_lo = ThermalMaterial(
            name="low-rho", thermal_conductivity=1.5,
            thermal_resistivity=1.0 / 1.5,
            volumetric_heat_capacity=SOIL_STANDARD.volumetric_heat_capacity,
        )
        soil_hi = ThermalMaterial(
            name="high-rho", thermal_conductivity=0.5,
            thermal_resistivity=1.0 / 0.5,
            volumetric_heat_capacity=SOIL_STANDARD.volumetric_heat_capacity,
        )
        ground = ConstantGroundTemperature(GROUND_TEMP)
        d = 0.2
        load = LoadProfile.constant(I_PHASE, DURATION)

        def run(soil):
            inst = CableInstallation(soil=soil, ground_temp_model=ground)
            for x in [-d, 0.0, d]:
                inst.add_cable(cable, x=x, depth=DEPTH, load=load)
            return ThermalSimulation(inst).run_transient(dt=DT, duration=DURATION)

        res_lo = run(soil_lo)
        res_hi = run(soil_hi)
        assert res_hi.max_conductor_temp(1) > res_lo.max_conductor_temp(1)


# ═══════════════════════════════════════════════════════════════════════
#  12. Mutual heating matrix stored correctly in network
# ═══════════════════════════════════════════════════════════════════════


class TestNetworkRmMatrix:
    """The Rm matrix stored in CableThermalNetwork must match the
    standalone mutual_heating_resistance function."""

    def test_matrix_matches_function(self, cable):
        xs = [-0.1, 0.0, 0.1]
        ds = [DEPTH, DEPTH, DEPTH]
        net = CableThermalNetwork(
            cables=[cable]*3, positions_x=xs, depths=ds, soil=SOIL_STANDARD,
        )
        for i in range(3):
            for j in range(3):
                expected = mutual_heating_resistance(
                    xs[i], ds[i], xs[j], ds[j], SOIL_STANDARD,
                )
                assert net._Rm[i, j] == pytest.approx(expected, rel=1e-12)

    def test_diagonal_is_zero(self, cable):
        xs = [0.0, 0.2]
        net = CableThermalNetwork(
            cables=[cable]*2, positions_x=xs,
            depths=[DEPTH]*2, soil=SOIL_STANDARD,
        )
        assert net._Rm[0, 0] == 0.0
        assert net._Rm[1, 1] == 0.0


# ═══════════════════════════════════════════════════════════════════════
#  13. No mutual heating for a single cable
# ═══════════════════════════════════════════════════════════════════════


class TestSingleCableNoMutual:
    """A single cable must produce the same result regardless of the
    mutual heating implementation (Rm matrix is 1×1 with zero diagonal)."""

    def test_rm_zero(self, cable):
        net = CableThermalNetwork(
            cables=[cable], positions_x=[0.0],
            depths=[DEPTH], soil=SOIL_STANDARD,
        )
        assert net._Rm[0, 0] == 0.0

    def test_temperature_unchanged(self, cable, ground):
        result = _run_transient(cable, [0.0], [DEPTH], ground)
        T = result.max_conductor_temp(0)
        assert T > GROUND_TEMP
        assert T < cable.max_conductor_temp


# ═══════════════════════════════════════════════════════════════════════
#  14. Flat formation spread is reasonable order of magnitude
# ═══════════════════════════════════════════════════════════════════════


class TestFlatSpreadOrderOfMagnitude:
    """The centre–outer spread in a flat touching formation should be
    on the order of 0.5–5 °C for typical MV cables at rated current."""

    def test_spread_range(self, cable, ground):
        d = 2.0 * cable.outer_radius
        result = _run_transient(cable, [-d, 0.0, d], [DEPTH]*3, ground)
        spread = (
            result.max_conductor_temp(1) - result.max_conductor_temp(0)
        )
        assert 0.3 < spread < 10.0, (
            f"Spread {spread:.2f} °C outside plausible range"
        )
