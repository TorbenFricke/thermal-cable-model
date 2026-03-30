"""Tests that transient simulations converge to steady-state results.

For each scenario (isolated cable, parallel cables, crossing cables),
a long transient with constant load and constant ground temperature is
run.  The final transient temperatures must match the steady-state
solution within a specified tolerance.
"""

import math

import numpy as np
import pytest

from thermal_cable_model.cable import Cable
from thermal_cable_model.crossing import CableCrossing
from thermal_cable_model.ground import ConstantGroundTemperature
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import SOIL_STANDARD
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation


# ── Shared constants ──────────────────────────────────────────────────

GROUND_TEMP = 15.0
DT = 300  # 5-min time step [s]


# ── Shared fixtures ───────────────────────────────────────────────────


@pytest.fixture
def ground():
    return ConstantGroundTemperature(GROUND_TEMP)


@pytest.fixture
def mv_cable():
    return Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)


@pytest.fixture
def lv_cable():
    return Cable.three_core_xlpe_cu(150, voltage_class="LV", voltage_kv=0.6)


# ═══════════════════════════════════════════════════════════════════════
#  1. Isolated cable: transient → steady-state
# ═══════════════════════════════════════════════════════════════════════


class TestIsolatedCableConvergence:
    """A single cable with constant load must reach the steady-state
    temperature when the transient is run long enough."""

    CURRENT = 400.0
    DURATION = 72 * 3600  # 72 hours
    TOL_K = 0.05  # 50 mK tolerance

    def test_conductor_temperature(self, mv_cable, ground):
        duration = self.DURATION
        load = LoadProfile.constant(self.CURRENT, duration)

        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(mv_cable, x=0.0, depth=1.2, load=load)
        sim = ThermalSimulation(inst)

        ss = sim.run_steady_state()
        tr = sim.run_transient(dt=DT, duration=duration)

        ss_cond = ss.conductor_temps[0, 0]
        tr_cond = tr.conductor_temps[-1, 0]
        assert abs(tr_cond - ss_cond) < self.TOL_K, (
            f"Transient conductor temp {tr_cond:.4f} °C != "
            f"steady-state {ss_cond:.4f} °C (Δ={tr_cond - ss_cond:.4f})"
        )

    def test_all_node_temperatures(self, mv_cable, ground):
        """Every node (insulation, sheath, surface, soil) must also
        converge to the steady-state."""
        duration = self.DURATION
        load = LoadProfile.constant(self.CURRENT, duration)

        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(mv_cable, x=0.0, depth=1.2, load=load)
        sim = ThermalSimulation(inst)

        ss = sim.run_steady_state()
        tr = sim.run_transient(dt=DT, duration=duration)

        for label, ss_arr, tr_arr in [
            ("insulation", ss.insulation_temps, tr.insulation_temps),
            ("sheath", ss.sheath_temps, tr.sheath_temps),
            ("armour", ss.armour_temps, tr.armour_temps),
            ("surface", ss.surface_temps, tr.surface_temps),
            ("soil", ss.soil_temps, tr.soil_temps),
        ]:
            diff = abs(tr_arr[-1, 0] - ss_arr[0, 0])
            assert diff < self.TOL_K, (
                f"{label}: transient {tr_arr[-1, 0]:.4f} != "
                f"steady-state {ss_arr[0, 0]:.4f} (Δ={diff:.4f})"
            )

    def test_three_core_cable(self, lv_cable, ground):
        """Multi-core cable (n > 1) must also converge."""
        duration = self.DURATION
        load = LoadProfile.constant(280.0, duration)

        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(lv_cable, x=0.0, depth=1.0, load=load)
        sim = ThermalSimulation(inst)

        ss = sim.run_steady_state()
        tr = sim.run_transient(dt=DT, duration=duration)

        ss_cond = ss.conductor_temps[0, 0]
        tr_cond = tr.conductor_temps[-1, 0]
        assert abs(tr_cond - ss_cond) < self.TOL_K, (
            f"3-core transient {tr_cond:.4f} != "
            f"steady-state {ss_cond:.4f} (Δ={tr_cond - ss_cond:.4f})"
        )


# ═══════════════════════════════════════════════════════════════════════
#  2. Parallel cables: transient → steady-state
# ═══════════════════════════════════════════════════════════════════════


class TestParallelCableConvergence:
    """Three cables in flat formation with constant load.  Each cable's
    transient temperature must converge to its steady-state value,
    including the mutual heating effect."""

    CURRENT = 400.0
    DURATION = 72 * 3600
    TOL_K = 0.1  # slightly larger tolerance for coupled system

    def test_flat_formation(self, mv_cable, ground):
        duration = self.DURATION
        load = LoadProfile.constant(self.CURRENT, duration)

        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(mv_cable, x=-0.3, depth=1.0, load=load)
        inst.add_cable(mv_cable, x=0.0, depth=1.0, load=load)
        inst.add_cable(mv_cable, x=0.3, depth=1.0, load=load)
        sim = ThermalSimulation(inst)

        ss = sim.run_steady_state()
        tr = sim.run_transient(dt=DT, duration=duration)

        for k in range(3):
            ss_cond = ss.conductor_temps[0, k]
            tr_cond = tr.conductor_temps[-1, k]
            assert abs(tr_cond - ss_cond) < self.TOL_K, (
                f"Cable {k}: transient {tr_cond:.4f} != "
                f"steady-state {ss_cond:.4f} (Δ={tr_cond - ss_cond:.4f})"
            )

    def test_centre_hotter_than_outer(self, mv_cable, ground):
        """Sanity check: the centre cable in flat formation is always
        the hottest, both in steady-state and at the end of transient."""
        duration = self.DURATION
        load = LoadProfile.constant(self.CURRENT, duration)

        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(mv_cable, x=-0.3, depth=1.0, load=load)
        inst.add_cable(mv_cable, x=0.0, depth=1.0, load=load)
        inst.add_cable(mv_cable, x=0.3, depth=1.0, load=load)
        sim = ThermalSimulation(inst)

        ss = sim.run_steady_state()
        tr = sim.run_transient(dt=DT, duration=duration)

        assert ss.conductor_temps[0, 1] > ss.conductor_temps[0, 0]
        assert tr.conductor_temps[-1, 1] > tr.conductor_temps[-1, 0]

    def test_symmetric_outer_cables(self, mv_cable, ground):
        """With identical loads the two outer cables must reach the
        same temperature (symmetry)."""
        duration = self.DURATION
        load = LoadProfile.constant(self.CURRENT, duration)

        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(mv_cable, x=-0.3, depth=1.0, load=load)
        inst.add_cable(mv_cable, x=0.0, depth=1.0, load=load)
        inst.add_cable(mv_cable, x=0.3, depth=1.0, load=load)
        sim = ThermalSimulation(inst)

        ss = sim.run_steady_state()
        tr = sim.run_transient(dt=DT, duration=duration)

        assert abs(ss.conductor_temps[0, 0] - ss.conductor_temps[0, 2]) < 0.01
        assert abs(tr.conductor_temps[-1, 0] - tr.conductor_temps[-1, 2]) < 0.01


# ═══════════════════════════════════════════════════════════════════════
#  3. Crossing cables: transient → steady-state
# ═══════════════════════════════════════════════════════════════════════


class TestCrossingConvergence:
    """The transient crossing ΔT must converge to the steady-state ΔT
    as time increases.  Tests are performed on the crossing module
    directly (much faster than full simulation)."""

    @pytest.fixture
    def crossing(self, mv_cable, lv_cable):
        return CableCrossing(
            cable_upper=mv_cable,
            cable_lower=lv_cable,
            depth_upper=0.9,
            depth_lower=1.2,
            crossing_angle_deg=60,
            soil=SOIL_STANDARD,
        )

    def test_upper_cable_convergence(self, crossing, lv_cable):
        """Transient ΔT at upper cable must approach the steady-state
        to within 1% at a sufficiently large time."""
        W_lv = lv_cable.total_heat_per_length(280.0, 50.0)
        ss = crossing.temperature_rise_at_upper(W_lv)
        tr = crossing.transient_temperature_rise_at_upper(W_lv, 50_000 * 3600)
        assert abs(tr - ss) / ss < 0.01, (
            f"transient {tr:.4f} vs steady-state {ss:.4f} "
            f"(error = {abs(tr-ss)/ss*100:.2f}%)"
        )

    def test_lower_cable_convergence(self, crossing, mv_cable):
        W_mv = mv_cable.total_heat_per_length(350.0, 70.0)
        ss = crossing.temperature_rise_at_lower(W_mv)
        tr = crossing.transient_temperature_rise_at_lower(W_mv, 50_000 * 3600)
        assert abs(tr - ss) / ss < 0.01, (
            f"transient {tr:.4f} vs steady-state {ss:.4f} "
            f"(error = {abs(tr-ss)/ss*100:.2f}%)"
        )

    def test_transient_always_below_steady_state(self, crossing, lv_cable):
        """At every time, the transient ΔT must be strictly less than
        the steady-state value."""
        W = lv_cable.total_heat_per_length(280.0, 50.0)
        ss = crossing.temperature_rise_at_upper(W)
        for t_h in [1, 10, 100, 1_000, 10_000]:
            tr = crossing.transient_temperature_rise_at_upper(W, t_h * 3600)
            assert tr <= ss * 1.001, (
                f"t={t_h}h: transient {tr:.4f} > steady-state {ss:.4f}"
            )

    def test_monotonic_approach(self, crossing, lv_cable):
        """The transient ΔT must increase monotonically with time."""
        W = lv_cable.total_heat_per_length(280.0, 50.0)
        prev = 0.0
        for t_h in [1, 4, 12, 48, 200, 1_000, 10_000]:
            dT = crossing.transient_temperature_rise_at_upper(W, t_h * 3600)
            assert dT >= prev, (
                f"Non-monotonic: dT({t_h}h)={dT:.4f} < previous {prev:.4f}"
            )
            prev = dT

    @pytest.mark.parametrize("angle_deg", [30, 45, 60, 90])
    def test_convergence_all_angles(self, mv_cable, lv_cable, angle_deg):
        """Convergence must hold for every crossing angle."""
        crossing = CableCrossing(
            cable_upper=mv_cable,
            cable_lower=lv_cable,
            depth_upper=0.9,
            depth_lower=1.2,
            crossing_angle_deg=angle_deg,
            soil=SOIL_STANDARD,
        )
        W = lv_cable.total_heat_per_length(280.0, 50.0)
        ss = crossing.temperature_rise_at_upper(W)
        tr = crossing.transient_temperature_rise_at_upper(W, 50_000 * 3600)
        assert abs(tr - ss) / ss < 0.02, (
            f"angle={angle_deg}°: transient {tr:.4f} vs "
            f"steady-state {ss:.4f} (error = {abs(tr-ss)/ss*100:.2f}%)"
        )
