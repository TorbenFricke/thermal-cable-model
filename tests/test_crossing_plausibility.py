"""Plausibility tests for the cable crossing thermal model.

Each test encodes a physical invariant or monotonicity relationship that
must hold for any valid crossing configuration.  These are *not* regression
tests against specific numbers — they verify that the model behaves in a
physically consistent way.
"""

import math

import numpy as np
import pytest

from thermal_cable_model.cable import Cable
from thermal_cable_model.crossing import CableCrossing, crossing_derating_factor
from thermal_cable_model.materials import SOIL_STANDARD, ThermalMaterial
from thermal_cable_model.thermal_network import mutual_heating_resistance


# ── Shared fixtures ───────────────────────────────────────────────────


@pytest.fixture
def mv_cable():
    return Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)


@pytest.fixture
def lv_cable():
    return Cable.three_core_xlpe_cu(150, voltage_class="LV", voltage_kv=0.6)


@pytest.fixture
def standard_crossing(mv_cable, lv_cable):
    return CableCrossing(
        cable_upper=mv_cable,
        cable_lower=lv_cable,
        depth_upper=0.9,
        depth_lower=1.2,
        crossing_angle_deg=60,
        soil=SOIL_STANDARD,
    )


DEPTH_UPPER = 0.9
DEPTH_LOWER = 1.2
I_MV = 350.0
I_LV = 280.0


# ═══════════════════════════════════════════════════════════════════════
#  1. Crossing ΔT is always positive
# ═══════════════════════════════════════════════════════════════════════


class TestPositiveTemperatureRise:
    def test_steady_state_upper(self, standard_crossing, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        assert standard_crossing.temperature_rise_at_upper(W) > 0

    def test_steady_state_lower(self, standard_crossing, mv_cable):
        W = mv_cable.total_heat_per_length(I_MV, 70.0)
        assert standard_crossing.temperature_rise_at_lower(W) > 0

    @pytest.mark.parametrize("time_s", [60, 600, 3600, 86400, 7 * 86400])
    def test_transient_upper(self, standard_crossing, lv_cable, time_s):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        assert standard_crossing.transient_temperature_rise_at_upper(W, time_s) > 0

    @pytest.mark.parametrize("time_s", [60, 600, 3600, 86400, 7 * 86400])
    def test_transient_lower(self, standard_crossing, mv_cable, time_s):
        W = mv_cable.total_heat_per_length(I_MV, 70.0)
        assert standard_crossing.transient_temperature_rise_at_lower(W, time_s) > 0


# ═══════════════════════════════════════════════════════════════════════
#  2. Shallower crossing angle → larger ΔT
# ═══════════════════════════════════════════════════════════════════════


class TestAngleMonotonicity:
    ANGLES = [15, 30, 45, 60, 75, 90]

    def test_steady_state_dT_decreases_with_angle(self, mv_cable, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dTs = []
        for angle in self.ANGLES:
            cx = CableCrossing(
                cable_upper=mv_cable, cable_lower=lv_cable,
                depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
                crossing_angle_deg=angle, soil=SOIL_STANDARD,
            )
            dTs.append(cx.temperature_rise_at_upper(W))

        for i in range(1, len(dTs)):
            assert dTs[i] < dTs[i - 1], (
                f"ΔT at {self.ANGLES[i]}° ({dTs[i]:.4f}) should be < "
                f"ΔT at {self.ANGLES[i-1]}° ({dTs[i-1]:.4f})"
            )

    def test_transient_dT_decreases_with_angle(self, mv_cable, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        t = 24 * 3600.0
        dTs = []
        for angle in self.ANGLES:
            cx = CableCrossing(
                cable_upper=mv_cable, cable_lower=lv_cable,
                depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
                crossing_angle_deg=angle, soil=SOIL_STANDARD,
            )
            dTs.append(cx.transient_temperature_rise_at_upper(W, t))

        for i in range(1, len(dTs)):
            assert dTs[i] < dTs[i - 1]


# ═══════════════════════════════════════════════════════════════════════
#  3. Transient ΔT monotonically approaches but never exceeds
#     steady-state
# ═══════════════════════════════════════════════════════════════════════


class TestTransientConvergence:
    TIMES_S = [60, 600, 3600, 6 * 3600, 86400, 7 * 86400, 30 * 86400, 365 * 86400]

    def test_monotonically_increasing(self, standard_crossing, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dTs = [
            standard_crossing.transient_temperature_rise_at_upper(W, t)
            for t in self.TIMES_S
        ]
        for i in range(1, len(dTs)):
            assert dTs[i] >= dTs[i - 1] - 1e-6, (
                f"ΔT at t={self.TIMES_S[i]}s ({dTs[i]:.6f}) decreased from "
                f"t={self.TIMES_S[i-1]}s ({dTs[i-1]:.6f})"
            )

    def test_bounded_by_steady_state(self, standard_crossing, lv_cable):
        """Physically, a step-on transient must stay ≤ the steady-state."""
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dT_ss = standard_crossing.temperature_rise_at_upper(W)
        for t in self.TIMES_S:
            dT_t = standard_crossing.transient_temperature_rise_at_upper(W, t)
            assert dT_t <= dT_ss, (
                f"Transient ΔT at t={t}s ({dT_t:.4f}) exceeds "
                f"steady-state ({dT_ss:.4f})"
            )

    def test_converges_to_steady_state(self, standard_crossing, lv_cable):
        """After 1 year the transient should be within 5 % of steady-state."""
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dT_ss = standard_crossing.temperature_rise_at_upper(W)
        dT_1y = standard_crossing.transient_temperature_rise_at_upper(W, 365 * 86400)
        assert dT_1y >= 0.95 * dT_ss


# ═══════════════════════════════════════════════════════════════════════
#  4. Smaller vertical separation → larger ΔT
# ═══════════════════════════════════════════════════════════════════════


class TestVerticalSeparation:
    # Depths chosen so h_sum = 2.1 m is constant, isolating the dh effect.
    DEPTH_PAIRS = [
        (0.75, 1.35),  # sep = 0.60 m
        (0.85, 1.25),  # sep = 0.40 m
        (0.90, 1.20),  # sep = 0.30 m
        (1.00, 1.10),  # sep = 0.10 m
    ]

    def test_closer_cables_higher_dT(self, mv_cable, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dTs = []
        for d_up, d_low in self.DEPTH_PAIRS:
            cx = CableCrossing(
                cable_upper=mv_cable, cable_lower=lv_cable,
                depth_upper=d_up, depth_lower=d_low,
                crossing_angle_deg=60, soil=SOIL_STANDARD,
            )
            dTs.append(cx.temperature_rise_at_upper(W))

        seps = [abs(d_low - d_up) for d_up, d_low in self.DEPTH_PAIRS]
        for i in range(1, len(dTs)):
            assert dTs[i] > dTs[i - 1], (
                f"ΔT at sep={seps[i]:.2f} m ({dTs[i]:.4f}) should be > "
                f"ΔT at sep={seps[i-1]:.2f} m ({dTs[i-1]:.4f})"
            )


# ═══════════════════════════════════════════════════════════════════════
#  5. ΔT scales linearly with heat rate
# ═══════════════════════════════════════════════════════════════════════


class TestLinearScaling:
    def test_steady_state_linearity(self, standard_crossing, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dT_1x = standard_crossing.temperature_rise_at_upper(W)
        dT_2x = standard_crossing.temperature_rise_at_upper(2.0 * W)
        dT_half = standard_crossing.temperature_rise_at_upper(0.5 * W)

        assert dT_2x == pytest.approx(2.0 * dT_1x, abs=1e-10)
        assert dT_half == pytest.approx(0.5 * dT_1x, abs=1e-10)

    def test_transient_linearity(self, standard_crossing, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        t = 86400.0
        dT_1x = standard_crossing.transient_temperature_rise_at_upper(W, t)
        dT_2x = standard_crossing.transient_temperature_rise_at_upper(2.0 * W, t)

        assert dT_2x == pytest.approx(2.0 * dT_1x, abs=1e-8)


# ═══════════════════════════════════════════════════════════════════════
#  6. Derating factor ∈ (0, 1]
# ═══════════════════════════════════════════════════════════════════════


class TestDeratingFactor:
    def test_range(self, mv_cable):
        for dT in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]:
            df = crossing_derating_factor(mv_cable, DEPTH_UPPER, dT, SOIL_STANDARD)
            assert 0.0 <= df <= 1.0, f"df={df:.4f} out of [0,1] for ΔT={dT}"

    def test_zero_crossing_gives_unity(self, mv_cable):
        df = crossing_derating_factor(mv_cable, DEPTH_UPPER, 0.0, SOIL_STANDARD)
        assert df == pytest.approx(1.0)

    def test_monotonically_decreasing(self, mv_cable):
        dTs = [0.0, 1.0, 3.0, 5.0, 10.0, 20.0]
        dfs = [
            crossing_derating_factor(mv_cable, DEPTH_UPPER, dT, SOIL_STANDARD)
            for dT in dTs
        ]
        for i in range(1, len(dfs)):
            assert dfs[i] <= dfs[i - 1], (
                f"df at ΔT={dTs[i]} ({dfs[i]:.4f}) should be <= "
                f"df at ΔT={dTs[i-1]} ({dfs[i-1]:.4f})"
            )


# ═══════════════════════════════════════════════════════════════════════
#  7. Crossing ΔT ≈ parallel ΔT / sin(α)
# ═══════════════════════════════════════════════════════════════════════


class TestCrossingVsParallelRelationship:
    """The analytical relationship (for L → ∞) is:
        ΔT_crossing = ΔT_parallel / sin(α)
    At 90° they are equal; at smaller angles the crossing ΔT is LARGER
    because the source cable lingers near the crossing point."""

    @pytest.mark.parametrize("angle_deg", [30, 45, 60, 75])
    def test_crossing_exceeds_parallel_below_90(self, mv_cable, lv_cable, angle_deg):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)

        Rm = mutual_heating_resistance(
            0.0, DEPTH_UPPER, 0.0, DEPTH_LOWER, SOIL_STANDARD,
        )
        dT_parallel = W * Rm

        cx = CableCrossing(
            cable_upper=mv_cable, cable_lower=lv_cable,
            depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
            crossing_angle_deg=angle_deg, soil=SOIL_STANDARD,
        )
        dT_crossing = cx.temperature_rise_at_upper(W)

        assert dT_crossing > dT_parallel, (
            f"Crossing ΔT at {angle_deg}° ({dT_crossing:.4f}) should be > "
            f"parallel ΔT ({dT_parallel:.4f})"
        )

    def test_crossing_at_90_approx_equals_parallel(self, mv_cable, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)

        Rm = mutual_heating_resistance(
            0.0, DEPTH_UPPER, 0.0, DEPTH_LOWER, SOIL_STANDARD,
        )
        dT_parallel = W * Rm

        cx = CableCrossing(
            cable_upper=mv_cable, cable_lower=lv_cable,
            depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
            crossing_angle_deg=90, soil=SOIL_STANDARD,
        )
        dT_crossing = cx.temperature_rise_at_upper(W)

        assert dT_crossing == pytest.approx(dT_parallel, rel=0.02)

    @pytest.mark.parametrize("angle_deg", [30, 45, 60, 90])
    def test_scales_as_inverse_sin(self, mv_cable, lv_cable, angle_deg):
        """ΔT_crossing should be approximately ΔT_parallel / sin(α)."""
        W = lv_cable.total_heat_per_length(I_LV, 55.0)

        Rm = mutual_heating_resistance(
            0.0, DEPTH_UPPER, 0.0, DEPTH_LOWER, SOIL_STANDARD,
        )
        dT_parallel = W * Rm
        expected = dT_parallel / math.sin(math.radians(angle_deg))

        cx = CableCrossing(
            cable_upper=mv_cable, cable_lower=lv_cable,
            depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
            crossing_angle_deg=angle_deg, soil=SOIL_STANDARD,
        )
        dT_crossing = cx.temperature_rise_at_upper(W)

        assert dT_crossing == pytest.approx(expected, rel=0.03)


# ═══════════════════════════════════════════════════════════════════════
#  8. Higher soil thermal conductivity → lower ΔT
# ═══════════════════════════════════════════════════════════════════════


class TestSoilConductivityEffect:
    SOILS = [
        ThermalMaterial.from_conductivity("Low λ", 0.5, 1.5e6),
        SOIL_STANDARD,  # λ = 1.0
        ThermalMaterial.from_conductivity("High λ", 2.0, 2.5e6),
    ]

    def test_higher_conductivity_lower_dT(self, mv_cable, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dTs = []
        for soil in self.SOILS:
            cx = CableCrossing(
                cable_upper=mv_cable, cable_lower=lv_cable,
                depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
                crossing_angle_deg=60, soil=soil,
            )
            dTs.append(cx.temperature_rise_at_upper(W))

        for i in range(1, len(dTs)):
            assert dTs[i] < dTs[i - 1], (
                f"ΔT in λ={self.SOILS[i].thermal_conductivity} ({dTs[i]:.4f}) "
                f"should be < ΔT in λ={self.SOILS[i-1].thermal_conductivity} "
                f"({dTs[i-1]:.4f})"
            )


# ═══════════════════════════════════════════════════════════════════════
#  9. Superposition: isolated + crossing ΔT > isolated  (end-to-end)
# ═══════════════════════════════════════════════════════════════════════


class TestSuperpositionEndToEnd:
    def test_crossing_always_raises_temperature(self, mv_cable, lv_cable):
        """Full workflow: isolated sims + analytical crossing superposition."""
        from thermal_cable_model.ground import KasudaModel
        from thermal_cable_model.loads import LoadProfile
        from thermal_cable_model.simulation import CableInstallation, ThermalSimulation

        crossing = CableCrossing(
            cable_upper=mv_cable, cable_lower=lv_cable,
            depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
            crossing_angle_deg=60, soil=SOIL_STANDARD,
        )

        ground = KasudaModel(mean_surface_temp=12.0, annual_amplitude=10.0)
        duration_s = 6 * 3600
        load_mv = LoadProfile.constant(I_MV, duration_s)
        load_lv = LoadProfile.constant(I_LV, duration_s)

        inst_mv = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst_mv.add_cable(mv_cable, x=0.0, depth=DEPTH_UPPER, load=load_mv)
        result_mv = ThermalSimulation(inst_mv).run_transient(dt=600, duration=duration_s)

        inst_lv = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst_lv.add_cable(lv_cable, x=0.0, depth=DEPTH_LOWER, load=load_lv)
        result_lv = ThermalSimulation(inst_lv).run_transient(dt=600, duration=duration_s)

        for i, t in enumerate(result_mv.times):
            if t <= 0:
                continue
            W_lv_t = lv_cable.total_heat_per_length(I_LV, result_lv.conductor_temps[i, 0])
            W_mv_t = mv_cable.total_heat_per_length(I_MV, result_mv.conductor_temps[i, 0])

            dT_mv = crossing.transient_temperature_rise_at_upper(W_lv_t, t)
            dT_lv = crossing.transient_temperature_rise_at_lower(W_mv_t, t)

            assert dT_mv > 0, f"ΔT at MV should be > 0 at t={t:.0f} s, got {dT_mv}"
            assert dT_lv > 0, f"ΔT at LV should be > 0 at t={t:.0f} s, got {dT_lv}"

        max_mv_isolated = result_mv.max_conductor_temp(0)
        max_lv_isolated = result_lv.max_conductor_temp(0)

        mv_with_crossing = result_mv.conductor_temps[:, 0].copy()
        lv_with_crossing = result_lv.conductor_temps[:, 0].copy()
        for i, t in enumerate(result_mv.times):
            t_eff = max(t, 1.0)
            W_lv_t = lv_cable.total_heat_per_length(I_LV, result_lv.conductor_temps[i, 0])
            W_mv_t = mv_cable.total_heat_per_length(I_MV, result_mv.conductor_temps[i, 0])
            mv_with_crossing[i] += crossing.transient_temperature_rise_at_upper(W_lv_t, t_eff)
            lv_with_crossing[i] += crossing.transient_temperature_rise_at_lower(W_mv_t, t_eff)

        assert float(np.max(mv_with_crossing)) > max_mv_isolated
        assert float(np.max(lv_with_crossing)) > max_lv_isolated


# ═══════════════════════════════════════════════════════════════════════
#  10. Reciprocity: same heat rate → same mutual ΔT regardless of role
# ═══════════════════════════════════════════════════════════════════════


class TestReciprocity:
    def test_symmetric_geometry_gives_equal_mutual_rise(self):
        """For the same W the integral depends only on |dh| and h_sum,
        both of which are symmetric when swapping source/target roles."""
        cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
        W = cable.total_heat_per_length(400, 70.0)

        cx = CableCrossing(
            cable_upper=cable, cable_lower=cable,
            depth_upper=0.8, depth_lower=1.4,
            crossing_angle_deg=45, soil=SOIL_STANDARD,
        )

        dT_at_upper = cx.temperature_rise_at_upper(W)
        dT_at_lower = cx.temperature_rise_at_lower(W)

        assert dT_at_upper == pytest.approx(dT_at_lower, rel=1e-10), (
            f"ΔT at upper ({dT_at_upper:.6f}) should equal "
            f"ΔT at lower ({dT_at_lower:.6f}) for same W"
        )

    def test_transient_reciprocity(self):
        cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
        W = cable.total_heat_per_length(400, 70.0)

        cx = CableCrossing(
            cable_upper=cable, cable_lower=cable,
            depth_upper=0.7, depth_lower=1.3,
            crossing_angle_deg=60, soil=SOIL_STANDARD,
        )

        for t in [3600, 86400, 7 * 86400]:
            dT_up = cx.transient_temperature_rise_at_upper(W, t)
            dT_lo = cx.transient_temperature_rise_at_lower(W, t)
            assert dT_up == pytest.approx(dT_lo, rel=1e-6), (
                f"Transient reciprocity failed at t={t}s: "
                f"upper={dT_up:.6f}, lower={dT_lo:.6f}"
            )


# ═══════════════════════════════════════════════════════════════════════
#  11. Transient strictly below steady-state for any finite time
# ═══════════════════════════════════════════════════════════════════════


class TestTransientStrictlyBelowSteadyState:
    """The transient solution must stay strictly below the steady-state
    for any finite time.  Both use matched 3-D Green's functions
    (erfc/r for transient, 1/r for steady-state)."""

    CONFIGS = [
        (0.9, 1.2, 60),
        (0.5, 1.5, 90),
        (0.7, 0.9, 30),
        (1.0, 1.1, 45),
        (0.6, 2.0, 75),
    ]

    @pytest.mark.parametrize("d_up, d_low, angle", CONFIGS)
    def test_never_exceeds_steady_state(self, mv_cable, lv_cable, d_up, d_low, angle):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        cx = CableCrossing(
            cable_upper=mv_cable, cable_lower=lv_cable,
            depth_upper=d_up, depth_lower=d_low,
            crossing_angle_deg=angle, soil=SOIL_STANDARD,
        )
        dT_ss = cx.temperature_rise_at_upper(W)
        assert dT_ss > 0

        times = [60, 600, 3600, 86400, 7 * 86400, 30 * 86400, 365 * 86400]
        for t in times:
            dT_t = cx.transient_temperature_rise_at_upper(W, t)
            assert dT_t <= dT_ss, (
                f"Transient ΔT ({dT_t:.6f}) exceeds steady-state ({dT_ss:.6f}) "
                f"at t={t}s for depths=({d_up},{d_low}), angle={angle}°"
            )

    @pytest.mark.parametrize("d_up, d_low, angle", CONFIGS)
    def test_strictly_less_at_short_time(self, mv_cable, lv_cable, d_up, d_low, angle):
        """At 1 hour the transient must be well below steady-state."""
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        cx = CableCrossing(
            cable_upper=mv_cable, cable_lower=lv_cable,
            depth_upper=d_up, depth_lower=d_low,
            crossing_angle_deg=angle, soil=SOIL_STANDARD,
        )
        dT_ss = cx.temperature_rise_at_upper(W)
        dT_1h = cx.transient_temperature_rise_at_upper(W, 3600.0)
        assert dT_1h < dT_ss, (
            f"1-hour transient ({dT_1h:.6f}) should be strictly < "
            f"steady-state ({dT_ss:.6f})"
        )


# ═══════════════════════════════════════════════════════════════════════
#  12. ΔT at t = 0 is exactly zero
# ═══════════════════════════════════════════════════════════════════════


class TestTransientAtTimeZero:
    def test_upper_at_t0(self, standard_crossing, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        assert standard_crossing.transient_temperature_rise_at_upper(W, 0.0) == 0.0

    def test_lower_at_t0(self, standard_crossing, mv_cable):
        W = mv_cable.total_heat_per_length(I_MV, 70.0)
        assert standard_crossing.transient_temperature_rise_at_lower(W, 0.0) == 0.0

    def test_negative_time_returns_zero(self, standard_crossing, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        assert standard_crossing.transient_temperature_rise_at_upper(W, -1.0) == 0.0


# ═══════════════════════════════════════════════════════════════════════
#  13. Zero heat rate → zero ΔT
# ═══════════════════════════════════════════════════════════════════════


class TestZeroHeatRate:
    def test_steady_state(self, standard_crossing):
        assert standard_crossing.temperature_rise_at_upper(0.0) == 0.0
        assert standard_crossing.temperature_rise_at_lower(0.0) == 0.0

    def test_transient(self, standard_crossing):
        assert standard_crossing.transient_temperature_rise_at_upper(0.0, 86400) == 0.0
        assert standard_crossing.transient_temperature_rise_at_lower(0.0, 86400) == 0.0


# ═══════════════════════════════════════════════════════════════════════
#  14. Steady-state ΔT scales exactly as 1/λ
# ═══════════════════════════════════════════════════════════════════════


class TestInverseConductivityScaling:
    """The steady-state formula is ΔT = W/(4πλ) · ∫(…)ds, so ΔT ∝ 1/λ
    when all other parameters (including the integrand) are held constant."""

    def test_doubling_lambda_halves_dT(self, mv_cable, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)

        soil_1 = ThermalMaterial.from_conductivity("λ=1", 1.0, 1.8e6)
        soil_2 = ThermalMaterial.from_conductivity("λ=2", 2.0, 1.8e6)

        cx1 = CableCrossing(
            cable_upper=mv_cable, cable_lower=lv_cable,
            depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
            crossing_angle_deg=60, soil=soil_1,
        )
        cx2 = CableCrossing(
            cable_upper=mv_cable, cable_lower=lv_cable,
            depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
            crossing_angle_deg=60, soil=soil_2,
        )

        dT_1 = cx1.temperature_rise_at_upper(W)
        dT_2 = cx2.temperature_rise_at_upper(W)

        assert dT_2 == pytest.approx(dT_1 / 2.0, rel=1e-10)


# ═══════════════════════════════════════════════════════════════════════
#  15. Deeper burial (constant separation) increases crossing ΔT
# ═══════════════════════════════════════════════════════════════════════


class TestDepthEffect:
    """With constant vertical separation, pushing both cables deeper
    increases h_sum while keeping dh the same.  The image correction
    (ground surface cooling) weakens, so the net ΔT rises."""

    DEPTH_PAIRS = [
        (0.5, 0.8),   # shallow, h_sum = 1.3 — strong image correction
        (0.9, 1.2),   # standard, h_sum = 2.1
        (1.5, 1.8),   # deep, h_sum = 3.3
        (2.5, 2.8),   # very deep, h_sum = 5.3 — weak image correction
    ]

    def test_deeper_burial_higher_dT(self, mv_cable, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dTs = []
        for d_up, d_low in self.DEPTH_PAIRS:
            cx = CableCrossing(
                cable_upper=mv_cable, cable_lower=lv_cable,
                depth_upper=d_up, depth_lower=d_low,
                crossing_angle_deg=60, soil=SOIL_STANDARD,
            )
            dTs.append(cx.temperature_rise_at_upper(W))

        for i in range(1, len(dTs)):
            assert dTs[i] > dTs[i - 1], (
                f"ΔT at depths {self.DEPTH_PAIRS[i]} ({dTs[i]:.4f}) should be "
                f"> ΔT at depths {self.DEPTH_PAIRS[i-1]} ({dTs[i-1]:.4f})"
            )


# ═══════════════════════════════════════════════════════════════════════
#  16. Crossing angle ≈ 0° (parallel) must be rejected
# ═══════════════════════════════════════════════════════════════════════


class TestInvalidAngleRejected:
    def test_zero_angle_raises(self, mv_cable, lv_cable):
        with pytest.raises(ValueError, match="parallel"):
            CableCrossing(
                cable_upper=mv_cable, cable_lower=lv_cable,
                depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
                crossing_angle_deg=0.0, soil=SOIL_STANDARD,
            )

    def test_near_zero_angle_raises(self, mv_cable, lv_cable):
        with pytest.raises(ValueError, match="parallel"):
            CableCrossing(
                cable_upper=mv_cable, cable_lower=lv_cable,
                depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
                crossing_angle_deg=1e-5, soil=SOIL_STANDARD,
            )


# ═══════════════════════════════════════════════════════════════════════
#  17. Integration domain convergence — result is stable as L grows
# ═══════════════════════════════════════════════════════════════════════


class TestIntegrationConvergence:
    def test_steady_state_converges_with_half_length(self, standard_crossing, lv_cable):
        """Increasing integration_half_length beyond 50 m should change
        the result by less than 0.1 %."""
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dT_50 = standard_crossing.temperature_rise_at_upper(W, integration_half_length=50.0)
        dT_200 = standard_crossing.temperature_rise_at_upper(W, integration_half_length=200.0)

        assert dT_200 == pytest.approx(dT_50, rel=1e-3)

    def test_short_domain_underestimates(self, standard_crossing, lv_cable):
        """A very short integration domain should give a lower ΔT."""
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dT_5 = standard_crossing.temperature_rise_at_upper(W, integration_half_length=5.0)
        dT_50 = standard_crossing.temperature_rise_at_upper(W, integration_half_length=50.0)

        assert dT_5 < dT_50


# ═══════════════════════════════════════════════════════════════════════
#  18. Derating factor consistent with crossing ΔT
# ═══════════════════════════════════════════════════════════════════════


class TestDeratingConsistency:
    """The derating factor for a real crossing ΔT should always be < 1,
    and a shallower crossing angle (more heating) should give a lower
    derating factor."""

    def test_real_crossing_derates(self, mv_cable, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        cx = CableCrossing(
            cable_upper=mv_cable, cable_lower=lv_cable,
            depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
            crossing_angle_deg=60, soil=SOIL_STANDARD,
        )
        dT = cx.temperature_rise_at_upper(W)

        df = crossing_derating_factor(mv_cable, DEPTH_UPPER, dT, SOIL_STANDARD)
        assert 0.0 < df < 1.0

    def test_larger_crossing_dT_gives_lower_derating(self, mv_cable, lv_cable):
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dfs = []
        for angle in [30, 60, 90]:
            cx = CableCrossing(
                cable_upper=mv_cable, cable_lower=lv_cable,
                depth_upper=DEPTH_UPPER, depth_lower=DEPTH_LOWER,
                crossing_angle_deg=angle, soil=SOIL_STANDARD,
            )
            dT = cx.temperature_rise_at_upper(W)
            dfs.append(crossing_derating_factor(mv_cable, DEPTH_UPPER, dT, SOIL_STANDARD))

        for i in range(1, len(dfs)):
            assert dfs[i] > dfs[i - 1], (
                "Derating factor should increase with angle (less heating)"
            )


# ═══════════════════════════════════════════════════════════════════════
#  19. Transient ΔT monotonically increasing for the lower cable too
# ═══════════════════════════════════════════════════════════════════════


class TestTransientMonotonicBothCables:
    """Test 3 only checked the upper cable.  Verify the lower cable
    independently to catch asymmetric bugs."""

    TIMES_S = [60, 600, 3600, 6 * 3600, 86400, 7 * 86400, 30 * 86400]

    def test_lower_cable_monotonic(self, standard_crossing, mv_cable):
        W = mv_cable.total_heat_per_length(I_MV, 70.0)
        dTs = [
            standard_crossing.transient_temperature_rise_at_lower(W, t)
            for t in self.TIMES_S
        ]
        for i in range(1, len(dTs)):
            assert dTs[i] >= dTs[i - 1] - 1e-6

    def test_lower_cable_bounded(self, standard_crossing, mv_cable):
        W = mv_cable.total_heat_per_length(I_MV, 70.0)
        dT_ss = standard_crossing.temperature_rise_at_lower(W)
        for t in self.TIMES_S:
            dT_t = standard_crossing.transient_temperature_rise_at_lower(W, t)
            assert dT_t <= dT_ss, (
                f"Transient ΔT at lower ({dT_t:.6f}) exceeds "
                f"steady-state ({dT_ss:.6f}) at t={t}s"
            )


# ═══════════════════════════════════════════════════════════════════════
#  20. Order-of-magnitude sanity check
# ═══════════════════════════════════════════════════════════════════════


class TestOrderOfMagnitude:
    """The crossing ΔT for typical configurations should be in the
    single-digit °C range — not millikelvins and not hundreds of °C."""

    def test_typical_steady_state_range(self, standard_crossing, lv_cable, mv_cable):
        W_lv = lv_cable.total_heat_per_length(I_LV, 55.0)
        W_mv = mv_cable.total_heat_per_length(I_MV, 70.0)

        dT_upper = standard_crossing.temperature_rise_at_upper(W_lv)
        dT_lower = standard_crossing.temperature_rise_at_lower(W_mv)

        assert 0.1 < dT_upper < 50.0, f"ΔT at upper = {dT_upper:.2f} °C"
        assert 0.1 < dT_lower < 50.0, f"ΔT at lower = {dT_lower:.2f} °C"

    def test_typical_derating_is_moderate(self, standard_crossing, mv_cable, lv_cable):
        """For a 60° crossing the derating should be noticeable but not
        catastrophic — roughly in the 0.80–0.99 range."""
        W = lv_cable.total_heat_per_length(I_LV, 55.0)
        dT = standard_crossing.temperature_rise_at_upper(W)
        df = crossing_derating_factor(mv_cable, DEPTH_UPPER, dT, SOIL_STANDARD)
        assert 0.80 < df < 1.0, f"Derating factor = {df:.3f}"
