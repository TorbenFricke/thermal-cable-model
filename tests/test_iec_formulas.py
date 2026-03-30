"""Tier-1 validation: IEC 60287 building blocks vs closed-form reference.

Each test compares library output to hand-computed values from the same
equations documented in IEC 60287-1-1 / IEC 60287-2-1.  No proprietary
benchmark data is required.
"""

from __future__ import annotations

import math

import pytest

from thermal_cable_model.cable import Cable, CableLayer
from thermal_cable_model.cable import _dielectric_loss
from thermal_cable_model.ground import ConstantGroundTemperature
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import COPPER, SOIL_STANDARD, XLPE
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation
from thermal_cable_model.thermal_network import (
    InternalThermalResistances,
    _van_wormer,
    external_thermal_resistance,
    mutual_heating_resistance,
)


def _iec60287_1_1_conductor_rise(
    Wc: float,
    Wd: float,
    T1: float,
    T2: float,
    T3: float,
    T4: float,
    n: int,
    lambda1: float,
    lambda2: float,
) -> float:
    """Temperature rise of the conductor above ambient [K], IEC 60287-1-1."""
    return (
        (Wc + 0.5 * Wd) * T1
        + (Wc * (1.0 + lambda1) + Wd) * T2
        + (Wc * (1.0 + lambda1 + lambda2) + Wd) * T3
        + n * (Wc * (1.0 + lambda1 + lambda2) + Wd) * T4
    )


class TestThermalResistanceT1T2T3:
    """Cylindrical layer: T = rho/(2 pi) ln(r2/r1)."""

    def test_single_layer_matches_formula(self):
        r1, r2 = 0.010, 0.020
        layer = CableLayer(XLPE, r1, r2)
        expected = XLPE.thermal_resistivity / (2.0 * math.pi) * math.log(r2 / r1)
        assert layer.thermal_resistance_per_length == pytest.approx(expected, rel=0, abs=1e-12)

    def test_internal_split_two_layers(self):
        """from_cable: two layers → T1 first, T3 last, T2 = 0."""
        A = 50e-6
        r_c = math.sqrt(A / math.pi)
        layers = [
            CableLayer(XLPE, r_c, r_c + 2e-3),
            CableLayer(XLPE, r_c + 2e-3, r_c + 3e-3),
        ]
        cable = Cable(
            name="test",
            voltage_class="LV",
            n_conductors=1,
            conductor_area=A,
            conductor_material=COPPER,
            layers=layers,
            ac_resistance_20c=1e-3,
        )
        tr = InternalThermalResistances.from_cable(cable)
        assert tr.T2 == 0.0
        assert tr.T1 == pytest.approx(layers[0].thermal_resistance_per_length, rel=0, abs=1e-12)
        assert tr.T3 == pytest.approx(layers[1].thermal_resistance_per_length, rel=0, abs=1e-12)

    def test_internal_split_three_plus_layers(self):
        """Three layers: T1 first, T2 sum of middle, T3 last."""
        r0 = 0.005
        layers = [
            CableLayer(XLPE, r0, r0 + 1e-3),
            CableLayer(XLPE, r0 + 1e-3, r0 + 2e-3),
            CableLayer(XLPE, r0 + 2e-3, r0 + 3e-3),
        ]
        cable = Cable(
            name="test",
            voltage_class="LV",
            n_conductors=1,
            conductor_area=math.pi * r0**2,
            conductor_material=COPPER,
            layers=layers,
            ac_resistance_20c=1e-3,
        )
        tr = InternalThermalResistances.from_cable(cable)
        assert tr.T1 == pytest.approx(layers[0].thermal_resistance_per_length, rel=0, abs=1e-12)
        assert tr.T2 == pytest.approx(
            layers[1].thermal_resistance_per_length, rel=0, abs=1e-12
        )
        assert tr.T3 == pytest.approx(layers[2].thermal_resistance_per_length, rel=0, abs=1e-12)


class TestExternalResistanceT4:
    """IEC 60287-2-1 §2.2.7: T4 = rho/(2 pi) acosh(2L/De)."""

    def test_known_depth_and_diameter(self):
        depth = 1.0
        outer_r = 0.025
        rho = 1.0
        from thermal_cable_model.materials import ThermalMaterial

        soil = ThermalMaterial("ref", 1.0 / rho, rho, 1.8e6)
        u = 2.0 * depth / (2.0 * outer_r)
        expected = rho / (2.0 * math.pi) * math.log(u + math.sqrt(u**2 - 1.0))
        t4 = external_thermal_resistance(depth, outer_r, soil)
        assert t4 == pytest.approx(expected, rel=0, abs=1e-12)

    def test_invalid_geometry_raises(self):
        from thermal_cable_model.materials import ThermalMaterial

        soil = ThermalMaterial("ref", 1.0, 1.0, 1.8e6)
        with pytest.raises(ValueError):
            external_thermal_resistance(0.02, 0.02, soil)


class TestVanWormerCoefficient:
    """p = 1/(2 ln(r2/r1)) - 1/((r2/r1)^2 - 1)."""

    def test_radius_ratio_two(self):
        r1, r2 = 0.01, 0.02
        layer = CableLayer(XLPE, r1, r2)
        ratio = r2 / r1
        expected = 1.0 / (2.0 * math.log(ratio)) - 1.0 / (ratio**2 - 1.0)
        assert _van_wormer(layer) == pytest.approx(expected, rel=0, abs=1e-12)

    def test_thin_layer_defaults_to_half(self):
        layer = CableLayer(XLPE, 0.01, 0.0100005)
        assert _van_wormer(layer) == pytest.approx(0.5, abs=0.01)


class TestACResistance:
    def test_linear_temperature_correction(self):
        cable = Cable.single_core_xlpe_cu(240, voltage_kv=20.0)
        r20 = cable.ac_resistance_20c
        alpha = cable.temp_coeff_resistance
        for t in (20.0, 50.0, 90.0):
            expected = r20 * (1.0 + alpha * (t - 20.0))
            assert cable.ac_resistance(t) == pytest.approx(expected, rel=0, abs=1e-15)


class TestMutualHeatingResistance:
    """Rm = rho/(2 pi) ln(d'/d)."""

    def test_horizontal_spacing(self):
        rho = 1.0
        from thermal_cable_model.materials import ThermalMaterial

        soil = ThermalMaterial("ref", 1.0 / rho, rho, 1.8e6)
        x1, y1 = 0.0, 1.0
        x2, y2 = 0.3, 1.0
        d = 0.3
        d_image = math.sqrt(0.3**2 + (y1 + y2) ** 2)
        expected = rho / (2.0 * math.pi) * math.log(d_image / d)
        rm = mutual_heating_resistance(x1, y1, x2, y2, soil)
        assert rm == pytest.approx(expected, rel=0, abs=1e-12)

    def test_same_position_zero(self):
        from thermal_cable_model.materials import ThermalMaterial

        soil = ThermalMaterial("ref", 1.0, 1.0, 1.8e6)
        assert mutual_heating_resistance(0.0, 1.0, 0.0, 1.0, soil) == 0.0


class TestDielectricLoss:
    def test_matches_manual_xlpe_formula(self):
        voltage_kv = 20.0
        t_ins = 5.5e-3
        r_c = 0.01
        wd = _dielectric_loss(voltage_kv, t_ins, r_c)
        eps_r = 2.5
        tan_delta = 4e-4
        omega = 2.0 * math.pi * 50.0
        u0 = voltage_kv * 1e3 / math.sqrt(3.0)
        r_i = r_c + t_ins
        c_lin = 2.0 * math.pi * eps_r * 8.854e-12 / math.log(r_i / r_c)
        expected = omega * c_lin * u0**2 * tan_delta
        assert wd == pytest.approx(expected, rel=0, abs=1e-18)

    def test_below_3kv_zero(self):
        assert _dielectric_loss(1.0, 2e-3, 0.005) == 0.0


class TestSteadyStateRatingEquation:
    """Network steady-state conductor temperature vs IEC 60287-1-1 rise equation."""

    @pytest.mark.parametrize("current_a", [100.0, 400.0, 600.0])
    def test_single_core_xlpe_matches_closed_form(self, current_a):
        cable = Cable.single_core_xlpe_cu(240, voltage_kv=20.0)
        ground = ConstantGroundTemperature(15.0)
        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(cable, 0.0, 1.2, LoadProfile.constant(current_a, 3600.0))
        sim = ThermalSimulation(inst)
        ss = sim.run_steady_state()
        tc = float(ss.conductor_temps[0, 0])
        tamb = float(ss.ambient_temps[0, 0])
        tr = InternalThermalResistances.from_cable(cable)
        t4 = external_thermal_resistance(1.2, cable.outer_radius, SOIL_STANDARD)
        wc = cable.conductor_loss(current_a, tc)
        wd = cable.dielectric_loss
        n = cable.n_conductors
        rise = _iec60287_1_1_conductor_rise(
            wc,
            wd,
            tr.T1,
            tr.T2,
            tr.T3,
            t4,
            n,
            cable.loss_factor_sheath,
            cable.loss_factor_armour,
        )
        assert tc == pytest.approx(tamb + rise, abs=0.05)

    def test_three_core_matches_closed_form(self):
        cable = Cable.three_core_xlpe_cu(150, voltage_kv=0.6)
        ground = ConstantGroundTemperature(15.0)
        I = 250.0
        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(cable, 0.0, 1.0, LoadProfile.constant(I, 3600.0))
        sim = ThermalSimulation(inst)
        ss = sim.run_steady_state()
        tc = float(ss.conductor_temps[0, 0])
        tamb = float(ss.ambient_temps[0, 0])
        tr = InternalThermalResistances.from_cable(cable)
        t4 = external_thermal_resistance(1.0, cable.outer_radius, SOIL_STANDARD)
        wc = cable.conductor_loss(I, tc)
        wd = cable.dielectric_loss
        n = cable.n_conductors
        rise = _iec60287_1_1_conductor_rise(
            wc,
            wd,
            tr.T1,
            tr.T2,
            tr.T3,
            t4,
            n,
            cable.loss_factor_sheath,
            cable.loss_factor_armour,
        )
        assert tc == pytest.approx(tamb + rise, abs=0.05)
