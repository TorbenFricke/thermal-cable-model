"""CIGRE TB 880 benchmarks vs `thermal_cable_model`.

Reference numbers are the printed values from the open CIGRE TB880 notebooks
(case 01 introductory, case 02 trefoil in HDPE ducts):

https://github.com/frdmendoza/cbl_CIGRE_TB880

Every asserted “computed” quantity is obtained from `thermal_cable_model`
(`CableLayer`, `iec_line_capacitance_per_length`, `dielectric_loss_per_length`,
`external_thermal_resistance`, duct helpers in `thermal_network`, etc.), not
from re-implementing the notebook’s formula bodies in this file.

The TB880 case 01 notebook used a different IEC clause for external thermal
resistance (IEC 60287-2-1 §4.2.4.3.2) than
`external_thermal_resistance` (§2.2.7 / ``ln(u + √(u²−1))``); a constant
`REF_CASE01_T4_NOTEBOOK_KMW` records the notebook value for documentation only.
"""

from __future__ import annotations

import math

import pytest

from thermal_cable_model.cable import (
    Cable,
    CableLayer,
    dielectric_loss_per_length,
    iec_line_capacitance_per_length,
)
from thermal_cable_model.ground import ConstantGroundTemperature
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import (
    ALUMINUM,
    COPPER,
    PE_JACKET,
    SEMI_CONDUCTING_SCREEN,
    SOIL_STANDARD,
    XLPE,
    ThermalMaterial,
)
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation
from thermal_cable_model.thermal_network import (
    InternalThermalResistances,
    cylindrical_shell_thermal_resistance,
    cylindrical_shell_thermal_resistance_diameters,
    external_buried_duct_thermal_resistance_tb880_form,
    external_thermal_resistance,
    plastic_duct_air_gap_thermal_resistance,
)

# ── Geometry (mm) from TB880 open notebooks ───────────────────────────────

TB880_THK_MM = [0.0, 1.5, 15.5, 1.3, 0.8, 3.5]
TB880_DIA0_MM = 30.3
TB880_L_MM = 1000.0
TB880_RHO_SOIL_KMW = 1.0

TB880_CASE02_D0_MM = 140.0
TB880_CASE02_DD_MM = 119.4
TB880_CASE02_U_AIR = 1.87
TB880_CASE02_V_AIR = 0.312
TB880_CASE02_Y_AIR = 0.0037
TB880_CASE02_THETA_OP_C = 70.0
TB880_RHO_DUCT_WALL_KMW = 3.5

# ── Oracle constants (notebook print-outs) ────────────────────────────────

REF_CASE01_T1_KMW = 0.4198714890
REF_CASE01_T3_KMW = 0.0867193748
REF_CASE01_T4_NOTEBOOK_KMW = 1.5946928925

REF_CASE02_T3_RAW_KMW = 0.0541996092
REF_CASE02_T4_PRIME_KMW = 0.3520961015
REF_CASE02_T4_DOUBLE_PRIME_KMW = 0.08866064719
REF_CASE02_T4_TRIPLE_PRIME_KMW = 1.38002093961
REF_CASE02_T4_TOTAL_KMW = 1.82077768833

REF_CASE02_CAP_F_PER_M = 2.1107662202e-10
REF_CASE02_WD_W_PER_M = 0.3851382172


def tb880_diameters_mm() -> list[float]:
    dia = [0.0] * 6
    dia[0] = TB880_DIA0_MM
    for i in range(1, 6):
        dia[i] = dia[i - 1] + 2.0 * TB880_THK_MM[i]
    return dia


def tb880_radii_m() -> list[float]:
    return [0.5 * d * 1e-3 for d in tb880_diameters_mm()]


def tb880_t1_cable_layers() -> list[CableLayer]:
    r = tb880_radii_m()
    return [
        CableLayer(SEMI_CONDUCTING_SCREEN, r[0], r[1]),
        CableLayer(XLPE, r[1], r[2]),
        CableLayer(SEMI_CONDUCTING_SCREEN, r[2], r[3]),
    ]


def tb880_oversheath_layer() -> CableLayer:
    r = tb880_radii_m()
    return CableLayer(PE_JACKET, r[4], r[5])


def tb880_soil_rho1() -> ThermalMaterial:
    return ThermalMaterial(
        name="TB880 soil ρ=1",
        thermal_conductivity=1.0 / TB880_RHO_SOIL_KMW,
        thermal_resistivity=TB880_RHO_SOIL_KMW,
        volumetric_heat_capacity=1.8e6,
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
    return (
        (Wc + 0.5 * Wd) * T1
        + (Wc * (1.0 + lambda1) + Wd) * T2
        + (Wc * (1.0 + lambda1 + lambda2) + Wd) * T3
        + n * (Wc * (1.0 + lambda1 + lambda2) + Wd) * T4
    )


def tb880_like_single_core_cable() -> Cable:
    """630 mm²–style stack: screens + XLPE + Al sheath + PE (TB880 diameters)."""
    dia = tb880_diameters_mm()
    r = tb880_radii_m()
    A = math.pi * r[0] ** 2
    layers = [
        CableLayer(SEMI_CONDUCTING_SCREEN, r[0], r[1]),
        CableLayer(XLPE, r[1], r[2]),
        CableLayer(SEMI_CONDUCTING_SCREEN, r[2], r[3]),
        CableLayer(ALUMINUM, r[3], r[4]),
        CableLayer(PE_JACKET, r[4], r[5]),
    ]
    c_lin = iec_line_capacitance_per_length(2.5, dia[1], dia[2])
    wd = dielectric_loss_per_length(c_lin, 132.0, 0.001, 50.0)
    return Cable(
        name="TB880-style 630 mm² single-core",
        voltage_class="MV",
        n_conductors=1,
        conductor_area=A,
        conductor_material=COPPER,
        layers=layers,
        ac_resistance_20c=28.3e-6,
        temp_coeff_resistance=3.93e-3,
        max_conductor_temp=90.0,
        loss_factor_sheath=0.0,
        loss_factor_armour=0.0,
        dielectric_loss=wd,
    )


class TestTb880Case01AgainstThermalCableModel:
    """Case 01 geometry: layer T1/T3 and package T4 vs TB880 reference constants."""

    def test_t1_sum_matches_reference(self):
        layers = tb880_t1_cable_layers()
        t1 = sum(L.thermal_resistance_per_length for L in layers)
        assert t1 == pytest.approx(REF_CASE01_T1_KMW, rel=0, abs=1e-6)

    def test_t1_each_layer_matches_cylindrical_shell_helper(self):
        layers = tb880_t1_cable_layers()
        for L in layers:
            expected = cylindrical_shell_thermal_resistance(
                L.material.thermal_resistivity,
                L.inner_radius,
                L.outer_radius,
            )
            assert L.thermal_resistance_per_length == pytest.approx(
                expected, rel=0, abs=1e-12
            )

    def test_t3_oversheath_raw_matches_case02_reference(self):
        L = tb880_oversheath_layer()
        assert L.thermal_resistance_per_length == pytest.approx(
            REF_CASE02_T3_RAW_KMW, rel=0, abs=1e-6
        )

    def test_t3_case01_iec_42443_factor_1_6_on_raw_only(self):
        """IEC 60287-2-1 §4.2.4.3.2: multiply jacket thermal resistance by 1.6."""
        L = tb880_oversheath_layer()
        assert 1.6 * L.thermal_resistance_per_length == pytest.approx(
            REF_CASE01_T3_KMW, rel=0, abs=1e-6
        )

    def test_external_t4_matches_documented_acosh_closed_form(self):
        dia = tb880_diameters_mm()
        depth_m = TB880_L_MM * 1e-3
        outer_r_m = 0.5 * dia[5] * 1e-3
        soil = tb880_soil_rho1()
        t4_pkg = external_thermal_resistance(depth_m, outer_r_m, soil)
        u = 2.0 * depth_m / (2.0 * outer_r_m)
        oracle = TB880_RHO_SOIL_KMW / (2.0 * math.pi) * math.log(
            u + math.sqrt(u**2 - 1.0)
        )
        assert t4_pkg == pytest.approx(oracle, rel=0, abs=1e-9)
        assert abs(t4_pkg - REF_CASE01_T4_NOTEBOOK_KMW) > 0.3


class TestTb880Case02DuctChainAgainstThermalCableModel:
    """Case 02 duct resistances and dielectric quantities via package APIs."""

    def test_t1_same_geometry_as_case01(self):
        t1 = sum(L.thermal_resistance_per_length for L in tb880_t1_cable_layers())
        assert t1 == pytest.approx(REF_CASE01_T1_KMW, rel=0, abs=1e-6)

    def test_t3_raw_matches_reference(self):
        L = tb880_oversheath_layer()
        assert L.thermal_resistance_per_length == pytest.approx(
            REF_CASE02_T3_RAW_KMW, rel=0, abs=1e-6
        )

    def test_t4_prime_plastic_duct_table4(self):
        de_mm = tb880_diameters_mm()[5]
        t4p = plastic_duct_air_gap_thermal_resistance(
            de_mm,
            TB880_CASE02_THETA_OP_C,
            TB880_CASE02_U_AIR,
            TB880_CASE02_V_AIR,
            TB880_CASE02_Y_AIR,
        )
        assert t4p == pytest.approx(REF_CASE02_T4_PRIME_KMW, rel=0, abs=1e-6)

    def test_t4_double_prime_duct_wall_ln_ratio(self):
        t4pp = cylindrical_shell_thermal_resistance_diameters(
            TB880_RHO_DUCT_WALL_KMW,
            TB880_CASE02_DD_MM,
            TB880_CASE02_D0_MM,
        )
        assert t4pp == pytest.approx(REF_CASE02_T4_DOUBLE_PRIME_KMW, rel=0, abs=1e-7)

    def test_t4_triple_prime_external_buried_duct(self):
        t4ppp = external_buried_duct_thermal_resistance_tb880_form(
            TB880_L_MM,
            TB880_CASE02_D0_MM,
            TB880_RHO_SOIL_KMW,
        )
        assert t4ppp == pytest.approx(
            REF_CASE02_T4_TRIPLE_PRIME_KMW, rel=0, abs=1e-6
        )

    def test_t4_total_sum(self):
        de_mm = tb880_diameters_mm()[5]
        t4p = plastic_duct_air_gap_thermal_resistance(
            de_mm,
            TB880_CASE02_THETA_OP_C,
            TB880_CASE02_U_AIR,
            TB880_CASE02_V_AIR,
            TB880_CASE02_Y_AIR,
        )
        t4pp = cylindrical_shell_thermal_resistance_diameters(
            TB880_RHO_DUCT_WALL_KMW,
            TB880_CASE02_DD_MM,
            TB880_CASE02_D0_MM,
        )
        t4ppp = external_buried_duct_thermal_resistance_tb880_form(
            TB880_L_MM,
            TB880_CASE02_D0_MM,
            TB880_RHO_SOIL_KMW,
        )
        assert (t4p + t4pp + t4ppp) == pytest.approx(
            REF_CASE02_T4_TOTAL_KMW, rel=0, abs=1e-6
        )

    def test_line_capacitance_per_metre(self):
        dia = tb880_diameters_mm()
        c = iec_line_capacitance_per_length(2.5, dia[1], dia[2])
        assert c == pytest.approx(REF_CASE02_CAP_F_PER_M, rel=0, abs=1e-15)

    def test_dielectric_loss_wd(self):
        dia = tb880_diameters_mm()
        c = iec_line_capacitance_per_length(2.5, dia[1], dia[2])
        wd = dielectric_loss_per_length(c, 132.0, 0.001, 50.0)
        assert wd == pytest.approx(REF_CASE02_WD_W_PER_M, rel=0, abs=1e-7)


class TestTb880SteadyStateVsIec6028711Rise:
    """Assembled network vs IEC 60287-1-1 conductor rise (same pattern as IEC tests)."""

    @pytest.mark.parametrize("current_a", [200.0, 450.0, 700.0])
    def test_tb880_like_cable_matches_closed_form(self, current_a: float):
        cable = tb880_like_single_core_cable()
        ground = ConstantGroundTemperature(20.0)
        depth_m = TB880_L_MM * 1e-3
        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(
            cable, 0.0, depth_m, LoadProfile.constant(current_a, 3600.0)
        )
        sim = ThermalSimulation(inst)
        ss = sim.run_steady_state()
        tc = float(ss.conductor_temps[0, 0])
        tamb = float(ss.ambient_temps[0, 0])
        tr = InternalThermalResistances.from_cable(cable)
        t4 = external_thermal_resistance(depth_m, cable.outer_radius, SOIL_STANDARD)
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
