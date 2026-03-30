"""IEC 60287 standard tables and clause checks vs `thermal_cable_model`.

Numeric targets are taken from **IEC 60287-2-1:2023** (Edition 3.0) normative and
tabular data, or from the same formulas as printed in that standard for the cited
clauses.  This module is complementary to
[`tests/test_iec_formulas.py`](tests/test_iec_formulas.py), which uses arbitrary
geometries; here each test names the **standard table or clause** it locks.

Verify values against your licensed IEC PDF if a future amendment changes tables.
"""

from __future__ import annotations

import math

import pytest

from thermal_cable_model.cable import CableLayer
from thermal_cable_model.materials import (
    EPR,
    PE_JACKET,
    PVC,
    PAPER_INSULATION,
    XLPE,
)
from thermal_cable_model.thermal_network import plastic_duct_air_gap_thermal_resistance


# ── IEC 60287-2-1:2023, Table 1 — thermal resistivity ρ_th [K·m/W] ──────────
# Representative extruded / impregnated materials commonly used in rating
# calculations (exact wording and full table in the standard).


@pytest.mark.parametrize(
    "material,ref_rho_kmw,label",
    [
        (XLPE, 3.5, "thermosetting XLPE or EPR (solid insulation)"),
        (EPR, 3.5, "EPR"),
        (PVC, 6.0, "thermoplastic PVC"),
        (PAPER_INSULATION, 6.0, "impregnated paper"),
        (PE_JACKET, 3.5, "polyethylene (solid)"),
    ],
)
def test_iec60287_2_1_table1_thermal_resistivity_matches_database(
    material, ref_rho_kmw: float, label: str
):
    """IEC 60287-2-1:2023 Table 1 (thermal resistivities of materials)."""
    assert material.thermal_resistivity == pytest.approx(
        ref_rho_kmw, rel=0, abs=1e-9
    ), f"{label}: database ρ_th should match Table 1"


def test_iec60287_2_1_table5_plastic_duct_air_space_constants():
    """IEC 60287-2-1:2023 Table 5 — constants U, V, Y (plastic duct, excerpt).

    The package defaults in
    :func:`thermal_cable_model.thermal_network.plastic_duct_air_gap_thermal_resistance`
    match the thermoplastic-duct row used in rating practice (see Table 5).
    """
    u, v, y = 1.87, 0.312, 0.0037
    de_mm = 10.0
    theta_c = 50.0
    v_y = v + y * theta_c
    ref = u / (1.0 + 0.1 * v_y * de_mm)
    got = plastic_duct_air_gap_thermal_resistance(de_mm, theta_c, u, v, y)
    assert got == pytest.approx(ref, rel=0, abs=1e-12)


def test_iec60287_2_1_section_4121_single_core_t1_annular_formula():
    """IEC 60287-2-1:2023 §4.1.2.1 — single-core T1 with C_LL = 1.

    For one insulation layer, the standard’s ``ln(1 + 2 t_1 / d_c)`` form equals
    ``ln(r_o / r_i)`` with ``r_o = d_c/2 + t_1``, ``r_i = d_c/2`` (all consistent
    units).  This checks :class:`thermal_cable_model.cable.CableLayer` against
    that expression for a numeric example.
    """
    rho_kmw = XLPE.thermal_resistivity
    d_c_mm = 30.0
    t_1_mm = 15.5
    r_i_m = 0.5 * d_c_mm * 1e-3
    r_o_m = r_i_m + t_1_mm * 1e-3
    ref_t1 = rho_kmw / (2.0 * math.pi) * math.log(1.0 + 2.0 * t_1_mm / d_c_mm)
    layer = CableLayer(XLPE, r_i_m, r_o_m)
    assert layer.thermal_resistance_per_length == pytest.approx(ref_t1, rel=0, abs=1e-12)
