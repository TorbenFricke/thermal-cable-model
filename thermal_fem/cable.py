"""Cable construction definitions for LV and MV applications.

Cables are described by concentric cylindrical layers from conductor outward.
Factory class-methods provide common cable constructions; users can also
build custom cables layer-by-layer.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from thermal_fem.materials import (
    ALUMINUM,
    COPPER,
    EPR,
    LEAD_SHEATH,
    PE_JACKET,
    PVC,
    PVC_JACKET,
    STEEL_ARMOUR,
    XLPE,
    ThermalMaterial,
)


@dataclass
class CableLayer:
    """Single concentric cylindrical layer.

    Parameters
    ----------
    material : ThermalMaterial
    inner_radius : float   [m]
    outer_radius : float   [m]
    """

    material: ThermalMaterial
    inner_radius: float
    outer_radius: float

    @property
    def thickness(self) -> float:
        return self.outer_radius - self.inner_radius

    @property
    def thermal_resistance_per_length(self) -> float:
        """Radial thermal resistance per unit length [K·m/W]."""
        if self.thickness < 1e-9:
            return 0.0
        return (
            self.material.thermal_resistivity
            / (2.0 * math.pi)
            * math.log(self.outer_radius / self.inner_radius)
        )

    @property
    def thermal_capacitance_per_length(self) -> float:
        """Thermal capacitance per unit length [J/(m·K)]."""
        area = math.pi * (self.outer_radius ** 2 - self.inner_radius ** 2)
        return self.material.volumetric_heat_capacity * area


@dataclass
class Cable:
    """Complete cable construction.

    Parameters
    ----------
    name : str
        Descriptive identifier.
    voltage_class : str
        ``"LV"`` or ``"MV"``.
    n_conductors : int
        Number of power conductors (1 for single-core, 3 for three-core).
    conductor_area : float
        Cross-sectional area of one conductor [m²].
    conductor_material : ThermalMaterial
    layers : list[CableLayer]
        Ordered list from innermost insulation outward (excluding the
        conductor itself).
    ac_resistance_20c : float
        AC resistance at 20 °C per unit length [Ω/m].
    temp_coeff_resistance : float
        Temperature coefficient of resistance [1/K] (≈ 3.93e-3 for Cu,
        4.03e-3 for Al).
    max_conductor_temp : float
        Maximum continuous conductor temperature [°C].
    loss_factor_sheath : float
        λ₁ — ratio of sheath/screen losses to conductor losses.
    loss_factor_armour : float
        λ₂ — ratio of armour losses to conductor losses.
    dielectric_loss : float
        W_d — dielectric loss per phase per unit length [W/m].
    """

    name: str
    voltage_class: str
    n_conductors: int
    conductor_area: float
    conductor_material: ThermalMaterial
    layers: list[CableLayer] = field(default_factory=list)
    ac_resistance_20c: float = 0.0
    temp_coeff_resistance: float = 3.93e-3
    max_conductor_temp: float = 90.0
    loss_factor_sheath: float = 0.0
    loss_factor_armour: float = 0.0
    dielectric_loss: float = 0.0

    # ── derived geometry ─────────────────────────────────────────────

    @property
    def conductor_radius(self) -> float:
        """Equivalent solid conductor radius [m]."""
        return math.sqrt(self.conductor_area / math.pi)

    @property
    def outer_radius(self) -> float:
        if not self.layers:
            return self.conductor_radius
        return self.layers[-1].outer_radius

    @property
    def outer_diameter(self) -> float:
        return 2.0 * self.outer_radius

    # ── electrical properties ────────────────────────────────────────

    def ac_resistance(self, temperature: float) -> float:
        """Temperature-corrected AC resistance [Ω/m]."""
        return self.ac_resistance_20c * (
            1.0 + self.temp_coeff_resistance * (temperature - 20.0)
        )

    def conductor_loss(self, current: float, temperature: float) -> float:
        """Joule loss per conductor per unit length W_c [W/m]."""
        return current ** 2 * self.ac_resistance(temperature)

    def total_heat_per_length(self, current: float, temperature: float) -> float:
        """Total heat generation per cable per unit length [W/m].

        Includes conductor, sheath, armour, and dielectric losses.
        """
        wc = self.conductor_loss(current, temperature)
        return self.n_conductors * (
            wc * (1.0 + self.loss_factor_sheath + self.loss_factor_armour)
            + self.dielectric_loss
        )

    # ── thermal properties ───────────────────────────────────────────

    @property
    def conductor_capacitance_per_length(self) -> float:
        """Thermal capacitance of the conductor [J/(m·K)]."""
        return (
            self.conductor_material.volumetric_heat_capacity
            * self.conductor_area
        )

    # ── factory methods for common constructions ─────────────────────

    @classmethod
    def single_core_xlpe_cu(
        cls,
        conductor_area_mm2: float,
        voltage_class: str = "MV",
        voltage_kv: float = 20.0,
    ) -> Cable:
        """Single-core copper XLPE cable (typical MV distribution cable).

        Approximate geometry generated from the conductor area; for exact
        modelling supply layer dimensions directly.
        """
        A = conductor_area_mm2 * 1e-6  # m²
        r_c = math.sqrt(A / math.pi)

        r_ins_inner = r_c + 0.5e-3  # semi-conducting screen
        insulation_thickness = _insulation_thickness_xlpe(voltage_kv)
        r_ins_outer = r_ins_inner + insulation_thickness

        r_screen_outer = r_ins_outer + 0.3e-3  # copper tape screen
        r_jacket_outer = r_screen_outer + 2.0e-3

        layers = [
            CableLayer(XLPE, r_c, r_ins_outer),
            CableLayer(COPPER, r_ins_outer, r_screen_outer),
            CableLayer(PE_JACKET, r_screen_outer, r_jacket_outer),
        ]

        rac20 = _rac_from_area_cu(conductor_area_mm2)

        return cls(
            name=f"1×{conductor_area_mm2:.0f} mm² Cu XLPE {voltage_kv:.0f} kV",
            voltage_class=voltage_class,
            n_conductors=1,
            conductor_area=A,
            conductor_material=COPPER,
            layers=layers,
            ac_resistance_20c=rac20,
            temp_coeff_resistance=3.93e-3,
            max_conductor_temp=90.0,
            loss_factor_sheath=0.05,
            loss_factor_armour=0.0,
            dielectric_loss=_dielectric_loss(voltage_kv, insulation_thickness, r_c),
        )

    @classmethod
    def single_core_xlpe_al(
        cls,
        conductor_area_mm2: float,
        voltage_class: str = "MV",
        voltage_kv: float = 20.0,
    ) -> Cable:
        """Single-core aluminium XLPE cable."""
        A = conductor_area_mm2 * 1e-6
        r_c = math.sqrt(A / math.pi)

        r_ins_inner = r_c + 0.5e-3
        insulation_thickness = _insulation_thickness_xlpe(voltage_kv)
        r_ins_outer = r_ins_inner + insulation_thickness

        r_screen_outer = r_ins_outer + 0.3e-3
        r_jacket_outer = r_screen_outer + 2.0e-3

        layers = [
            CableLayer(XLPE, r_c, r_ins_outer),
            CableLayer(ALUMINUM, r_ins_outer, r_screen_outer),
            CableLayer(PE_JACKET, r_screen_outer, r_jacket_outer),
        ]

        rac20 = _rac_from_area_al(conductor_area_mm2)

        return cls(
            name=f"1×{conductor_area_mm2:.0f} mm² Al XLPE {voltage_kv:.0f} kV",
            voltage_class=voltage_class,
            n_conductors=1,
            conductor_area=A,
            conductor_material=ALUMINUM,
            layers=layers,
            ac_resistance_20c=rac20,
            temp_coeff_resistance=4.03e-3,
            max_conductor_temp=90.0,
            loss_factor_sheath=0.05,
            loss_factor_armour=0.0,
            dielectric_loss=_dielectric_loss(voltage_kv, insulation_thickness, r_c),
        )

    @classmethod
    def three_core_xlpe_cu(
        cls,
        conductor_area_mm2: float,
        voltage_class: str = "LV",
        voltage_kv: float = 0.6,
    ) -> Cable:
        """Three-core copper XLPE cable (typical LV distribution cable)."""
        A = conductor_area_mm2 * 1e-6
        r_c = math.sqrt(A / math.pi)

        insulation_thickness = _insulation_thickness_xlpe(voltage_kv)
        r_ins_outer = r_c + insulation_thickness

        # Three-core equivalent: assume trefoil touching -> effective outer radius
        # enclosing three insulated cores in a round cross section
        r_trefoil = r_ins_outer * (1.0 + 2.0 / math.sqrt(3.0)) / 2.0 + 0.5e-3
        r_bedding_outer = r_trefoil + 1.5e-3
        r_armour_outer = r_bedding_outer + 2.0e-3
        r_jacket_outer = r_armour_outer + 2.0e-3

        layers = [
            CableLayer(XLPE, r_c, r_ins_outer),
            CableLayer(PVC, r_ins_outer, r_bedding_outer),
            CableLayer(STEEL_ARMOUR, r_bedding_outer, r_armour_outer),
            CableLayer(PVC_JACKET, r_armour_outer, r_jacket_outer),
        ]

        rac20 = _rac_from_area_cu(conductor_area_mm2)

        return cls(
            name=f"3×{conductor_area_mm2:.0f} mm² Cu XLPE {voltage_kv:.1f} kV",
            voltage_class=voltage_class,
            n_conductors=3,
            conductor_area=A,
            conductor_material=COPPER,
            layers=layers,
            ac_resistance_20c=rac20,
            temp_coeff_resistance=3.93e-3,
            max_conductor_temp=70.0,
            loss_factor_sheath=0.0,
            loss_factor_armour=0.03,
            dielectric_loss=0.0,
        )

    @classmethod
    def three_core_pvc_cu(
        cls,
        conductor_area_mm2: float,
        voltage_class: str = "LV",
        voltage_kv: float = 0.6,
    ) -> Cable:
        """Three-core copper PVC cable (LV)."""
        A = conductor_area_mm2 * 1e-6
        r_c = math.sqrt(A / math.pi)

        insulation_thickness = _insulation_thickness_pvc(voltage_kv)
        r_ins_outer = r_c + insulation_thickness

        r_trefoil = r_ins_outer * (1.0 + 2.0 / math.sqrt(3.0)) / 2.0 + 0.5e-3
        r_bedding_outer = r_trefoil + 1.5e-3
        r_armour_outer = r_bedding_outer + 2.0e-3
        r_jacket_outer = r_armour_outer + 2.0e-3

        layers = [
            CableLayer(PVC, r_c, r_ins_outer),
            CableLayer(PVC, r_ins_outer, r_bedding_outer),
            CableLayer(STEEL_ARMOUR, r_bedding_outer, r_armour_outer),
            CableLayer(PVC_JACKET, r_armour_outer, r_jacket_outer),
        ]

        rac20 = _rac_from_area_cu(conductor_area_mm2)

        return cls(
            name=f"3×{conductor_area_mm2:.0f} mm² Cu PVC {voltage_kv:.1f} kV",
            voltage_class=voltage_class,
            n_conductors=3,
            conductor_area=A,
            conductor_material=COPPER,
            layers=layers,
            ac_resistance_20c=rac20,
            temp_coeff_resistance=3.93e-3,
            max_conductor_temp=70.0,
            loss_factor_sheath=0.0,
            loss_factor_armour=0.03,
            dielectric_loss=0.0,
        )


# ── private helpers ──────────────────────────────────────────────────


def _insulation_thickness_xlpe(voltage_kv: float) -> float:
    """Approximate XLPE insulation thickness [m] from rated voltage."""
    if voltage_kv <= 1.0:
        return 1.0e-3
    if voltage_kv <= 6.0:
        return 3.4e-3
    if voltage_kv <= 10.0:
        return 3.4e-3
    if voltage_kv <= 20.0:
        return 5.5e-3
    if voltage_kv <= 30.0:
        return 8.0e-3
    return 10.0e-3


def _insulation_thickness_pvc(voltage_kv: float) -> float:
    if voltage_kv <= 1.0:
        return 1.5e-3
    return 2.5e-3


def _rac_from_area_cu(area_mm2: float) -> float:
    """Approximate 20 °C AC resistance for stranded copper conductor [Ω/m]."""
    rho_dc = 1.7241e-8  # Ω·m at 20 °C
    # AC effects: skin + proximity factor ≈ 1.02–1.10 for typical sizes
    ks = 1.0 + 0.02 * min(area_mm2 / 500.0, 1.0)
    return rho_dc / (area_mm2 * 1e-6) * ks


def _rac_from_area_al(area_mm2: float) -> float:
    """Approximate 20 °C AC resistance for aluminium conductor [Ω/m]."""
    rho_dc = 2.8264e-8
    ks = 1.0 + 0.02 * min(area_mm2 / 500.0, 1.0)
    return rho_dc / (area_mm2 * 1e-6) * ks


def _dielectric_loss(voltage_kv: float, insulation_thickness: float,
                     conductor_radius: float) -> float:
    """Approximate dielectric loss per phase [W/m] for XLPE."""
    if voltage_kv < 3.0:
        return 0.0
    eps_r = 2.5
    tan_delta = 4e-4  # XLPE typical
    omega = 2.0 * math.pi * 50.0
    U0 = voltage_kv * 1e3 / math.sqrt(3.0)
    r_c = conductor_radius
    r_i = r_c + insulation_thickness
    C = 2.0 * math.pi * eps_r * 8.854e-12 / math.log(r_i / r_c)
    return omega * C * U0 ** 2 * tan_delta
