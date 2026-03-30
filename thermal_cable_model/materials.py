"""Thermal material property database.

Properties sourced from IEC 60287-2-1 and standard reference tables.
All quantities are in SI units.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ThermalMaterial:
    """Isotropic thermal material.

    Parameters
    ----------
    name : str
        Human-readable identifier.
    thermal_conductivity : float
        λ in W/(m·K).
    thermal_resistivity : float
        ρ = 1/λ in (K·m)/W.
    volumetric_heat_capacity : float
        ρ_m · c_p in J/(m³·K).
    """

    name: str
    thermal_conductivity: float
    thermal_resistivity: float
    volumetric_heat_capacity: float

    @classmethod
    def from_conductivity(
        cls,
        name: str,
        conductivity: float,
        volumetric_heat_capacity: float,
    ) -> ThermalMaterial:
        return cls(
            name=name,
            thermal_conductivity=conductivity,
            thermal_resistivity=1.0 / conductivity,
            volumetric_heat_capacity=volumetric_heat_capacity,
        )

    @property
    def thermal_diffusivity(self) -> float:
        """α = λ / (ρ·c_p) in m²/s."""
        return self.thermal_conductivity / self.volumetric_heat_capacity


# ── Conductor materials ──────────────────────────────────────────────
COPPER = ThermalMaterial.from_conductivity(
    "Copper", conductivity=385.0, volumetric_heat_capacity=3.45e6
)
ALUMINUM = ThermalMaterial.from_conductivity(
    "Aluminium", conductivity=230.0, volumetric_heat_capacity=2.42e6
)

# ── Insulation materials (IEC 60287-2-1 Table 1) ────────────────────
XLPE = ThermalMaterial.from_conductivity(
    "XLPE", conductivity=1.0 / 3.5, volumetric_heat_capacity=2.4e6
)
# Semi-conducting screens (IEC 60949 / CIGRE guidance; typical ρ_th ≈ 2.5 K·m/W)
SEMI_CONDUCTING_SCREEN = ThermalMaterial.from_conductivity(
    "Semi-conducting screen", conductivity=1.0 / 2.5, volumetric_heat_capacity=2.0e6
)
PVC = ThermalMaterial.from_conductivity(
    "PVC", conductivity=1.0 / 6.0, volumetric_heat_capacity=1.7e6
)
EPR = ThermalMaterial.from_conductivity(
    "EPR", conductivity=1.0 / 3.5, volumetric_heat_capacity=2.4e6
)
PAPER_INSULATION = ThermalMaterial.from_conductivity(
    "Paper insulation", conductivity=1.0 / 6.0, volumetric_heat_capacity=2.0e6
)

# ── Jacket / bedding / filler ───────────────────────────────────────
PVC_JACKET = ThermalMaterial.from_conductivity(
    "PVC jacket", conductivity=1.0 / 6.0, volumetric_heat_capacity=1.7e6
)
PE_JACKET = ThermalMaterial.from_conductivity(
    "PE jacket", conductivity=1.0 / 3.5, volumetric_heat_capacity=2.3e6
)
POLYPROPYLENE_FILLER = ThermalMaterial.from_conductivity(
    "PP filler", conductivity=1.0 / 6.0, volumetric_heat_capacity=1.9e6
)

# ── Metallic layers ─────────────────────────────────────────────────
LEAD_SHEATH = ThermalMaterial.from_conductivity(
    "Lead sheath", conductivity=35.0, volumetric_heat_capacity=1.45e6
)
STEEL_ARMOUR = ThermalMaterial.from_conductivity(
    "Steel armour", conductivity=50.0, volumetric_heat_capacity=3.8e6
)

# ── Soil types ───────────────────────────────────────────────────────
SOIL_DRY = ThermalMaterial.from_conductivity(
    "Dry soil", conductivity=1.0 / 1.2, volumetric_heat_capacity=1.3e6
)
SOIL_WET = ThermalMaterial.from_conductivity(
    "Wet soil", conductivity=1.0 / 0.7, volumetric_heat_capacity=2.5e6
)
SOIL_STANDARD = ThermalMaterial.from_conductivity(
    "Standard soil (ρ=1.0)", conductivity=1.0, volumetric_heat_capacity=1.8e6
)
THERMAL_BACKFILL = ThermalMaterial.from_conductivity(
    "Thermal backfill", conductivity=1.0 / 0.5, volumetric_heat_capacity=2.2e6
)
CONCRETE_DUCT = ThermalMaterial.from_conductivity(
    "Concrete duct", conductivity=1.2, volumetric_heat_capacity=2.0e6
)

# ── Air (still, for duct calculations) ──────────────────────────────
AIR_STILL = ThermalMaterial.from_conductivity(
    "Still air", conductivity=0.025, volumetric_heat_capacity=1200.0
)

MATERIALS: dict[str, ThermalMaterial] = {
    m.name: m
    for m in [
        COPPER, ALUMINUM,
        XLPE, SEMI_CONDUCTING_SCREEN, PVC, EPR, PAPER_INSULATION,
        PVC_JACKET, PE_JACKET, POLYPROPYLENE_FILLER,
        LEAD_SHEATH, STEEL_ARMOUR,
        SOIL_DRY, SOIL_WET, SOIL_STANDARD, THERMAL_BACKFILL, CONCRETE_DUCT,
        AIR_STILL,
    ]
}
