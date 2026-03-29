Flat vs Trefoil Formation
=========================

**Script:** ``examples/flat_vs_trefoil.py``

Compares two common installation geometries for a symmetric 3-phase cable
system: flat formation and trefoil formation.

Setup
-----

- **Cable:** 1 × 240 mm\ :sup:`2` Cu XLPE, 20 kV (three identical single-core
  cables per formation)
- **Load:** balanced 400 A per phase, constant for 48 hours
- **Soil:** standard soil
- **Nominal depth:** 1.0 m (centroid for trefoil)

Flat formation
  Three cables side-by-side, touching, at 1.0 m depth.

Trefoil formation
  Three cables in an equilateral triangle (touching), centroid at 1.0 m depth.

What it demonstrates
--------------------

1. Setting up two different geometric configurations
2. Comparing peak conductor temperatures between formations
3. Assessing thermal asymmetry between phases (the centre cable in flat
   formation is hotter than the outer cables)
4. Custom matplotlib plots for direct phase-by-phase comparison

Output files
------------

- ``flat_vs_trefoil_cross_section.png`` — side-by-side cross-section comparison
- ``flat_vs_trefoil_temperatures.png`` — conductor temperature traces for all
  six phases
