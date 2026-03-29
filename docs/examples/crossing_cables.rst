Cable Crossing
==============

**Script:** ``examples/crossing_cables.py``

Analyses the mutual thermal effect when an MV cable crosses over an LV cable
at 60°.  Includes derating calculations and a comparison between isolated
and crossing-affected temperatures.

Setup
-----

- **Upper cable:** 1 × 240 mm\ :sup:`2` Cu XLPE 20 kV at 0.9 m depth, 350 A
- **Lower cable:** 3 × 150 mm\ :sup:`2` Cu XLPE 0.6 kV at 1.2 m depth, 280 A
- **Crossing angle:** 60°
- **Soil:** standard soil

What it demonstrates
--------------------

1. Creating a :class:`~thermal_fem.CableCrossing`
2. Computing steady-state mutual temperature rise
3. Computing derating factors with
   :func:`~thermal_fem.crossing.crossing_derating_factor`
4. Tabulating the transient temperature rise over time
5. Running isolated cable simulations and superimposing the crossing effect
6. Quantifying the temperature difference: isolated vs. with crossing

Output files
------------

- ``crossing_cables_temperatures.png`` — conductor temperature comparison
  (isolated vs. with crossing) for both MV and LV cables
