Parallel Cables
================

**Script:** ``examples/parallel_cables.py``

Three single-core 20 kV cables in flat formation, each carrying a different
load profile.  This example shows how mutual heating affects the temperature
of parallel cables.

Setup
-----

- **Cable:** 1 × 240 mm\ :sup:`2` Cu XLPE, 20 kV (three identical cables)
- **Formation:** flat, 0.3 m centre-to-centre spacing at 1.0 m depth
- **Soil:** wet soil (ρ = 0.7 K·m/W)
- **Ground temperature:** Kasuda model (starting in July)
- **Loads:**

  - Left cable: cyclic industrial (500/100 A, 8 h period, 60% duty)
  - Centre cable: constant 450 A
  - Right cable: residential daily pattern

What it demonstrates
--------------------

1. Adding multiple cables to a :class:`~thermal_cable_model.CableInstallation`
2. Using different :class:`~thermal_cable_model.LoadProfile` types for each cable
3. Mutual heating between parallel cables
4. Comparing peak temperatures across cables
5. Cross-section visualisation

Output files
------------

- ``parallel_cables_temperatures.png`` — temperature and current history
- ``parallel_cables_cross_section.png`` — cable arrangement with peak temperatures
