Seasonal Dynamic Simulation
===========================

**Script:** ``examples/seasonal_dynamic.py``

A full-year transient simulation combining a repeating daily load pattern
with seasonal ground temperature variation.  Shows how conductor temperature
varies with both load and ambient conditions.

Setup
-----

- **Cable:** 1 × 300 mm\ :sup:`2` Cu XLPE, 20 kV
- **Load:** 24-hour residential pattern (130–500 A), repeated 365 days
- **Ground temperature:** Kasuda model (mean 10.5 °C, amplitude 13 °C)
- **Burial depth:** 1.2 m
- **Time step:** 1 hour

What it demonstrates
--------------------

1. Using :meth:`LoadProfile.daily_pattern <thermal_cable_model.LoadProfile.daily_pattern>`
   for year-long simulations
2. Interaction between daily load cycling and seasonal ambient variation
3. Identifying the day and season when the peak conductor temperature occurs
4. Seasonal statistics (average and maximum per season)

Output files
------------

- ``seasonal_dynamic_temperatures.png`` — year-long temperature and current history
