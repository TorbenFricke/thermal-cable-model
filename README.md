# Thermal FEM - Cable Temperature Rating Model

Cable temperature rating model for buried LV and MV power cables.

## Features

- Dynamic load profiles (constant, cyclic, daily patterns, CSV import)
- Seasonal ground temperature via Kasuda model
- Parallel cables with mutual heating (IEC 60287 image method)
- Cable crossings at arbitrary angles with derating
- Multi-layer cable models with factory methods
- Temperature outputs: conductor, insulation, cable surface, near-cable soil
- 2-D FEM solver for soil temperature field visualization

## Installation

    pip install numpy scipy matplotlib

## Examples

- examples/single_cable.py: Single MV cable, steady-state and 1-year transient
- examples/parallel_cables.py: Three parallel cables with different load profiles
- examples/crossing_cables.py: MV over LV crossing at 60 degrees with derating
- examples/seasonal_dynamic.py: Year-long daily load pattern with seasonal variation
- examples/fem_temperature_field.py: 2-D soil temperature field around two cables

## Physics

- Cable thermal resistance T1-T3: IEC 60287-2-1 cylindrical layer model
- External soil resistance T4: image method for semi-infinite half-space
- Mutual heating: superposition of line heat sources with mirror images
- Thermal capacitance: IEC 60853 with Van Wormer splitting
- Transient solver: implicit Euler with Picard iteration for R(T)
- Cable crossings: line-source integration per CIGRE TB 640
- Ground temperature: Kasuda and Archenbach (1965) equation
- 2-D soil field: Q4 bilinear FEM on graded Cartesian mesh

## References

- IEC 60287-1-1: Electric cables, current rating equations
- IEC 60287-2-1: Thermal resistance calculations
- IEC 60853-2: Cyclic and emergency current rating
- CIGRE TB 640: Rating Calculations of Insulated Cables
- Kasuda and Archenbach (1965): Earth Temperature and Thermal Diffusivity
