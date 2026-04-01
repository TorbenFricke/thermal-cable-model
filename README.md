# Thermal FEM - Cable Temperature Rating Model

Cable temperature rating model for buried LV and MV power cables.

**Documentation:** [https://torbenfricke.github.io/thermal-cable-model/](https://torbenfricke.github.io/thermal-cable-model/) — installation, user guide, examples, theory, and API reference.

## Features

- Dynamic load profiles (constant, cyclic, daily patterns, CSV import)
- Seasonal ground temperature via Kasuda model
- Parallel cables with mutual heating (IEC 60287 image method)
- Cable crossings at arbitrary angles with derating
- Multi-layer cable models with factory methods
- Temperature outputs: conductor, insulation, cable surface, near-cable soil

## Installation

    pip install numpy scipy matplotlib

## Examples

- examples/single_cable.py: Single MV cable, steady-state and 1-year transient
- examples/parallel_cables.py: Three parallel cables with different load profiles
- examples/crossing_cables.py: MV over LV crossing at 60 degrees with derating
- examples/seasonal_dynamic.py: Year-long daily load pattern with seasonal variation
- examples/flat_vs_trefoil.py: Flat vs trefoil formation comparison

## Physics

- Cable thermal resistance T1-T3: IEC 60287-2-1 cylindrical layer model
- External soil resistance T4: image method for semi-infinite half-space
- Mutual heating: superposition of line heat sources with mirror images
- Thermal capacitance: IEC 60853 with Van Wormer splitting
- Transient solver: implicit Euler with Picard iteration for R(T)
- Cable crossings: line-source integration per CIGRE TB 640
- Ground temperature: Kasuda and Archenbach (1965) equation

## References

- IEC 60287-1-1: Electric cables, current rating equations
- IEC 60287-2-1: Thermal resistance calculations
- IEC 60853-2: Cyclic and emergency current rating
- CIGRE TB 640: Rating Calculations of Insulated Cables
- Kasuda and Archenbach (1965): Earth Temperature and Thermal Diffusivity
