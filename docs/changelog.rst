Changelog
=========

v0.1.0 (2025)
--------------

Initial release.

- Multi-layer cable model with factory methods for single-core and three-core
  LV/MV cables (Cu and Al conductors, XLPE and PVC insulation)
- Built-in thermal material database (IEC 60287-2-1)
- Dynamic load profiles: constant, cyclic, daily pattern, CSV import
- Kasuda ground temperature model for seasonal variation
- IEC 60287 lumped-parameter thermal network (T1–T4)
- Mutual heating for parallel cables via image method
- IEC 60853 thermal capacitance with Van Wormer splitting
- Implicit Euler transient solver with Picard iteration
- Cable crossing analysis (CIGRE TB 640) with derating factor
- 2-D finite element solver for soil temperature fields (Q4 elements)
- Plotting utilities for temperature history, cross-sections, and contour fields
- Six example scripts covering single cable, parallel cables, crossings,
  seasonal dynamics, FEM fields, and flat vs trefoil comparisons
