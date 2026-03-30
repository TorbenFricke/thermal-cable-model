Theory and Physics
==================

This section summarises the physical models and numerical methods implemented
in Thermal Cable Model.  For full details, consult the referenced standards.

Cable thermal network (IEC 60287 / IEC 60853)
----------------------------------------------

Each cable is represented as a lumped-parameter thermal circuit with six
nodes per cable:

.. list-table::
   :header-rows: 1
   :widths: 10 20 70

   * - Node
     - Symbol
     - Physical meaning
   * - 0
     - θ\ :sub:`c`
     - Conductor temperature
   * - 1
     - θ\ :sub:`i`
     - Insulation midpoint (Van Wormer split of T1)
   * - 2
     - θ\ :sub:`sh`
     - Sheath / screen temperature (boundary between T1 and T2)
   * - 3
     - θ\ :sub:`a`
     - Armour temperature (boundary between T2 and T3)
   * - 4
     - θ\ :sub:`s`
     - Cable surface (outer jacket) temperature
   * - 5
     - θ\ :sub:`soil`
     - Near-cable soil node temperature

The steady-state conductor temperature rise above ambient is the sum of
the temperature drops across all resistances in the chain (IEC 60287-1-1):

.. math::

   \theta_c - T_\text{amb}
   = \bigl(W_c + \tfrac{1}{2}\,W_d\bigr)\,T_1
   + \bigl[W_c\,(1+\lambda_1) + W_d\bigr]\,T_2
   + \bigl[W_c\,(1+\lambda_1+\lambda_2) + W_d\bigr]\,T_3
   + n\,\bigl[W_c\,(1+\lambda_1+\lambda_2) + W_d\bigr]\,T_4

where:

- *n* — number of load-carrying conductors (cores) in the cable
- *W*\ :sub:`c` = *I*\ ² · *R*\ :sub:`ac`\ (*T*) — AC conductor loss per
  conductor [W/m]
- *W*\ :sub:`d` — dielectric loss per conductor [W/m]
- *T*\ :sub:`1`, *T*\ :sub:`2`, *T*\ :sub:`3` — internal thermal
  resistances (insulation, bedding, jacket) per cable [K·m/W]
- *T*\ :sub:`4` — external soil thermal resistance per cable [K·m/W]
- λ\ :sub:`1`, λ\ :sub:`2` — sheath and armour loss factors (ratio of
  sheath/armour losses to conductor losses, per conductor)

The factor *n* appears only on the *T*\ :sub:`4` term because all *n*
conductors share a single external thermal path.  In the internal
resistances *T*\ :sub:`1` – *T*\ :sub:`3`, each conductor has its own
parallel path, so the per-conductor and per-cable factors cancel.

In the transient model, these resistances are resolved as individual
branches in a six-node thermal circuit.  The resistance between each pair
of adjacent nodes is:

.. math::

   \begin{aligned}
   R_{0 \to 1} &= p \;\frac{T_1}{n}
       & &\text{(conductor → insulation midpoint)} \\[4pt]
   R_{1 \to 2} &= (1 - p)\;\frac{T_1}{n}
       & &\text{(insulation midpoint → sheath)} \\[4pt]
   R_{2 \to 3} &= \frac{T_2}{n}
       & &\text{(sheath → armour)} \\[4pt]
   R_{3 \to 4} &= \frac{T_3}{n}
       & &\text{(armour → surface)} \\[4pt]
   R_{4 \to 5} &= \tfrac{1}{2}\,T_4
       & &\text{(surface → soil node)} \\[4pt]
   R_{5 \to \text{amb}} &= \tfrac{1}{2}\,T_4
       & &\text{(soil node → ambient)}
   \end{aligned}

where *p* is the **Van Wormer coefficient** that splits the insulation
resistance T1 so that a thermal capacitance node can be placed at the
optimal intermediate point.

Each resistance occupies its own branch in the circuit, so the transient
response correctly captures the distinct thermal time constants of the
insulation, bedding, and jacket layers.

The thermal resistances are:

- **T1** — insulation (IEC 60287-2-1, cylindrical layer)
- **T2** — bedding between insulation screen and armour
- **T3** — outer serving / jacket
- **T4** — external soil thermal resistance (image method)

Internal thermal resistances T1–T3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a cylindrical layer with inner radius *r*\ :sub:`1` and outer radius
*r*\ :sub:`2`:

.. math::

   T = \frac{\rho_\text{th}}{2\pi} \ln\!\frac{r_2}{r_1}

where ρ\ :sub:`th` is the material's thermal resistivity [(K·m)/W].

Van Wormer coefficient
~~~~~~~~~~~~~~~~~~~~~~

The Van Wormer coefficient *p* determines where to split the insulation
resistance T1 to place the capacitance node.  For a cylindrical layer with
radius ratio *r*\ :sub:`2` / *r*\ :sub:`1`:

.. math::

   p = \frac{1}{2\,\ln(r_2 / r_1)}
     - \frac{1}{(r_2 / r_1)^2 - 1}

The value of *p* is always in the range [0, 0.5].  For thin layers
(*r*\ :sub:`2` ≈ *r*\ :sub:`1`) it approaches 0.5, placing the node at the
midpoint.  For thick insulation layers (typical of MV cables), *p* is
smaller, shifting the node closer to the conductor where the thermal
gradient is steepest.

The inner portion of the resistance (*p* · T1) carries the conductor heat
flux, while the outer portion ((1 − *p*) · T1) connects to the sheath.
This splitting ensures that the lumped-parameter model reproduces the
correct thermal time constant of the insulation layer as derived in
IEC 60853-2.

External thermal resistance T4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a cable buried at depth *L* with external diameter *D*\ :sub:`e` in a
soil with thermal resistivity ρ\ :sub:`soil`:

.. math::

   T_4 = \frac{\rho_\text{soil}}{2\pi}
   \ln\!\left(\frac{2L}{D_e} + \sqrt{\left(\frac{2L}{D_e}\right)^2 - 1}\right)

This is the IEC 60287-2-1 §2.2.7 formula for a single isolated cable in a
semi-infinite homogeneous medium with an isothermal ground surface.

Thermal capacitances (IEC 60853)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Thermal capacitance per unit length for a cylindrical layer:

.. math::

   Q = \rho_m c_p \cdot \pi (r_2^2 - r_1^2)

Capacitances are allocated to the six circuit nodes:

- **Node 0** (conductor): conductor thermal mass
- **Node 1** (insulation midpoint): insulation layer capacitance
- **Node 2** (sheath): screen / sheath layer capacitance
- **Node 3** (armour): bedding and armour layer capacitance
- **Node 4** (surface): jacket layer capacitance
- **Node 5** (soil): effective soil annulus around the cable

Mutual heating (image method)
-----------------------------

When multiple cables are installed in proximity, the temperature rise at
cable *i* due to heat emission from cable *j* is:

.. math::

   \Delta T_{4,ij} = \frac{\rho_\text{soil}}{2\pi}
   \ln\!\frac{d'_{ij}}{d_{ij}}

where:

- *d*\ :sub:`ij` is the real distance between cable centres
- *d'*\ :sub:`ij` is the distance from cable *i* to the **image** of cable *j*
  (reflected about the ground surface)

The image method assumes an isothermal ground surface and a homogeneous,
semi-infinite soil domain.

State-space formulation
-----------------------

The complete thermal network is expressed as a system of ordinary
differential equations:

.. math::

   \mathbf{C}\,\frac{d\boldsymbol{\theta}}{dt} +
   \mathbf{G}\,\boldsymbol{\theta} = \mathbf{P}(t)

where:

- **C** is the diagonal capacitance matrix [J/(m·K)]
- **G** is the conductance matrix [W/(m·K)] (sparse, coupling parallel cables
  through mutual heating)
- **P** is the forcing vector containing heat sources and ambient coupling
- **θ** is the state vector of node temperatures [°C]

Transient solver
----------------

The system is integrated in time using the **implicit (backward) Euler**
method:

.. math::

   \left(\frac{\mathbf{C}}{\Delta t} + \mathbf{G}\right)
   \boldsymbol{\theta}^{n+1}
   = \frac{\mathbf{C}}{\Delta t}\,\boldsymbol{\theta}^n
   + \mathbf{P}(t^{n+1})

Because the AC resistance (and hence the heat source) depends on conductor
temperature, the system is mildly non-linear.  This is handled by **Picard
iteration**: at each time step, the forcing vector is re-evaluated using the
latest conductor temperature estimate, and the linear system is re-solved.
Typically 2–3 iterations suffice.

Cable crossings (CIGRE TB 640)
------------------------------

When two cable routes cross at angle α and different depths, the steady-state
temperature rise at the target cable due to a line source at the crossing
cable is:

.. math::

   \Delta T = \frac{W}{4\pi\lambda}
   \int_{-L}^{L}
   \left[
       \frac{1}{\sqrt{s^2 \sin^2\!\alpha + \Delta h^2}}
     - \frac{1}{\sqrt{s^2 \sin^2\!\alpha + (h_1 + h_2)^2}}
   \right] ds

The first term is the contribution of the real source; the second term is the
image correction for the isothermal ground surface.

For transient analysis, the exponential integral E\ :sub:`1` formulation is
used:

.. math::

   \Delta T(t) = \frac{W}{4\pi\lambda}
   \int_{-L}^{L}
   \left[
       \mathrm{E}_1\!\left(\frac{r^2}{4\alpha t}\right)
     - \mathrm{E}_1\!\left(\frac{r'^2}{4\alpha t}\right)
   \right] ds

The **derating factor** gives the permissible current reduction:

.. math::

   f = \sqrt{\frac{\Delta T_\text{max} - \Delta T_\text{crossing}}
                   {\Delta T_\text{max}}}

Kasuda ground temperature model
-------------------------------

The Kasuda & Archenbach (1965) equation models the undisturbed ground
temperature as a function of depth *z* and time *t*:

.. math::

   T(z, t) = T_\text{mean} - T_\text{amp}\,
   \exp\!\left(-z\sqrt{\frac{\pi}{P\alpha}}\right)
   \cos\!\left(\frac{2\pi}{P}\left(t - t_0
   - \frac{z}{2}\sqrt{\frac{P}{\pi\alpha}}\right)\right)

Key features:

- **Exponential attenuation** of the annual temperature wave with depth
- **Phase lag** — temperature peaks arrive later at greater depth
- At sufficient depth (≈ 8–10 m for typical soils) the temperature approaches
  the annual mean

References
----------

- IEC 60287-1-1: *Electric cables — Calculation of the current rating —
  Part 1-1: Current rating equations (100% load factor) and calculation of
  losses — General*
- IEC 60287-2-1: *Electric cables — Calculation of the current rating —
  Part 2-1: Thermal resistance — Calculation of thermal resistance*
- IEC 60853-2: *Calculation of the cyclic and emergency current rating of
  cables — Part 2: Cyclic rating of cables greater than 18/30 (36) kV and
  emergency ratings for cables of all voltages*
- CIGRÉ Technical Brochure 640 (2015): *A Guide for Rating Calculations of
  Insulated Cables*
- Kasuda, T. and Archenbach, P.R. (1965): *Earth Temperature and Thermal
  Diffusivity at Selected Stations in the United States*, NBS Report 8972
