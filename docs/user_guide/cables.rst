Cable Models
============

Thermal Cable Model represents cables as concentric cylindrical layers from the
conductor outward.  Each layer has a material, inner radius, and outer
radius.  The :class:`~thermal_cable_model.Cable` class computes thermal resistances
and heat losses from this geometry.

Factory methods
---------------

The easiest way to create a cable is through one of the built-in factory
methods, which generate realistic layer geometries from the conductor
cross-sectional area:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Method
     - Construction
     - Default voltage
   * - :meth:`Cable.single_core_xlpe_cu <thermal_cable_model.Cable.single_core_xlpe_cu>`
     - 1-core Cu / XLPE / PE jacket
     - MV (20 kV)
   * - :meth:`Cable.single_core_xlpe_al <thermal_cable_model.Cable.single_core_xlpe_al>`
     - 1-core Al / XLPE / PE jacket
     - MV (20 kV)
   * - :meth:`Cable.three_core_xlpe_cu <thermal_cable_model.Cable.three_core_xlpe_cu>`
     - 3-core Cu / XLPE / SWA
     - LV (0.6 kV)
   * - :meth:`Cable.three_core_pvc_cu <thermal_cable_model.Cable.three_core_pvc_cu>`
     - 3-core Cu / PVC / SWA
     - LV (0.6 kV)

Example:

.. code-block:: python

   from thermal_cable_model import Cable

   mv_cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
   lv_cable = Cable.three_core_xlpe_cu(150, voltage_class="LV", voltage_kv=0.6)

   print(mv_cable.name)            # "1×240 mm² Cu XLPE 20 kV"
   print(mv_cable.outer_diameter * 1e3)  # outer diameter in mm

Custom cable construction
-------------------------

For cable types not covered by the factory methods, build a cable layer by
layer:

.. code-block:: python

   from thermal_cable_model import Cable, CableLayer
   from thermal_cable_model.materials import COPPER, XLPE, LEAD_SHEATH, PVC_JACKET

   layers = [
       CableLayer(XLPE, inner_radius=8.74e-3, outer_radius=16.74e-3),
       CableLayer(LEAD_SHEATH, inner_radius=16.74e-3, outer_radius=18.24e-3),
       CableLayer(PVC_JACKET, inner_radius=18.24e-3, outer_radius=20.24e-3),
   ]

   cable = Cable(
       name="Custom 240 mm² cable",
       voltage_class="MV",
       n_conductors=1,
       conductor_area=240e-6,          # m²
       conductor_material=COPPER,
       layers=layers,
       ac_resistance_20c=7.56e-5,      # Ω/m at 20 °C
       temp_coeff_resistance=3.93e-3,  # 1/K (copper)
       max_conductor_temp=90.0,        # °C
       loss_factor_sheath=0.05,        # λ₁
       loss_factor_armour=0.0,         # λ₂
       dielectric_loss=0.3,            # W/m per phase
   )

Cable parameters
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Unit
     - Description
   * - ``conductor_area``
     - m\ :sup:`2`
     - Cross-sectional area of one conductor
   * - ``ac_resistance_20c``
     - Ω/m
     - AC resistance at 20 °C per unit length
   * - ``temp_coeff_resistance``
     - 1/K
     - Temperature coefficient of resistance (3.93 × 10\ :sup:`-3` for Cu)
   * - ``max_conductor_temp``
     - °C
     - Maximum continuous conductor temperature (90 °C for XLPE)
   * - ``loss_factor_sheath``
     - —
     - λ\ :sub:`1` — sheath-to-conductor loss ratio
   * - ``loss_factor_armour``
     - —
     - λ\ :sub:`2` — armour-to-conductor loss ratio
   * - ``dielectric_loss``
     - W/m
     - Dielectric loss per phase per unit length

Electrical properties
---------------------

The temperature-corrected AC resistance is computed as:

.. math::

   R_\text{ac}(T) = R_\text{ac,20} \bigl[1 + \alpha_{20}(T - 20)\bigr]

Joule (conductor) loss per unit length:

.. math::

   W_c = I^2 \cdot R_\text{ac}(T)

Total cable heat per unit length (including sheath, armour, and dielectric
losses):

.. math::

   W_\text{total} = n \bigl[W_c (1 + \lambda_1 + \lambda_2) + W_d\bigr]

where *n* is the number of conductors.
