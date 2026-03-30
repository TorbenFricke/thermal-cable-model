Materials
=========

Thermal FEM includes a built-in database of thermal material properties
sourced from IEC 60287-2-1 and standard reference tables.  All quantities
are in SI units.

The :class:`~thermal_cable_model.ThermalMaterial` class
-------------------------------------------------------

Each material is described by:

.. list-table::
   :header-rows: 1
   :widths: 35 20 45

   * - Property
     - Unit
     - Description
   * - ``thermal_conductivity``
     - W/(m·K)
     - λ — thermal conductivity
   * - ``thermal_resistivity``
     - (K·m)/W
     - ρ = 1/λ
   * - ``volumetric_heat_capacity``
     - J/(m\ :sup:`3`\ ·K)
     - ρ\ :sub:`m` · c\ :sub:`p`

The derived property :attr:`thermal_diffusivity` is computed as
α = λ / (ρ\ :sub:`m` · c\ :sub:`p`).

You can create a custom material using the
:meth:`~thermal_cable_model.ThermalMaterial.from_conductivity` class method:

.. code-block:: python

   from thermal_cable_model import ThermalMaterial

   my_soil = ThermalMaterial.from_conductivity(
       name="Sandy soil",
       conductivity=1.5,                # W/(m·K)
       volumetric_heat_capacity=2.0e6,  # J/(m³·K)
   )

Built-in materials
------------------

Conductor materials
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 20

   * - Constant
     - λ [W/(m·K)]
     - ρ [(K·m)/W]
     - ρc\ :sub:`p` [MJ/(m³·K)]
   * - ``COPPER``
     - 385
     - 0.0026
     - 3.45
   * - ``ALUMINUM``
     - 230
     - 0.0043
     - 2.42

Insulation materials (IEC 60287-2-1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 20 20 20

   * - Constant
     - λ [W/(m·K)]
     - ρ [(K·m)/W]
     - ρc\ :sub:`p` [MJ/(m³·K)]
   * - ``XLPE``
     - 0.286
     - 3.5
     - 2.4
   * - ``PVC``
     - 0.167
     - 6.0
     - 1.7
   * - ``EPR``
     - 0.286
     - 3.5
     - 2.4
   * - ``PAPER_INSULATION``
     - 0.167
     - 6.0
     - 2.0

Jacket / bedding / filler
~~~~~~~~~~~~~~~~~~~~~~~~~~

``PVC_JACKET``, ``PE_JACKET``, ``POLYPROPYLENE_FILLER``

Metallic layers
~~~~~~~~~~~~~~~

``LEAD_SHEATH`` (λ = 35 W/(m·K)), ``STEEL_ARMOUR`` (λ = 50 W/(m·K))

Soil types
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 20 20 20

   * - Constant
     - λ [W/(m·K)]
     - ρ [(K·m)/W]
     - ρc\ :sub:`p` [MJ/(m³·K)]
   * - ``SOIL_DRY``
     - 0.83
     - 1.2
     - 1.3
   * - ``SOIL_WET``
     - 1.43
     - 0.7
     - 2.5
   * - ``SOIL_STANDARD``
     - 1.0
     - 1.0
     - 1.8
   * - ``THERMAL_BACKFILL``
     - 2.0
     - 0.5
     - 2.2
   * - ``CONCRETE_DUCT``
     - 1.2
     - 0.83
     - 2.0

Other
~~~~~

``AIR_STILL`` (λ = 0.025 W/(m·K)) — for duct installations.

The :data:`~thermal_cable_model.MATERIALS` dictionary
------------------------------------------------------

All built-in materials are also available as a name-keyed dictionary:

.. code-block:: python

   from thermal_cable_model import MATERIALS

   soil = MATERIALS["Standard soil (ρ=1.0)"]
