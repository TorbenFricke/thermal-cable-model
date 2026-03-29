Cable Crossings
===============

When two cable circuits cross at different depths and at an angle, the heat
from one cable raises the temperature of the other at the crossing point.
The :class:`~thermal_fem.CableCrossing` class computes this mutual
temperature rise using the analytical line-source integration method from
CIGRE Technical Brochure 640.

Defining a crossing
-------------------

.. code-block:: python

   from thermal_fem import Cable, CableCrossing
   from thermal_fem.materials import SOIL_STANDARD

   mv_cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
   lv_cable = Cable.three_core_xlpe_cu(150, voltage_class="LV", voltage_kv=0.6)

   crossing = CableCrossing(
       cable_upper=mv_cable,
       cable_lower=lv_cable,
       depth_upper=0.9,          # m
       depth_lower=1.2,          # m
       crossing_angle_deg=60.0,  # degrees
       soil=SOIL_STANDARD,
   )

.. note::

   The crossing angle is measured between the cable axes projected onto the
   horizontal plane.  90° is a perpendicular crossing; 0° (parallel) is
   rejected — use the parallel cable model instead.

Steady-state temperature rise
------------------------------

Compute the temperature rise at each cable due to the heat emission of the
other:

.. code-block:: python

   W_mv = mv_cable.total_heat_per_length(350.0, 70.0)  # W/m
   W_lv = lv_cable.total_heat_per_length(280.0, 55.0)

   dT_at_mv = crossing.temperature_rise_at_upper(W_lv)
   dT_at_lv = crossing.temperature_rise_at_lower(W_mv)
   print(f"ΔT at MV cable: {dT_at_mv:.2f} °C")
   print(f"ΔT at LV cable: {dT_at_lv:.2f} °C")

The underlying integral is:

.. math::

   \Delta T = \frac{W}{4\pi\lambda}
   \int_{-L}^{L}
   \left[
       \frac{1}{\sqrt{s^2 \sin^2\!\alpha + \Delta h^2}}
     - \frac{1}{\sqrt{s^2 \sin^2\!\alpha + (h_1 + h_2)^2}}
   \right] ds

The second term is the image correction for the isothermal ground surface.

Transient temperature rise
--------------------------

For time-dependent analysis, the transient Green's function uses the
exponential integral E\ :sub:`1`:

.. code-block:: python

   for hours in [1, 4, 12, 24, 48]:
       dT = crossing.transient_temperature_rise_at_upper(
           W_lv, time_s=hours * 3600
       )
       print(f"t = {hours:3d} h → ΔT = {dT:.2f} °C")

Derating factor
---------------

The :func:`~thermal_fem.crossing.crossing_derating_factor` computes how
much the cable's permissible current must be reduced due to the crossing:

.. code-block:: python

   from thermal_fem.crossing import crossing_derating_factor

   df = crossing_derating_factor(mv_cable, depth=0.9, delta_T_crossing=dT_at_mv,
                                 soil=SOIL_STANDARD)
   print(f"Derating factor: {df:.3f}")
   print(f"Derated current: {350 * df:.0f} A")

The factor is computed as:

.. math::

   f = \sqrt{\frac{\Delta T_\text{max} - \Delta T_\text{crossing}}{\Delta T_\text{max}}}

where :math:`\Delta T_\text{max} = T_\text{max,conductor} - 20\,°\text{C}`.

Integration with ThermalSimulation
-----------------------------------

Cable crossings can be registered on a :class:`~thermal_fem.CableInstallation`
and are automatically applied during both steady-state and transient
simulations:

.. code-block:: python

   from thermal_fem import CableInstallation, ThermalSimulation

   inst = CableInstallation(soil=SOIL_STANDARD)
   inst.add_cable(mv_cable, x=0.0, depth=0.9, load=load_mv)
   inst.add_cable(lv_cable, x=0.0, depth=1.2, load=load_lv)
   inst.add_crossing(crossing)

   sim = ThermalSimulation(inst)
   result = sim.run_transient(dt=300, duration=24 * 3600)
