Ground Temperature Models
=========================

The undisturbed ground temperature serves as the far-field thermal boundary
condition for the cable thermal network.  Thermal Cable Model provides two models.

Constant temperature
--------------------

The simplest boundary — a spatially and temporally uniform ambient temperature:

.. code-block:: python

   from thermal_cable_model.ground import ConstantGroundTemperature

   ground = ConstantGroundTemperature(temperature_c=15.0)

This is the default when no ground model is specified in
:class:`~thermal_cable_model.CableInstallation`.

Kasuda model
------------

The :class:`~thermal_cable_model.KasudaModel` implements the Kasuda & Archenbach
(1965) equation for the undisturbed ground temperature as a function of
depth and time of year:

.. math::

   T(z, t) = T_\text{mean} - T_\text{amp} \,
   \exp\!\left(-z\sqrt{\frac{\pi}{P\,\alpha}}\right)
   \cos\!\left(\frac{2\pi}{P}
   \left(t - t_0 - \frac{z}{2}\sqrt{\frac{P}{\pi\,\alpha}}\right)\right)

where:

- :math:`T_\text{mean}` — annual mean surface temperature [°C]
- :math:`T_\text{amp}` — half the peak-to-peak annual surface temperature swing [°C]
- :math:`P` — period (one year = 365.25 days)
- :math:`t_0` — day of minimum surface temperature
- :math:`z` — depth below surface [m]
- :math:`\alpha` — soil thermal diffusivity [m\ :sup:`2`/s]

Example:

.. code-block:: python

   from thermal_cable_model import KasudaModel

   ground = KasudaModel(
       mean_surface_temp=10.0,      # °C — central European annual mean
       annual_amplitude=12.0,       # °C — half of the annual range
       day_of_min_surface_temp=35,  # early February
       soil_diffusivity=0.5e-6,     # m²/s — typical soil
   )

   # Temperature at 1.2 m depth in mid-summer (day 200)
   T = ground.temperature(depth=1.2, time_s=200 * 86400)
   print(f"{T:.1f} °C")

The model captures two key physical effects:

1. **Amplitude attenuation** — the annual temperature swing decays
   exponentially with depth.
2. **Phase lag** — the peak temperature at depth occurs later in the year
   than at the surface.

Choosing parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 25 40

   * - Parameter
     - Typical range
     - Notes
   * - ``mean_surface_temp``
     - 5–20 °C
     - From local climate data
   * - ``annual_amplitude``
     - 5–15 °C
     - Half the difference between warmest and coldest monthly means
   * - ``day_of_min_surface_temp``
     - 20–50
     - Day of year; ≈ 35 for central Europe (early February)
   * - ``soil_diffusivity``
     - 0.3–0.8 × 10\ :sup:`-6` m\ :sup:`2`/s
     - Depends on soil type and moisture content

Custom ground models
--------------------

Implement the :class:`~thermal_cable_model.GroundTemperatureModel` abstract base
class to define your own boundary condition:

.. code-block:: python

   from thermal_cable_model.ground import GroundTemperatureModel

   class MyGroundModel(GroundTemperatureModel):
       def temperature(self, depth: float, time_s: float) -> float:
           # your implementation here
           ...
