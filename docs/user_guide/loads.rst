Load Profiles
=============

The :class:`~thermal_fem.LoadProfile` class defines the time-varying RMS
current flowing through a cable.  Internally it stores time-stamped
current values and performs piecewise-linear interpolation.

Constant load
-------------

.. code-block:: python

   from thermal_fem import LoadProfile

   load = LoadProfile.constant(current_a=400.0, duration_s=48 * 3600)

Cyclic (on/off) load
---------------------

A rectangular waveform alternating between a peak and base current:

.. code-block:: python

   load = LoadProfile.cyclic(
       peak_current=500.0,
       base_current=100.0,
       period_s=8 * 3600,    # 8-hour cycle
       duty_cycle=0.6,       # 60% at peak
       n_cycles=6,
   )

Daily pattern
-------------

Specify 24 hourly current values (one per hour, starting at 00:00) and
repeat for multiple days:

.. code-block:: python

   hourly = [
       120, 100, 90, 85, 80, 90, 150, 280,
       350, 320, 300, 290, 280, 270, 260, 270,
       310, 380, 420, 400, 350, 280, 200, 150,
   ]
   load = LoadProfile.daily_pattern(hourly, n_days=365)

CSV import
----------

Load a two-column CSV file (time in seconds, current in amperes):

.. code-block:: python

   load = LoadProfile.from_csv(
       "measured_load.csv",
       time_col=0,
       current_col=1,
       delimiter=",",
       skip_header=1,
   )

Custom profile
--------------

For full control, pass arrays directly:

.. code-block:: python

   import numpy as np

   times = np.array([0, 3600, 7200, 10800])
   currents = np.array([0, 400, 400, 200])
   load = LoadProfile(times, currents)

Querying the profile
--------------------

.. code-block:: python

   load.current_at(5000.0)       # interpolated current at t = 5000 s
   load.current_array(np.arange(0, 3600, 60))  # vectorised
   load.duration                  # total duration in seconds
