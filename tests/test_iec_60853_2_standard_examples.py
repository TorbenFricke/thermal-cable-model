"""IEC 60853-2 style transient checks vs `thermal_cable_model`.

IEC 60853-2 describes cyclic and emergency rating using the same lumped thermal
network philosophy as IEC 60287.  Full reproduction of **Appendix F** worked
examples (e.g. Example F6) requires the complete cable data and tabulated
thermal resistances/capacitances from a licensed copy of **IEC 60853-2**; those
tables are not embedded here.

The tests below validate behaviours that the standard’s methods assume:

- After a **step increase** in load (emergency-style), the transient solution must
  converge to the **steady-state** conductor temperature predicted by the same
  network (IEC 60287-1-1 steady-state rise at the new current).

A skipped placeholder documents where a strict Appendix F regression would go
once inputs are transcribed from the standard.
"""

from __future__ import annotations

import numpy as np
import pytest

from thermal_cable_model.cable import Cable
from thermal_cable_model.ground import ConstantGroundTemperature
from thermal_cable_model.loads import LoadProfile
from thermal_cable_model.materials import SOIL_STANDARD
from thermal_cable_model.simulation import CableInstallation, ThermalSimulation


@pytest.mark.skip(
    reason=(
        "IEC 60853-2 Appendix F Example F6: needs full cable T/Q and ambient "
        "data transcribed from licensed IEC 60853-2; compare θ_c(t) at listed "
        "times to standard tables."
    )
)
def test_iec60853_2_appendix_f_example_f6_conductor_temperature_trace():
    """Placeholder for strict Appendix F regression (see module docstring)."""
    pass  # pragma: no cover


class TestIec60853_2EmergencyStepLoad:
    """Step load transient approaches IEC steady-state at the elevated current."""

    # Soil-node thermal mass in the 6-lump network relaxes slowly; a multi-day tail
    # at I₂ is required before θ_c matches the 60287 steady-state at I₂ (see also
    # IEC 60853-2 discussion of long emergency durations).
    DT_S = 3600.0
    TOL_K = 0.05

    def test_step_up_then_converges_to_steady_state_rating(self):
        cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)
        ground = ConstantGroundTemperature(15.0)
        depth_m = 1.2
        i_before = 200.0
        i_after = 450.0
        t_step_s = 4.0 * 3600.0
        tail_s = 40.0 * 24.0 * 3600.0
        duration = t_step_s + tail_s

        times = np.array([0.0, t_step_s - 1e-3, t_step_s, duration])
        currents = np.array([i_before, i_before, i_after, i_after])
        load = LoadProfile(times, currents)

        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(cable, 0.0, depth_m, load)
        sim = ThermalSimulation(inst)
        tr = sim.run_transient(dt=self.DT_S, duration=duration)

        inst_ref = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst_ref.add_cable(
            cable,
            0.0,
            depth_m,
            LoadProfile.constant(i_after, duration),
        )
        sim_ref = ThermalSimulation(inst_ref)
        ss = sim_ref.run_steady_state()

        tc_final = float(tr.conductor_temps[-1, 0])
        tc_ss = float(ss.conductor_temps[0, 0])
        assert tc_final == pytest.approx(tc_ss, abs=self.TOL_K)

    def test_conductor_rises_after_step(self):
        """Immediately after the step, conductor temperature increases monotonically
        over a short post-step window (physical expectation for higher losses)."""
        cable = Cable.single_core_xlpe_cu(150, voltage_class="MV", voltage_kv=20.0)
        ground = ConstantGroundTemperature(18.0)
        depth_m = 1.0
        t_step_s = 6.0 * 3600.0
        duration = t_step_s + 8.0 * 3600.0
        times = np.array([0.0, t_step_s - 1e-3, t_step_s, duration])
        currents = np.array([80.0, 80.0, 300.0, 300.0])
        load = LoadProfile(times, currents)

        inst = CableInstallation(soil=SOIL_STANDARD, ground_temp_model=ground)
        inst.add_cable(cable, 0.0, depth_m, load)
        sim = ThermalSimulation(inst)
        tr = sim.run_transient(dt=self.DT_S, duration=duration)

        idx_step = int(round(t_step_s / self.DT_S))
        idx_end_short = min(idx_step + 12, tr.n_steps - 1)
        segment = tr.conductor_temps[idx_step : idx_end_short + 1, 0]
        assert np.all(np.diff(segment) >= -1e-6), (
            "expected conductor temperature not to decrease right after a load step up"
        )
