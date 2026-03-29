"""Example: 2-D soil temperature field using the FEM solver.

Computes and plots the steady-state temperature distribution around
two parallel MV cables using the finite element method.
"""

import matplotlib
matplotlib.use("Agg")

from thermal_fem.cable import Cable
from thermal_fem.ground import KasudaModel
from thermal_fem.materials import SOIL_STANDARD
from thermal_fem.fem.solver import FEMSolver
from thermal_fem.visualization import plot_soil_temperature_field


def main():
    cable = Cable.single_core_xlpe_cu(240, voltage_class="MV", voltage_kv=20.0)

    ground = KasudaModel(
        mean_surface_temp=10.0, annual_amplitude=12.0,
        day_of_min_surface_temp=35.0, soil_diffusivity=0.5e-6,
    )

    # Summer conditions
    I = 400.0
    T_cond_est = 70.0
    W = cable.total_heat_per_length(I, T_cond_est)
    print(f"Cable heat rate at {I:.0f} A: {W:.2f} W/m")

    fem = FEMSolver(
        domain_x=(-1.5, 1.5),
        domain_y=(0.0, 3.0),
        soil=SOIL_STANDARD,
        ground_model=ground,
        base_nx=50,
        base_ny=50,
    )

    fem.add_cable(cable, x=-0.2, depth=1.2, heat_rate=W)
    fem.add_cable(cable, x=0.2,  depth=1.2, heat_rate=W)

    t_summer = 200 * 86400
    print("Solving steady-state FEM (summer conditions)...")
    result = fem.solve_steady_state(time_s=t_summer)

    field = result.field_at(0)
    print(f"  Temperature range: {field.min():.1f} – {field.max():.1f} °C")

    fig = plot_soil_temperature_field(
        field, result.x, result.y,
        cable_positions=[(-0.2, 1.2), (0.2, 1.2)],
        show_mesh=True,
    )
    fig.savefig("fem_temperature_field.png", dpi=150)
    print("\nSaved: fem_temperature_field.png")


if __name__ == "__main__":
    main()
