import numpy as np
import matplotlib.pyplot as plt
import ansys.rocky.core as pyrocky
import pytest

rocky = pyrocky.launch_rocky(server_port=50361,rocky_version=261,close_existing=True, headless=False)

project = rocky.api.OpenProject(r"C:\\Ansys\\Mine\\FlatPlate\\iisph_cubic_GREAT.rocky")
study = rocky.api.GetStudy()
inlets_outlets = study.GetInletsOutletsCollection()
user_process_collection = project.GetUserProcessCollection()
sph_settings = study.GetSphSettings()
materials = study.GetMaterialCollection()
fluid_material = materials.GetDefaultFluidMaterial()
density = fluid_material.GetDensity()
viscosity = fluid_material.GetViscosity()
distance = 0.95

inlet_config = {
    'inlet_center': (0, 0.0075, 0.1),
    'inlet_length': 0.015,
    'inlet_width': 0.2,
    'inlet_orientation': (0, 0, 90),
    'boundary_condition': 'velocity',
    'velocity': 0.01,
}
solver_config = {
    'simulation_duration': 300,
    'output_frequency': 0.5,
    'simulation_target': 'GPU'
}

'''Post-processing'''

def blasius(z_min=0.005, z_max=0.025, n_positions=20):
    z_positions = np.linspace(z_min, z_max, n_positions)

    n_values = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0]
    
    derivative_values = [0, 0.0332, 0.0664, 0.0996, 0.1328, 0.1659, 0.1989, 0.2647, 0.3298, 0.3938,
                        0.4563, 0.5168, 0.5748, 0.6298, 0.6813, 0.7290, 0.7725, 0.8115, 0.8460,
                        0.9130, 0.9555, 0.9795, 0.9915]
    
    derivatives = []
    for z in z_positions:
        n = z * np.sqrt(inlet_config['velocity'] / ((viscosity/density) * distance))
        derivative_interpolation = np.interp(n, n_values, derivative_values)
        derivatives.append(derivative_interpolation)

    blasius_velocities = [inlet_config['velocity'] * derivative for derivative in derivatives]
    return np.array(z_positions), np.array(blasius_velocities)

def rocky_postprocessing(user_process_collection, sph_settings, num_time_steps=10, z_min=0.005, z_max=0.025, n_positions=20):
    heights = np.linspace(z_min, z_max, n_positions)
    cube = user_process_collection.CreateCubeProcess(sph_settings)
    cube.SetSize(0.1, 0, 0.002)

    output_interval = solver_config['output_frequency']
    total_duration = solver_config['simulation_duration']
    n_output_steps = int(round(total_duration / output_interval)) + 1
    start_step = max(0, n_output_steps - num_time_steps)
    time_steps = list(range(start_step, start_step + num_time_steps))

    rocky_velocities = []
    for z in heights:
        cube.SetCenter(1.0, 0.0, float(z))
        cube_velocity_samples = []
        for time_step in time_steps:
            velocity = cube.GetGridFunction('Velocity : X').GetAverage(time_step=time_step)
            if velocity is None:
                raise RuntimeError(f"No velocity data at time step {time_step} for cube at z={z:.4f}")
            cube_velocity_samples.append(velocity)
        rocky_velocities.append(np.mean(cube_velocity_samples))

    return heights, np.array(rocky_velocities, dtype=float)

def assert_results(blasius_velocities, rocky_velocities, tolerance):
    try:
        assert rocky_velocities == pytest.approx(blasius_velocities, rel=tolerance)
        return True
    except AssertionError:
        return False

def plot(z_positions, blasius_velocities, heights, rocky_velocities):
    plt.figure(figsize=(10, 6))
    plt.plot(rocky_velocities, heights, color='purple', linewidth=2, label='Rocky time-average')
    plt.plot(blasius_velocities, z_positions, linestyle='--', color='orange', linewidth=2, label='Blasius')
    plt.xlabel('Velocity X [m/s]')
    plt.ylabel('Wall-normal position z [m]')
    plt.title('Time-averaged X-velocity profile vs. Blasius solution')
    plt.legend()
    plt.grid(True)
    plt.show()

z_positions, blasius_velocities = blasius(z_min=0.005, z_max=0.025, n_positions=20)
heights, rocky_velocities = rocky_postprocessing(user_process_collection, sph_settings, num_time_steps=10, z_min=0.005, z_max=0.025, n_positions=20)
assert_results(blasius_velocities, rocky_velocities, tolerance=0.05)
plot(z_positions, blasius_velocities, heights, rocky_velocities)

project.SaveProject(False)
rocky.close()