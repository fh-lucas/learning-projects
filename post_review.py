import numpy as np
import matplotlib.pyplot as plt
import ansys.rocky.core as pyrocky

# =============================================================================
# CONFIGURAÇÕES DO USUÁRIO
# =============================================================================
X_LOCATION = 1.0                # posição x ao longo da placa [m]
Z_MIN = 0.005                   # altura mínima (início da camada limite)
Z_MAX = 0.025                   # altura máxima (próximo à corrente livre)
N_POSITIONS = 20                # número de pontos em z
NUM_TIME_STEPS = 10             # últimos passos de saída para média temporal

inlet_config = {
    'velocity': 0.01,           # velocidade de corrente livre [m/s]
    # ... outros parâmetros do inlet (não usados diretamente aqui)
}

solver_config = {
    'simulation_duration': 300,
    'output_frequency': 0.5,
}

# Tamanho do cubo de amostragem (deve ser maior que o diâmetro da partícula SPH = 0.0025 m)
CUBE_SIZE = (0.1, 0.01, 0.005)   # (dx, dy, dz) em metros

# Nome correto da grid function da velocidade X (verifique no Rocky)
VELOCITY_FUNC_NAME = 'Velocity X'   # ou 'Velocity : X' / 'X Velocity'

# =============================================================================
# CONEXÃO COM O ROCKY
# =============================================================================
rocky = pyrocky.launch_rocky(server_port=50361, rocky_version=261,
                             close_existing=True, headless=False)
project = rocky.api.OpenProject(r"C:\\Ansys\\Mine\\FlatPlate\\iisph_cubic_GREAT.rocky")
study = rocky.api.GetStudy()

user_process_collection = project.GetUserProcessCollection()
sph_settings = study.GetSphSettings()
materials = study.GetMaterialCollection()
fluid_material = materials.GetDefaultFluidMaterial()
density = fluid_material.GetDensity()
viscosity = fluid_material.GetViscosity()
nu = viscosity / density   # viscosidade cinemática [m²/s]

# =============================================================================
# FUNÇÕES
# =============================================================================
def blasius_profile(x, u_inf, nu, z_min, z_max, n_points):
    """Perfil de velocidade de Blasius: u(z) para placa plana."""
    z = np.linspace(z_min, z_max, n_points)
    eta_table = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0]
    fprime_table = [0, 0.0332, 0.0664, 0.0996, 0.1328, 0.1659, 0.1989, 0.2647, 0.3298, 0.3938,
                    0.4563, 0.5168, 0.5748, 0.6298, 0.6813, 0.7290, 0.7725, 0.8115, 0.8460,
                    0.9130, 0.9555, 0.9795, 0.9915]
    factor = np.sqrt(u_inf / (nu * x))
    u = np.zeros_like(z)
    for i, zi in enumerate(z):
        eta = zi * factor
        u[i] = u_inf * np.interp(eta, eta_table, fprime_table)
    return z, u

def rocky_average_velocity(user_proc_collection, x_center, z_positions,
                           num_time_steps, vel_func_name, cube_size):
    """
    Retorna (z_positions, u_mean) com média temporal ignorando NaN.
    """
    # Criar cubo
    cube = user_proc_collection.CreateCubeProcess()
    cube.SetSize(*cube_size)
    
    # Determinar índices dos passos de tempo
    total_time = solver_config['simulation_duration']
    out_freq = solver_config['output_frequency']
    n_steps_total = int(round(total_time / out_freq)) + 1
    start_step = max(0, n_steps_total - num_time_steps)
    time_steps = list(range(start_step, start_step + num_time_steps))
    
    u_mean = []
    for z in z_positions:
        cube.SetCenter(x_center, 0.0, float(z))
        vel_samples = []
        for t in time_steps:
            try:
                vel = cube.GetGridFunction(vel_func_name).GetAverage(time_step=t)
                if vel is not None and not np.isnan(vel):
                    vel_samples.append(vel)
            except Exception:
                # Se a chamada falhar, apenas ignora (pode ser que não haja partículas no cubo)
                continue
        if len(vel_samples) == 0:
            print(f"Aviso: nenhum valor válido para z = {z:.4f} m. Usando NaN.")
            u_mean.append(np.nan)
        else:
            u_mean.append(np.nanmean(vel_samples))
    return np.array(z_positions), np.array(u_mean)

def compare_and_plot(z_theory, u_theory, z_sim, u_sim, tolerance=0.05):
    """Comparação com tolerância relativa e gráfico."""
    # Interpolar o perfil teórico nos mesmos z da simulação (para comparação direta)
    u_theory_interp = np.interp(z_sim, z_theory, u_theory)
    
    # Ignorar NaNs na simulação
    valid = ~np.isnan(u_sim)
    if not np.any(valid):
        print("Nenhum dado válido na simulação para comparar.")
        return False
    
    u_sim_valid = u_sim[valid]
    u_theory_valid = u_theory_interp[valid]
    
    max_rel_err = np.max(np.abs((u_sim_valid - u_theory_valid) / u_theory_valid))
    within_tol = np.allclose(u_sim_valid, u_theory_valid, rtol=tolerance, atol=1e-8)
    
    print(f"Máximo erro relativo: {max_rel_err:.3%}")
    if within_tol:
        print("✅ Perfil dentro da tolerância.")
    else:
        print(f"❌ Perfil fora da tolerância ({tolerance:.0%}).")
    
    # Gráfico
    plt.figure(figsize=(10,6))
    plt.plot(u_sim, z_sim, 'o-', color='purple', label='Rocky (média temporal)')
    plt.plot(u_theory, z_theory, '--', color='orange', label='Blasius')
    plt.xlabel('Velocidade X [m/s]')
    plt.ylabel('Posição normal z [m]')
    plt.title(f'Perfil de velocidade em x = {X_LOCATION} m')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    return within_tol

# =============================================================================
# EXECUÇÃO PRINCIPAL
# =============================================================================
# 1. Perfil teórico
z_theory, u_theory = blasius_profile(X_LOCATION, inlet_config['velocity'], nu,
                                     Z_MIN, Z_MAX, N_POSITIONS)

# 2. Dados da simulação
z_sim = np.linspace(Z_MIN, Z_MAX, N_POSITIONS)
_, u_sim = rocky_average_velocity(user_process_collection, X_LOCATION, z_sim,
                                  NUM_TIME_STEPS, VELOCITY_FUNC_NAME, CUBE_SIZE)

# 3. Comparação
compare_and_plot(z_theory, u_theory, z_sim, u_sim, tolerance=0.05)

# 4. Finalização
project.SaveProject(False)
rocky.close()