def shear_stress_pump(gamma_dot, eta):
    """
    Calculate shear stress in the pump
    τ = γ_dot * η

    :param gamma_dot: Shear rate (s^-1)
    :param eta: Apparent viscosity (Pa s)
    :return: Shear stress (Pa)
    """
    return gamma_dot * eta

def reynolds_number(rho, u, L, eta, v):
    """
    Calculate Reynolds number
    Re = (ρuL)/η or Re = uL/v

    :param rho: Density of fluid (kg/m^3)
    :param u: Flow speed (m/s)
    :param L: Characteristic length (m)
    :param eta: Viscosity of fluid (Pa s)
    :param v: Kinematic viscosity of fluid (m^2/s)
    :return: Reynolds number
    """
    return (rho * u * L) / eta if eta else u * L / v

def blasius_shear_stress(Cf, rho, u_bar):
    """
    Calculate shear stress using Blasius equation for laminar flow
    τ = Cf * 1/2 * ρ * u_bar^2

    :param Cf: Fanning friction factor
    :param rho: Density of fluid (kg/m^3)
    :param u_bar: Average flow velocity (m/s)
    :return: Shear stress (Pa)
    """
    return Cf * 0.5 * rho * u_bar**2

def shear_stress_impeller(eta, N, Re_L):
    """
    Calculate shear stress adjacent to rotating impeller
    τ = 6.30 * η * N * Re_L^0.5

    :param eta: Viscosity of fluid (Pa s)
    :param N: Rotational speed of impeller
    :param Re_L: Local Reynolds number
    :return: Shear stress (Pa)
    """
    return 6.30 * eta * N * Re_L**0.5

def local_reynolds_number(N, d_L, rho, eta):
    """
    Calculate the local Reynolds number
    Re_L = (Nd_L^2 * ρ) / η

    :param N: Rotational speed of impeller
    :param d_L: Local diameter of the impeller (m)
    :param rho: Density of fluid (kg/m^3)
    :param eta: Viscosity of fluid (Pa s)
    :return: Local Reynolds number
    """
    return (N * d_L**2 * rho) / eta

def estimate_damage(Nr, N0, phi, n):
    """
    Estimate the magnitude of damage to cells
    N_r/N_0 = (1 - φ)^n

    :param Nr: Cell concentration remaining viable post exposure
    :param N0: Cell concentration before exposure
    :param phi: Proportion of cells that pass through high shear stress zone
    :param n: Number of passages through high shear stress zone
    :return: Ratio of viable cell concentration post exposure to pre exposure
    """
    return (1 - phi)**n

# Example usage
# shear_stress = shear_stress_pump(gamma_dot=1.0, eta=0.001)
# re = reynolds_number(rho=1000, u=1, L=1, eta=0.001, v=None)
# tau_blasius = blasius_shear_stress(Cf=0.0791, rho=1000, u_bar=1)
# tau_impeller = shear_stress_impeller(eta=0.001, N=100, Re_L=10000)
# re_local = local_reynolds_number(N=100, d_L=0.1, rho=1000, eta=0.001)
# damage_estimate = estimate_damage(Nr=5000, N0=10000, phi=0.1, n=3)
