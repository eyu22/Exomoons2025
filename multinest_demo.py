import numpy as np
import pymultinest
import batman
from tqdm import tqdm
import os
import time
import threading
import glob
import time

start_time = time.time()

# Load your time, flux, and error arrays from AVG_SAP.dat
data = np.loadtxt('/pscratch/sd/e/emmayu/kipping_exomoons/KIC-3239945/photometry/planet1/AVG_SAP.dat')
time_array = data[:,0]
observed_flux = data[:,1]
flux_err = data[:,2]

# Clean old outputs
output_base = 'your_fit_'
for f in glob.glob(output_base + '*'):
    os.remove(f)

# Estimated max likelihood calls
max_calls = 50000

# Background thread for real-time progress
def monitor_progress(output_base, max_calls):
    log_file = output_base + 'log.txt'
    last_calls = 0
    pbar = tqdm(total=max_calls, ncols=80, position=0, leave=False,
            bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}')
    while not os.path.exists(output_base + 'post_equal_weights.dat'):
        if os.path.existsa(log_file):
            try:
                with open(log_file) as f:
                    lines = f.readlines()
                # Find most recent line with "Calls:"
                for line in reversed(lines):
                    if "Calls:" in line:
                        calls = int(line.split("Calls:")[1].split()[0])
                        break
                else:
                    calls = last_calls
                if calls > last_calls:
                    pbar.update(calls - last_calls)
                    last_calls = calls
            except Exception:
                pass
        time.sleep(0.5)
    pbar.update(max_calls - last_calls)
    pbar.close()

# Start the monitor thread
progress_thread = threading.Thread(target=monitor_progress, args=(output_base, max_calls))
progress_thread.start()

def prior(cube, ndim, nparams):
    cube[0] = 5 + cube[0] * 5           # P: 5-10 days (example)
    cube[1] = cube[1] * 1               # tau: 0-1 phase
    cube[2] = cube[2] * 0.2             # p (Rp/Rs)
    cube[3] = cube[3] * 1               # b: 0-1
    cube[4] = 10**(-3 + cube[4]*6)      # rho_star log-uniform 1e-3-1e3
    cube[5] = cube[5] * 1               # q1 limb darkening
    cube[6] = cube[6] * 1               # q2 limb darkening

def generate_mandel_agol_model(P, tau, p, b, rho_star, q1, q2, time_array):
    """
    Parameters:
    - P: orbital period [days]
    - tau: transit time (mid-transit time) [same units as time_array]
    - p: Rp/Rs (planet radius / star radius)
    - b: impact parameter (unitless)
    - rho_star: stellar density [g/cm^3 or solar units, converted below]
    - q1, q2: quadratic limb darkening coefficients in Kipping's parameterization
    - time_array: 1D numpy array of time stamps [days]
    
    Returns:
    - flux: model flux at each time point
    """

    # Constants
    G = 6.67430e-11       # gravitational constant (m^3 kg^-1 s^-2)
    Msun = 1.98847e30     # solar mass (kg)
    Rsun = 6.957e8        # solar radius (m)
    day = 86400.0         # seconds in a day
    
    # Convert stellar density rho_star to SI (kg/m^3) if needed
    # Here assume rho_star is in g/cm^3:
    rho_star_SI = rho_star * 1000  # 1 g/cm3 = 1000 kg/m3
    
    # Calculate semi-major axis in stellar radii (a/Rs) using Kepler's Third Law:
    # a^3 = (G * M_star * P^2) / (4 * pi^2)
    # M_star can be estimated from rho_star and Rs:
    # rho_star = M_star / (4/3 pi Rs^3)
    # => M_star = rho_star * (4/3) * pi * Rs^3
    
    M_star = rho_star_SI * (4/3) * np.pi * Rsun**3  # kg
    a_meters = ((G * M_star * (P * day)**2) / (4 * np.pi**2))**(1/3)
    a_rs = a_meters / Rsun  # semi-major axis in stellar radii
    
    # Calculate inclination from impact parameter: b = (a/Rs) * cos(i)
    incl_rad = np.arccos(b / a_rs)
    incl_deg = np.degrees(incl_rad)
    
    # Convert limb darkening parameters q1, q2 (Kipping 2013) to u1, u2
    u1 = 2 * np.sqrt(q1) * q2
    u2 = np.sqrt(q1) * (1 - 2 * q2)
    
    # Setup batman parameters
    params = batman.TransitParams()
    params.t0 = tau         # time of inferior conjunction (mid-transit)
    params.per = P          # orbital period in days
    params.rp = p           # planet radius (in stellar radii)
    params.a = a_rs         # semi-major axis (in stellar radii)
    params.inc = incl_deg   # orbital inclination in degrees
    params.ecc = 0.0        # eccentricity (assume circular)
    params.w = 90.0         # longitude of periastron (degrees)
    params.u = [u1, u2]     # quadratic limb darkening coefficients
    params.limb_dark = "quadratic"
    
    # Create model
    m = batman.TransitModel(params, time_array)
    flux = m.light_curve(params)
    
    return flux

# Estimate how many iterations you expect
max_calls = 50000  # adjust based on n_live_points, complexity
progress = tqdm(total=max_calls)

def loglike(cube, ndim, nparams):
    progress.update(1)
    P, tau, p, b, rho_star, q1, q2 = cube[:7]
    model_flux = generate_mandel_agol_model(P, tau, p, b, rho_star, q1, q2, time_array)
    chi2 = np.sum(((model_flux - observed_flux)/flux_err)**2)
    return -0.5 * chi2

pymultinest.run(loglike, prior, n_dims=7, n_params=7,
                outputfiles_basename='your_fit_',
                resume=False, verbose=True, n_live_points=50)

end_time = time.time()
elapsed = end_time - start_time  # elapsed time in seconds

# Optionally, convert to hours, minutes, seconds
hours = int(elapsed // 3600)
minutes = int((elapsed % 3600) // 60)
seconds = int(elapsed % 60)

print(f"MultiNest run completed in {hours}h {minutes}m {seconds}s.")