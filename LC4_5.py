"""
Eclipsing Binary Stars Light Curve Simulation

Major Assumptions:
- Uniform Luminosity Distribution (No limb darkening, luminosity/projected area = const)
- Circular Orbit (e = 0)
- Celestial bodies m1 and m2 is perfectly spherical with projected areas being perfect circles
"""

import numpy as np
import matplotlib.pyplot as plt

# General Scientific Constants (DO NOT CHANGE)
AU = 1.5e11  # Astronomical Unit
G = 6.67e-11  # Gravitational Constant

# Astrophysics (DO NOT CHANGE)
M_Sol = 1.99e30  # Solar Mass
L_Sol = 3.9e26  # Solar Luminosity
R_Sol = 6.96e8  # Solar Radii
yr = 365 * 24 * 60 * 60  # Years in seconds

# INPUT PARAMETERS
m1 = 0.5  # Solar Masses
m2 = 0.25  # Solar Masses
r1 = 0.6  # Solar Radii
r2 = 0.5  # r2 < r1 Solar Radii
L1 = 1  # Solar Luminosity
L2 = 0.5  # Solar Luminosity
P = 0.154  # Orbital Period in Days
i = 85  # Inclination Angle in Deg

# Calculations (DO NOT CHANGE)
A1 = np.pi * r1**2  # Area of Star 1 in R_Solar^2
A2 = np.pi * r2**2  # Area of Star 2 in R_Solar^2
A_total = A1 + A2  # Total Area R_Solar^2
P = P * 24 * 60 * 60
sma = ((P/yr)**2 * (m1 + m2))**(1/3) * AU  # SMA in m
v_rel = np.sqrt((G*(m1 + m2) * M_Sol)/sma)  # Relative Orbital Velocity in m/s
b_h = sma*np.cos(i * np.pi/180)/R_Sol  # Impact Parameter Physical Length in Solar Radii
w = (2 * np.pi)/P  # Angular Velocity of system in rad/s
phi = np.arcsin(np.sqrt((r1 + r2)**2 - b_h**2) * R_Sol/sma)  # Phase Shift for graph in rad
t_total = P/np.pi * np.arcsin(np.sqrt((r1 + r2)**2 - b_h**2) * R_Sol/sma)  # Transit Duration in s
A1 = np.pi * r1**2  # Projected Area of Star 1 in R_Sol^2
A2 = np.pi * r2**2  # Projected Area of Star 2 in R_Sol^2
L_total = L1 + L2  # Total Luminosity of System in L_Sol

# Boundary Limits (DO NOT CHANGE)
P_min = np.sqrt((((r1+r2)*R_Sol)/AU)**3/(m1+m2)) * yr  # Minimum Orbital Period in Seconds
i_min = np.arccos(((r1 + r2)*R_Sol)/sma) * (180/np.pi)  # Minimum Inclination for eclipse in Deg
i_grazing = np.arccos(((r1 - r2)*R_Sol)/sma) * (180/np.pi)  # Maximum Inclination for grazing eclipse

# Equations (DO NOT CHANGE)
def v(t):
    """Radial Velocity as a function of t"""
    return v_rel*np.cos(w*t - phi)

def p(t):
    """Starting Position of center of Star 2, = 0 at first contact"""
    return ((v_rel/w)/R_Sol) * np.sin(w*t - phi)

def d(t):
    """Distance function"""
    return np.sqrt(b_h**2 + p(t)**2)

# Projected Area Eclipsed in R_Sol^2 as a function of time (DO NOT CHANGE)
def A_c(t):
    """Projected Area Eclipsed"""
    d_val = d(t)
    term1 = r1**2 * np.arccos((d_val**2 + r1**2 - r2**2)/(2*d_val*r1))
    term2 = r2**2 * np.arccos((d_val**2 + r2**2 - r1**2)/(2*d_val*r2))
    term3 = 0.5 * np.sqrt((d_val**2 - (r2 - r1)**2) * ((r1 + r2)**2 - d_val**2))
    return term1 + term2 - term3

# Times (DO NOT CHANGE)
t1 = 0  # First Contact

# Primary Eclipse Functions (DO NOT CHANGE)
def L_PE1(t):
    """Primary Eclipse luminosity (partial)"""
    return L2 + ((A1 - A_c(t))/A1) * L1

def L_PE2(t):
    """Primary Eclipse luminosity (total)"""
    if np.isscalar(t):
        return L2 + ((A1 - A2)/A1) * L1
    else:
        return np.full_like(t, L2 + ((A1 - A2)/A1) * L1)

# Secondary Eclipse Functions (DO NOT CHANGE)
def L_SE1(t):
    """Secondary Eclipse luminosity (partial)"""
    return L1 + ((A2 - A_c(t))/A2) * L2

def L_SE2(t):
    """Secondary Eclipse luminosity (total)"""
    if np.isscalar(t):
        return L1
    else:
        return np.full_like(t, L1)

# Full Flux (DO NOT CHANGE)
def L_full(t):
    """Full flux"""
    if np.isscalar(t):
        return L_total
    else:
        return np.full_like(t, L_total)

# Graph Generation
plt.figure(1, figsize=(10, 6))
plt.grid(True)

# Create time arrays for plotting
t_primary = np.linspace(t1, t_total, 1000)
t_full1 = np.linspace(t_total, P/2, 1000)
t_secondary = np.linspace(P/2, P/2 + t_total, 1000)
t_full2 = np.linspace(P/2 + t_total, P, 1000)

# Plot Primary Eclipse
plt.plot(t_primary, L_PE1(t_primary), 'red', label='Primary Eclipse')

# Plot Full Flux
plt.plot(t_full1, L_full(t_full1), 'black', label='Full Flux')
plt.plot(t_full2, L_full(t_full2), 'black')

# Plot Secondary Eclipse
plt.plot(t_secondary, L_SE1(t_secondary), 'blue', label='Secondary Eclipse')

# Check for total eclipse condition
if b_h <= (r1 - r2):
    val = np.sqrt((((r1 - r2)**2 - b_h**2) * w**2 * R_Sol**2)/v_rel**2)
    t2 = (-np.arcsin(val) + phi)/w
    t3 = (np.arcsin(val) + phi)/w
    
    # Primary Eclipse Bottom Flux
    t_pe_total = np.linspace(t2, t3, 1000)
    plt.plot(t_pe_total, L_PE2(t_pe_total), 'red')
    
    # Secondary Eclipse Bottom Flux
    t_se_total = np.linspace(P/2 + t2, P/2 + t3, 1000)
    plt.plot(t_se_total, L_SE2(t_se_total), 'blue')

# Set plot properties
plt.title(f"m₁ = {m1} M☉, r₁ = {r1} R☉, L₁ = {L1} L☉\n" +
          f"m₂ = {m2} M☉, r₂ = {r2} R☉, L₂ = {L2} L☉\n" +
          f"P = {P/(24*60**2):.3f} d, a = {sma/AU:.4f} AU, i = {i}°")
plt.xlabel("Seconds")
plt.ylabel("Solar Luminosities")
plt.xlim([0, P])
plt.ylim([0, 1.2*L_total])
plt.legend()

# Save figure
plt.savefig('/mnt/user-data/outputs/binarycurve.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"Simulation complete!")
print(f"Orbital Period: {P/(24*60*60):.3f} days")
print(f"Semi-major axis: {sma/AU:.4f} AU")
print(f"Transit Duration: {t_total/60:.2f} minutes")
print(f"Impact Parameter: {b_h:.3f} R☉")
print(f"Minimum Inclination: {i_min:.2f}°")
print(f"Grazing Eclipse Inclination: {i_grazing:.2f}°")
