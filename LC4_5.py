"""
Eclipsing Binary Stars Light Curve Simulation

Major Assumptions:
- Uniform Luminosity Distribution (No limb darkening, luminosity/projected area = const)
- Circular Orbit (e = 0)
- Celestial bodies m1 and m2 is perfectly spherical with projected areas being perfect circles
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle

# General Scientific Constants (DO NOT CHANGE)
AU = 1.5e11  # Astronomical Unit
G = 6.67e-11  # Gravitational Constant

# Astrophysics (DO NOT CHANGE)
M_Sol = 1.99e30  # Solar Mass
L_Sol = 3.9e26  # Solar Luminosity
R_Sol = 6.96e8  # Solar Radii
yr = 365 * 24 * 60 * 60  # Years in seconds

# INPUT PARAMETERS: Primary Star
target = "Algol AB (Beta Persei)"  # Name of system
m1 = 3.17  # Solar Masses
r1 = 2.73  # Solar Radii
L1 = 182 # Solar Luminosity
primary_color = 'deepskyblue'  # Color for primary star in graph

# INPUT PARAMETERS: Secondary Star
m2 = 0.7  # Solar Masses
r2 = 3.48  # r2 < r1 Solar Radii
L2 = 6.92  # Solar Luminosity
secondary_color = 'orange'  # Color for secondary star in graph

# INPUT PARAMETERS: Orbital Parameters
P = 2.867328  # Orbital Period in Days
i = 98.7 # Inclination Angle in Deg
r_L1_color = "red"  # Color for L1 point in graph


# --------------------------------------------------
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
i_min2 = 180 - i_min
i_grazing = np.arccos(((r1 - r2)*R_Sol)/sma) * (180/np.pi)  # Minimum Inclination for grazing eclipse
i_grazing2 = 180 - i_grazing


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
    if np.isscalar(d_val):
        if d_val >= (r1 + r2):
            return 0.0      
        if d_val <= abs(r1 - r2):
            return np.pi * min(r1, r2)**2
        arg1 = np.clip((d_val**2 + r1**2 - r2**2)/(2*d_val*r1), -1.0, 1.0)
        arg2 = np.clip((d_val**2 + r2**2 - r1**2)/(2*d_val*r2), -1.0, 1.0)
        term1 = r1**2 * np.arccos(arg1)
        term2 = r2**2 * np.arccos(arg2)
        term3 = 0.5 * np.sqrt(max(0.0, (d_val**2 - (r2 - r1)**2) * ((r1 + r2)**2 - d_val**2)))
        return term1 + term2 - term3
    d_val = np.asarray(d_val, dtype=float)
    area = np.zeros_like(d_val)
    outside = d_val >= (r1 + r2)
    inside = d_val <= abs(r1 - r2)
    overlap = (~outside) & (~inside)
    area[inside] = np.pi * min(r1, r2)**2
    if np.any(overlap):
        d_ov = d_val[overlap]
        arg1 = np.clip((d_ov**2 + r1**2 - r2**2)/(2*d_ov*r1), -1.0, 1.0)
        arg2 = np.clip((d_ov**2 + r2**2 - r1**2)/(2*d_ov*r2), -1.0, 1.0)
        term1 = r1**2 * np.arccos(arg1)
        term2 = r2**2 * np.arccos(arg2)
        term3 = 0.5 * np.sqrt(np.clip((d_ov**2 - (r2 - r1)**2) * ((r1 + r2)**2 - d_ov**2), 0.0, None))
        area[overlap] = term1 + term2 - term3
    return area

# Times (DO NOT CHANGE)
t1 = 0  # First Contact

# Primary Eclipse Functions (DO NOT CHANGE)
def L_PE1(t):
    """Primary Eclipse luminosity (partial)"""
    return L2 + ((A1 - A_c(t))/A1) * L1

def L_PE2(t):
    """Primary Eclipse luminosity (total)"""
    if np.isscalar(t):
        if r1 >= r2:
            return L2 + ((A1 - A2)/A1) * L1
        return L2
    else:
        if r1 >= r2:
            return np.full_like(t, L2 + ((A1 - A2)/A1) * L1)
        return np.full_like(t, L2)

# Secondary Eclipse Functions (DO NOT CHANGE)
def L_SE1(t):
    # Secondary Eclipse luminosity (partial)
    return L1 + ((A2 - A_c(t))/A2) * L2

def L_SE2(t):
    # Secondary Eclipse luminosity (total)
    if np.isscalar(t):
        if r1 >= r2:
            return L1
        return L1 + ((A2 - A1)/A2) * L2
    else:
        if r1 >= r2:
            return np.full_like(t, L1)
        return np.full_like(t, L1 + ((A2 - A1)/A2) * L2)

# Full Flux (DO NOT CHANGE)
def L_full(t):
    """Full flux"""
    if np.isscalar(t):
        return L_total
    else:
        return np.full_like(t, L_total)

# Graph Generation
fig = plt.figure(figsize=(15, 12))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=0, wspace=0.3)
ax_top = fig.add_subplot(gs[0, :])
ax_primaryeclipse = fig.add_subplot(gs[1, 0], aspect='equal')
ax_orbit = fig.add_subplot(gs[1, 1], aspect='equal')
ax_secondaryeclipse = fig.add_subplot(gs[1, 2], aspect='equal')
ax_orbit2 = fig.add_subplot(gs[1, 3], aspect='equal')
ax_top.grid(True)
ax_orbit.grid(True, alpha=0.3)
ax_primaryeclipse.grid(True, alpha=0.3)
ax_secondaryeclipse.grid(True, alpha=0.3)
ax_orbit2.grid(True, alpha=0.3)

if (r1 + r2)*R_Sol <= sma or P >= P_min:
    # Create time arrays for plotting
    t_primary = np.linspace(t1, t_total, 1000)
    t_full1 = np.linspace(t_total, P/2, 1000)
    t_secondary = np.linspace(P/2, P/2 + t_total, 1000)
    t_full2 = np.linspace(P/2 + t_total, P, 1000)

    # Plot Primary Eclipse
    ax_top.plot(t_primary, L_PE1(t_primary), 'red', label='Primary Eclipse')

    # Plot Full Flux
    ax_top.plot(t_full1, L_full(t_full1), 'black', label='Full Flux')
    ax_top.plot(t_full2, L_full(t_full2), 'black')

    # Plot Secondary Eclipse
    ax_top.plot(t_secondary, L_SE1(t_secondary), 'blue', label='Secondary Eclipse')

    # Check for total eclipse condition
    if b_h <= abs(r1 - r2):
        val = np.sqrt((((r1 - r2)**2 - b_h**2) * w**2 * R_Sol**2)/v_rel**2)
        t2 = (-np.arcsin(val) + phi)/w
        t3 = (np.arcsin(val) + phi)/w
        
        # Primary Eclipse Bottom Flux
        t_pe_total = np.linspace(t2, t3, 1000)
        ax_top.plot(t_pe_total, L_PE2(t_pe_total), 'red')
        
        # Secondary Eclipse Bottom Flux
        t_se_total = np.linspace(P/2 + t2, P/2 + t3, 1000)
        ax_top.plot(t_se_total, L_SE2(t_se_total), 'blue')


    # Determine if Roche Lobe Overflow Occurs
    def r_L1(m1, m2, sma):
        # Calculate Roche Lobe Radius for m1
        m1 = m1 * M_Sol
        m2 = m2 * M_Sol
        return  sma*np.sqrt(m1)/(np.sqrt(m1) + np.sqrt(m2)) / R_Sol
    # Plot L1 point
    r_l1 = r_L1(m1, m2, sma)
    ax_orbit.plot(r_l1, 0, 'x', color=r_L1_color, label='L1 Point', zorder=4)

    ax_orbit2.plot(-r_l1, 0, 'x', color=r_L1_color, label='L1 Point', zorder=4)

    overflow = False
    if r_l1 < r1 or (sma/R_Sol - r_l1) < r2:
        overflow = True
        print("Roche Lobe Overflow Occurs: Mass transfer likely.")

    # Set plot properties
    ax_top.set_title(f"{target}\nm₁ = {m1} M☉, r₁ = {r1} R☉, L₁ = {L1} L☉," +
                    f"    m₂ = {m2} M☉, r₂ = {r2} R☉, L₂ = {L2} L☉\n" +
                    f"P = {P/(24*60**2):.3f} d, a = {sma/AU:.4f} AU, e = 0.000, i = {i}°\nEclipse Duration: {t_total/60:.3f} min, Impact Parameter b = {b_h/r1:.3f}, Roche Lobe Overflow: {'Yes' if overflow else 'No'}")
    ax_top.set_xlabel("Seconds")
    ax_top.set_ylabel("Solar Luminosities")
    ax_top.set_xlim([0, P])

    primary_eclipse_luminosity = L_PE1(t_total/2)
    secondary_eclipse_luminosity = L_SE1(P/2 + t_total/2)
    if (primary_eclipse_luminosity <= secondary_eclipse_luminosity):
        ax_top.set_ylim([primary_eclipse_luminosity - 0.05*primary_eclipse_luminosity, 1.05*L_total])
    elif (secondary_eclipse_luminosity < primary_eclipse_luminosity):
        ax_top.set_ylim([secondary_eclipse_luminosity - 0.05*secondary_eclipse_luminosity, 1.05*L_total])
    else:
        ax_top.set_ylim([0, 1.05*L_total])
    # Plot orbital diagram
    ax_orbit.set_title(f"Full Flux 1\n{t_total/60:.2f} < t < {(P/2/60):.2f} min\nL = {L_total:.2f} L☉")
    ax_orbit.set_xlabel("Solar Radii R☉")
    ax_orbit.set_ylabel("Solar Radii R☉")
    ax_orbit.set_xlim([-1.5*sma/R_Sol, 1.5*sma/R_Sol])
    ax_orbit.set_ylim([-1.5*sma/R_Sol, 1.5*sma/R_Sol])

    # Plot circular orbit projected as ellipse
    orbit = Ellipse(xy=(0, 0), width=2*sma/R_Sol, height=2*sma/R_Sol*np.sin(np.radians(90-i)), angle=0, 
                edgecolor='black', fc='None', lw=1, label='Sky-Projected Orbit', zorder=3)
    ax_orbit.add_patch(orbit)

    # Plot stars
    primarystar = plt.Circle((0, 0), r1, color=primary_color, fill=True, label='m1', zorder=2)
    secondarystar = plt.Circle((sma/R_Sol, 0), r2, color=secondary_color, fill=True, label='m2', zorder=1)
    ax_orbit.add_patch(primarystar)
    ax_orbit.add_patch(secondarystar)

    # Plot the primary eclipse
    ax_primaryeclipse.set_title(f"Primary Eclipse\n0 < t < {(t_total)/60:.2f} min\nL = {primary_eclipse_luminosity:.2f} L☉")
    ax_primaryeclipse.set_xlabel("Solar Radii R☉")
    ax_primaryeclipse.set_ylabel("Solar Radii R☉")
    ax_primaryeclipse.set_xlim([-1.5*sma/R_Sol, 1.5*sma/R_Sol])
    ax_primaryeclipse.set_ylim([-1.5*sma/R_Sol, 1.5*sma/R_Sol])

    ax_primaryeclipse.add_patch(Ellipse(xy=(0, 0), width=2*sma/R_Sol, height=2*sma/R_Sol*np.sin(np.radians(90-i)), angle=0, 
                edgecolor='black', fc='None', lw=1, label='Sky-Projected Orbit', zorder=3))
    ax_primaryeclipse.add_patch(Circle((0, 0), r1, color=primary_color, label='m1', zorder=1))
    ax_primaryeclipse.add_patch(Circle((0, -sma/R_Sol*np.sin(np.radians(90-i))), r2, color=secondary_color, label='m2', zorder=2))

    # Plot the secondary eclipse
    ax_secondaryeclipse.set_title(f"Secondary Eclipse\n{(P/2)/60:.2f} < t < {(P/2 + t_total)/60:.2f} min\nL = {secondary_eclipse_luminosity:.2f} L☉")
    ax_secondaryeclipse.set_xlabel("Solar Radii R☉")
    ax_secondaryeclipse.set_ylabel("Solar Radii R☉")
    ax_secondaryeclipse.set_xlim([-1.5*sma/R_Sol, 1.5*sma/R_Sol])
    ax_secondaryeclipse.set_ylim([-1.5*sma/R_Sol, 1.5*sma/R_Sol])
    ax_secondaryeclipse.add_patch(Ellipse(xy=(0, 0), width=2*sma/R_Sol, height=2*sma/R_Sol*np.sin(np.radians(90-i)), angle=0, 
                edgecolor='black', fc='None', lw=1, label='Sky-Projected Orbit', zorder=3))
    ax_secondaryeclipse.add_patch(Circle((0, 0), r1, color=primary_color, label='m1', zorder=2))
    ax_secondaryeclipse.add_patch(Circle((0, -sma/R_Sol*np.sin(np.radians(90-i))), r2, color=secondary_color, label='m2', zorder=1))

    # Plot Full Flux 2
    ax_orbit2.set_title(f"Full Flux 2\n{(P/2 + t_total)/60:.2f} < t < {(P)/60:.2f} min\nL = {L_total:.2f} L☉")
    ax_orbit2.set_xlabel("Solar Radii R☉")
    ax_orbit2.set_ylabel("Solar Radii R☉")
    ax_orbit2.set_xlim([-1.5*sma/R_Sol, 1.5*sma/R_Sol])
    ax_orbit2.set_ylim([-1.5*sma/R_Sol, 1.5*sma/R_Sol])
    ax_orbit2.add_patch(Ellipse(xy=(0, 0), width=2*sma/R_Sol, height=2*sma/R_Sol*np.sin(np.radians(90-i)), angle=0, 
                edgecolor='black', fc='None', lw=1, label='Sky-Projected Orbit', zorder=3))
    ax_orbit2.add_patch(Circle((0, 0), r1, color=primary_color, label='m1', zorder=2))
    ax_orbit2.add_patch(Circle((-sma/R_Sol, 0), r2, color=secondary_color, label='m2', zorder=1))
    ax_orbit2.legend(loc = "lower center",bbox_to_anchor=(-1.5, -0.35), ncol=4)



    print(f"Orbital Period: {P/(24*60*60):.3f} days")
    print(f"Semi-major axis: {sma/AU:.4f} AU")
    print(f"Transit Duration: {t_total/60:.3f} minutes")
    print(f"Impact Parameter: {b_h:.3f} R☉,    b/r₁ = {b_h/r1:.3f}")
    print(f"Minimum Inclination for Eclipse: {i_min:.2f}° < i < {i_min2:.2f}°")
    print(f"Minimum Possible Orbital Period: {P_min/(24*60*60):.3f} <= P < {P/(24*60*60):.3f} days")
    print(f"Minimum Grazing Eclipse Inclination: {i_grazing:.2f}° < i < {i_grazing2:.2f}°")
    print(f"Roche Lobe Radius for m₁: {r_L1(m1, m2, sma):.3f} R☉")
    print("Graph saved as 'binarycurve.png'")

elif (r1 + r2) > sma/R_Sol:
    ax_top.set_title(f"Invalid System\nStars are too large for orbital radius and orbital period is too short for eclipse to occur.\nr1 + r2: {r1+r2} R☉,    SMA: {sma/R_Sol:.3f} R☉,    Minimum Possible Period: {P_min/(24*60*60):.3f} days < P = {P/(24*60*60):.3f} days\nMinimum Eclipse Inclinations: {i_min:.2f}°  < i < {i_min2:.2f}°")

if ((i_min >= i) == True or (i >= i_min2) == True):
        print("No Eclipse Occurs: Inclination too low.")
        t_full = np.linspace(0, P, 1000)
        ax_top.set_ylim([0, 1.05*L_total])
        ax_top.plot(t_full, L_full(t_full), 'black', label='Full Flux')
        ax_primaryeclipse.set_title(f"Primary Eclipse\n0 < t < {(t_total)/60:.2f} min\nL = {L_total:.2f} L☉")
        ax_secondaryeclipse.set_title(f"Secondary Eclipse\n0 < t < {(t_total)/60:.2f} min\nL = {L_total:.2f} L☉")
        
fig.savefig('binarycurve.png', dpi=500, bbox_inches='tight')
