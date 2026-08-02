import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

# General Scientific Constants (DO NOT CHANGE)
AU = 1.5e11        # Astronomical Unit
G = 6.67e-11       # Gravitational Constant

# Astrophysics (DO NOT CHANGE)
M_Sol = 1.99e30    # Solar Mass
L_Sol = 3.9e26     # Solar Luminosity
R_Sol = 6.96e8     # Solar Radii
yr = 365 * 24 * 60 * 60  # Years in seconds

# INPUT PARAMETERS: Primary Star
target = "Algol AB (Beta Persei)"
m1 = 3.17           # Solar Masses
r1 = 2.73           # Solar Radii
L1 = 182            # Solar Luminosity
primary_color = 'deepskyblue'

# INPUT PARAMETERS: Secondary Star
m2 = 0.7            # Solar Masses
r2 = 3.48           # Solar Radii
L2 = 6.92           # Solar Luminosity
secondary_color = 'orange'

# INPUT PARAMETERS: Orbital Elements
ORBIT_INPUT = "P"  # a
P = 2.867328        # Orbital Period in Days (used when ORBIT_INPUT == "P")
a_AU = 0.062043     # Semi-Major Axis in AU (used when ORBIT_INPUT == "a")
i = 98.7            # Inclination Angle in Deg
e = 0.00            # Eccentricity (0 <= e < 1)
omega = 0.0        # Argument of Periastron in Deg (orientation of ellipse in orbital plane,
                     # measured from the ascending node; t=0 is defined as periastron passage)
r_L1_color = "red"

N_SAMPLES = 5 * 10**6   # Time resolution for the numerical scan (higher = smoother/sharper eclipses)
N_PERIODS = 1        # How many consecutive periods to display on the light curve plot (>= 1)

# --------------------------------------------------
# Derived Quantities (DO NOT CHANGE)
A1 = np.pi * r1**2
A2 = np.pi * r2**2
A_total = A1 + A2
L_total = L1 + L2

# P and a are made mutually dependent via Kepler's third law, given the masses above --
# whichever one is NOT chosen as ORBIT_INPUT is derived from the other, so they can never
# drift out of consistency with (m1 + m2).
if ORBIT_INPUT == "P":
    P = P * 24 * 60 * 60                          # seconds
    sma = ((P / yr)**2 * (m1 + m2))**(1 / 3) * AU  # semi-major axis, meters (derived, Kepler III)
elif ORBIT_INPUT == "a":
    sma = a_AU * AU                                # semi-major axis, meters (direct input)
    P = np.sqrt((sma / AU)**3 / (m1 + m2)) * yr    # Orbital period, seconds (derived, Kepler III)
else:
    raise ValueError("ORBIT_INPUT must be 'P' or 'a'")

a_Rsol = sma / R_Sol
w = (2 * np.pi) / P                                    # mean motion, rad/s


# --------------------------------------------------
# Kepler's Equation (DO NOT CHANGE)
def solve_kepler(M, e, tol=1e-12, max_iter=200):
    """Newton-Raphson solve of M = E - e*sin(E) for eccentric anomaly E."""
    M = np.atleast_1d(np.asarray(M, dtype=float))
    E = M.copy() if e < 0.8 else np.full_like(M, np.pi)
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        fp = 1 - e * np.cos(E)
        dE = f / fp
        E -= dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def true_anomaly(E, e):
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))


def orbit_state(t):
    """
    Relative position of star 2 w.r.t. star 1, in Solar Radii.
    Returns: x, y (sky-plane), z (line-of-sight, +z = star 2 in front of star 1),
             d (projected separation), r (true 3D separation), nu (true anomaly)
    """
    M = np.mod(w * np.asarray(t, dtype=float), 2 * np.pi)
    E = solve_kepler(M, e)
    nu = true_anomaly(E, e)
    r = a_Rsol * (1 - e**2) / (1 + e * np.cos(nu))
    om = np.radians(omega)
    inc = np.radians(i)
    x = r * np.cos(om + nu)
    y = r * np.sin(om + nu) * np.cos(inc)
    z = r * np.sin(om + nu) * np.sin(inc)
    d = np.sqrt(x**2 + y**2)
    return x, y, z, d, r, nu


# --------------------------------------------------
# Projected Area Eclipsed in R_Sol^2 (DO NOT CHANGE)
def A_c(d):
    d = np.atleast_1d(np.asarray(d, dtype=float))
    area = np.zeros_like(d)
    outside = d >= (r1 + r2)
    inside = d <= abs(r1 - r2)
    overlap = (~outside) & (~inside)
    area[inside] = np.pi * min(r1, r2)**2
    if np.any(overlap):
        dv = d[overlap]
        arg1 = np.clip((dv**2 + r1**2 - r2**2) / (2 * dv * r1), -1.0, 1.0)
        arg2 = np.clip((dv**2 + r2**2 - r1**2) / (2 * dv * r2), -1.0, 1.0)
        term1 = r1**2 * np.arccos(arg1)
        term2 = r2**2 * np.arccos(arg2)
        term3 = 0.5 * np.sqrt(np.clip((dv**2 - (r2 - r1)**2) * ((r1 + r2)**2 - dv**2), 0.0, None))
        area[overlap] = term1 + term2 - term3
    return area


def flux(d, z):
    """Total system luminosity given projected separation d and line-of-sight offset z."""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    z = np.atleast_1d(np.asarray(z, dtype=float))
    L = np.full_like(d, L_total)
    eclipsing = d < (r1 + r2)
    if np.any(eclipsing):
        dc = d[eclipsing]
        zc = z[eclipsing]
        Ac = A_c(dc)
        Lc = np.empty_like(dc)
        total = dc <= abs(r1 - r2)
        front2 = zc > 0   # star 2 in front -> occults star 1  ("primary" eclipse, if star1 brighter)
        front1 = ~front2  # star 1 in front -> occults star 2  ("secondary" eclipse)

        m = front2 & total
        Lc[m] = L2 if r2 >= r1 else L2 + (A1 - A2) / A1 * L1
        m = front2 & ~total
        Lc[m] = L2 + (A1 - Ac[m]) / A1 * L1

        m = front1 & total
        Lc[m] = L1 if r1 >= r2 else L1 + (A2 - A1) / A2 * L2
        m = front1 & ~total
        Lc[m] = L1 + (A2 - Ac[m]) / A2 * L2

        L[eclipsing] = Lc
    return L


def get_segments(mask):
    """Indices of contiguous True runs in a boolean array."""
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return []
    splits = np.where(np.diff(idx) > 1)[0]
    return np.split(idx, splits + 1)


def get_segments_circular(mask):
    """
    Like get_segments, but treats the array as periodic: if a True run touches both
    the first and last sample, it's really one continuous run that straddles the
    t=0/P seam, so the fragment at the end is stitched onto the fragment at the start
    (in chronological order: end-of-array piece, then start-of-array piece).
    Needed because internally t=0 is periastron, and an eclipse can happen to occur
    right around periastron -- without this, that eclipse would be seen as two
    unrelated fragments and only one would be picked up.
    """
    segs = [list(s) for s in get_segments(mask)]
    if len(segs) > 1 and mask[0] and mask[-1]:
        first, last = segs[0], segs[-1]
        segs = segs[1:-1] + [last + first]
    return segs


# --------------------------------------------------
# Numerical scan over one full period
t_arr = np.linspace(0, P, N_SAMPLES)
x_arr, y_arr, z_arr, d_arr, r_arr, nu_arr = orbit_state(t_arr)
L_arr = flux(d_arr, z_arr)

eclipse_mask = d_arr < (r1 + r2)
primary_mask = eclipse_mask & (z_arr > 0)     # star 2 occults star 1
secondary_mask = eclipse_mask & (z_arr < 0)   # star 1 occults star 2

primary_segs = get_segments_circular(primary_mask)
secondary_segs = get_segments_circular(secondary_mask)

eclipses_occur = len(primary_segs) > 0 or len(secondary_segs) > 0

# --------------------------------------------------
# Phase shift: re-anchor t=0 to the start (first contact) of the primary eclipse.
# t_arr up to this point is referenced to periastron passage (t=0); since the sample
# grid is uniform, "shifting the zero point" is just a roll of every sampled array by
# the number of samples between periastron and primary-eclipse first contact.
# --------------------------------------------------
# Phase shift: re-anchor t=0 to the start (first contact) of the primary eclipse.
# t_arr up to this point is referenced to periastron passage (t=0); since the sample
# grid is uniform, "shifting the zero point" is just a roll of every sampled array by
# the number of samples between periastron and primary-eclipse first contact.
#
# "Primary" is defined by DEPTH (lower L_min), not by which star happens to be in
# front -- so we compare the two candidate eclipses' minimum flux here, before
# rolling, rather than assuming z>0 (star 2 in front) is always the primary.
def _deepest_L(segs):
    if not segs:
        return None
    return L_arr[max(segs, key=len)].min()

L_min_primary_raw = _deepest_L(primary_segs)
L_min_secondary_raw = _deepest_L(secondary_segs)

if primary_segs and secondary_segs:
    idx0 = max(primary_segs, key=len)[0] if L_min_primary_raw <= L_min_secondary_raw \
        else max(secondary_segs, key=len)[0]
elif primary_segs:
    idx0 = max(primary_segs, key=len)[0]
elif secondary_segs:
    # No primary eclipse (grazing/geometry-limited case) -> anchor to secondary instead
    idx0 = max(secondary_segs, key=len)[0]
else:
    idx0 = 0

x_arr = np.roll(x_arr, -idx0)
y_arr = np.roll(y_arr, -idx0)
z_arr = np.roll(z_arr, -idx0)
d_arr = np.roll(d_arr, -idx0)
r_arr = np.roll(r_arr, -idx0)
nu_arr = np.roll(nu_arr, -idx0)
L_arr = np.roll(L_arr, -idx0)

eclipse_mask = d_arr < (r1 + r2)
primary_mask = eclipse_mask & (z_arr > 0)
secondary_mask = eclipse_mask & (z_arr < 0)
# Use the circular segmenter here too: the roll only guarantees the ANCHORED eclipse
# starts cleanly at the new index 0. The other eclipse (or, in extreme/high-e geometries,
# a secondary grazing crossing) could in principle still straddle the new t=0/P seam and
# get artificially split into two fragments if we used the plain (non-circular) segmenter.
primary_segs = get_segments_circular(primary_mask)
secondary_segs = get_segments_circular(secondary_mask)
# t_arr itself is unchanged (0..P, uniform) -- it now reads as "time since primary
# eclipse first contact" rather than "time since periastron".

r_peri = a_Rsol * (1 - e)  # periastron separation, R_Sol (used in summary printout below)

# --------------------------------------------------
# Graph Generation
fig = plt.figure(figsize=(15, 12))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=0, wspace=0.3)
ax_top = fig.add_subplot(gs[0, :])
ax_p1 = fig.add_subplot(gs[1, 0], aspect='equal')
ax_p2 = fig.add_subplot(gs[1, 1], aspect='equal')
ax_p3 = fig.add_subplot(gs[1, 2], aspect='equal')
ax_p4 = fig.add_subplot(gs[1, 3], aspect='equal')
for ax in (ax_top, ax_p1, ax_p2, ax_p3, ax_p4):
    ax.grid(True, alpha=0.3 if ax is not ax_top else 1.0)

if not eclipses_occur:
    # No eclipse: flat full-flux light curve, four identical "full system" panels
    for k in range(N_PERIODS):
        ax_top.plot(t_arr + k * P, L_arr, 'black', label='Full Flux' if k == 0 else None)
    ax_top.set_ylim([0, 1.05 * L_total])
    print("No Eclipse Occurs: inclination/eccentricity/geometry do not produce a transit.")
    panel_axes = [ax_p1, ax_p2, ax_p3, ax_p4]
    panel_times = [0, P / 4, P / 2, 3 * P / 4]
    panel_titles = ["Full System"] * 4
    pe = se = None
    nu_c_pe = nu_c_se = r_c_pe = r_c_se = None
else:
    # --- Identify primary & secondary eclipse windows (largest segment of each, if present) ---
    def eclipse_info(segs):
        if not segs:
            return None
        seg = max(segs, key=len)
        mid_idx = seg[np.argmin(L_arr[seg])]
        return {
            't_start': t_arr[seg[0]], 't_end': t_arr[seg[-1]],
            't_mid': t_arr[mid_idx], 'L_min': L_arr[mid_idx],
            'd_min': d_arr[mid_idx], 'duration': t_arr[seg[-1]] - t_arr[seg[0]],
        }

    pe_z_pos = eclipse_info(primary_segs)      # star 2 in front (occults star 1)
    se_z_neg = eclipse_info(secondary_segs)    # star 1 in front (occults star 2)

    # True anomaly and instantaneous separation AT each conjunction (x=0, i.e. om+nu=90 or
    # 270 deg) -- needed to properly generalize the analytic eclipse-geometry formulas
    # below to eccentric orbits, where the two conjunctions occur at different separations
    # unless omega = 0 or 180.
    nu_c_zpos = np.radians(90.0 - omega)
    nu_c_zneg = np.radians(270.0 - omega)
    r_c_zpos = a_Rsol * (1 - e**2) / (1 + e * np.cos(nu_c_zpos))
    r_c_zneg = a_Rsol * (1 - e**2) / (1 + e * np.cos(nu_c_zneg))

    # The "primary" eclipse is, by convention, whichever one is DEEPER (lower L_min) --
    # not necessarily the one where a particular star is in front. If star 2 is the
    # brighter star, occulting star 1 can actually produce the shallower dip, so we
    # compare depths directly rather than assuming z>0 is always "primary". The same
    # branching also carries along each conjunction's (nu_c, r_c) for the analytic
    # eclipse-geometry formulas further down.
    if pe_z_pos and se_z_neg:
        if pe_z_pos['L_min'] <= se_z_neg['L_min']:
            pe, pe_segs = pe_z_pos, primary_segs
            se, se_segs = se_z_neg, secondary_segs
            nu_c_pe, r_c_pe = nu_c_zpos, r_c_zpos
            nu_c_se, r_c_se = nu_c_zneg, r_c_zneg
        else:
            pe, pe_segs = se_z_neg, secondary_segs
            se, se_segs = pe_z_pos, primary_segs
            nu_c_pe, r_c_pe = nu_c_zneg, r_c_zneg
            nu_c_se, r_c_se = nu_c_zpos, r_c_zpos
    elif pe_z_pos:
        pe, pe_segs = pe_z_pos, primary_segs
        se, se_segs = None, []
        nu_c_pe, r_c_pe = nu_c_zpos, r_c_zpos
        nu_c_se, r_c_se = None, None
    elif se_z_neg:
        pe, pe_segs = se_z_neg, secondary_segs
        se, se_segs = None, []
        nu_c_pe, r_c_pe = nu_c_zneg, r_c_zneg
        nu_c_se, r_c_se = None, None
    else:
        pe, se, pe_segs, se_segs = None, None, [], []
        nu_c_pe = nu_c_se = r_c_pe = r_c_se = None

    # Plot light curve: black baseline, red over the (deeper) primary eclipse, blue over
    # the (shallower) secondary eclipse -- repeated across N_PERIODS consecutive periods,
    # since the underlying t_arr/L_arr/segments are one period's worth of data and simply
    # recur every P seconds.
    for k in range(N_PERIODS):
        offset = k * P
        ax_top.plot(t_arr + offset, L_arr, color='black', lw=1.2)
        for seg in pe_segs:
            ax_top.plot(t_arr[seg] + offset, L_arr[seg], color='red', lw=1.5)
        for seg in se_segs:
            ax_top.plot(t_arr[seg] + offset, L_arr[seg], color='blue', lw=1.5)
    legend_lines = [Line2D([0], [0], color='black', label='Full Flux'),
                     Line2D([0], [0], color='red', label='Primary Eclipse'),
                     Line2D([0], [0], color='blue', label='Secondary Eclipse')]
    ax_top.legend(handles=legend_lines, loc='lower right')

    y_candidates = [L_total]
    if pe: y_candidates.append(pe['L_min'])
    if se: y_candidates.append(se['L_min'])
    y_min = min(y_candidates)
    ax_top.set_ylim([y_min - 0.05 * y_min, 1.05 * L_total])

    # --- Choose 4 representative panel times in chronological order ---
    events = []
    if pe: events.append(('Primary Eclipse', pe['t_mid'], pe))
    if se: events.append(('Secondary Eclipse', se['t_mid'], se))
    events.sort(key=lambda ev: ev[1])

    panel_axes = [ax_p1, ax_p2, ax_p3, ax_p4]
    panel_times = []
    panel_titles = []
    if len(events) == 2:
        (name_a, t_a, info_a), (name_b, t_b, info_b) = events
        gap1_mid = 0.5 * (info_a['t_end'] + info_b['t_start'])
        gap2_start = info_b['t_end']
        gap2_end = info_a['t_start'] + P
        gap2_mid = 0.5 * (gap2_start + gap2_end)
        if gap2_mid > P:
            gap2_mid -= P
        panel_times = [t_a, gap1_mid, t_b, gap2_mid]
        panel_titles = [name_a, "Full Flux", name_b, "Full Flux"]
    else:
        # only one type of eclipse detected (grazing/geometry-limited case)
        name_a, t_a, info_a = events[0]
        gap_mid = 0.5 * (info_a['t_end'] + (info_a['t_start'] + P))
        if gap_mid > P:
            gap_mid -= P
        panel_times = [t_a, gap_mid, gap_mid, gap_mid]
        panel_titles = [name_a, "Full Flux", "Full Flux", "Full Flux"]

    ax_top.set_title(
        f"{target}\n"
        rf"m₁ = {m1} $M_\odot$, r₁ = {r1} $R_\odot$, L₁ = {L1} $L_\odot$,    "
        rf"m₂ = {m2} $M_\odot$, r₂ = {r2} $R_\odot$, L₂ = {L2} $L_\odot$" + "\n"
        f"P = {P/(24*60**2):.3f} d, a = {sma/AU:.4f} AU, e = {e:.3f}, ω = {omega:.1f}°, i = {i}°\n"
        + (f"Eclipse Duration {pe['duration']/60:.2f} min, b = {pe['d_min']/r1:.3f}   " if pe or se else "No Eclipse Occurs  ")
    )

ax_top.set_xlabel("Seconds")
ax_top.set_ylabel("Solar Luminosities")
ax_top.set_xlim([0, N_PERIODS * P])

# --------------------------------------------------
# Orbit / eclipse-geometry panels
lim = 1.5 * a_Rsol * (1 + e)  # generous view window that fits the full ellipse
nu_full = np.linspace(0, 2 * np.pi, 500)
r_full = a_Rsol * (1 - e**2) / (1 + e * np.cos(nu_full))
om = np.radians(omega)
inc = np.radians(i)
orbit_x = r_full * np.cos(om + nu_full)
orbit_y = r_full * np.sin(om + nu_full) * np.cos(inc)

dt = t_arr[1] - t_arr[0]
for ax, t_val, title in zip(panel_axes, panel_times, panel_titles):
    idx_t = int(round((t_val % P) / dt)) % N_SAMPLES
    x_t, y_t, z_t, d_t = float(x_arr[idx_t]), float(y_arr[idx_t]), float(z_arr[idx_t]), float(d_arr[idx_t])
    L_t = float(L_arr[idx_t])

    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])
    ax.set_xlabel(r"Solar Radii $R_\odot$")
    ax.set_ylabel(r"Solar Radii $R_\odot$")
    ax.set_title(rf"{title}" + "\n" + rf"t = {t_val/60:.2f} min" + "\n" + rf"L = {L_t:.2f} $L_\odot$")

    ax.plot(orbit_x, orbit_y, color='black', lw=1, zorder=3)
    star1 = Circle((0, 0), r1, color=primary_color, label='m1', zorder=(2 if z_t >= 0 else 4))
    star2 = Circle((x_t, y_t), r2, color=secondary_color, label='m2', zorder=(4 if z_t >= 0 else 2))
    ax.add_patch(star1)
    ax.add_patch(star2)

panel_axes[-1].legend(
    handles=[Line2D([0], [0], color='black', lw=1, label='Relative Orbit'),
             Line2D([0], [0], marker='o', color='w', markerfacecolor=primary_color, markersize=10, label='m1'),
             Line2D([0], [0], marker='o', color='w', markerfacecolor=secondary_color, markersize=10, label='m2')],
    loc="lower center", bbox_to_anchor=(-1.5, -0.35), ncol=3
)

print(f"Orbital Period: {P/(24*60*60):.3f} days")
print(f"Semi-major axis: {sma/AU:.4f} AU   (periastron: {r_peri:.3f} R☉, apastron: {a_Rsol*(1+e):.3f} R☉)")
print(f"Eccentricity: {e:.3f}   Argument of Periastron: {omega:.1f}°")
if eclipses_occur:
    if pe:
        print(f"Primary Eclipse:   duration {pe['duration']/60:.3f} min, min separation {pe['d_min']:.3f} R☉, L_min = {pe['L_min']:.3f} L☉")
    if se:
        print(f"Secondary Eclipse: duration {se['duration']/60:.3f} min, min separation {se['d_min']:.3f} R☉, L_min = {se['L_min']:.3f} L☉")
else:
    print("No eclipses occur for this geometry.")

# --------------------------------------------------
# Analytic eclipse-geometry metrics (generalized from LC_v4.5 to handle e > 0).
#
# LC_v4.5 assumed a circular orbit, so a single semi-major-axis-based impact parameter
# and transit duration applied to both eclipses. For an eccentric orbit, the two
# conjunctions generally sit at different separations (unless omega = 0 deg or 180 deg),
# so impact parameter, transit duration, and the eclipse/grazing inclination thresholds
# are computed separately for the primary and secondary conjunctions, using the
# already-identified (nu_c, r_c) for each from the depth-based pe/se assignment above.
#
# Note: while porting this over we also found and fixed a bug in the original transit
# duration formula -- it omitted a factor of 1/sin(i), which for the circular case alone
# undercounted duration by ~1.2% (597.5 min vs. the true 604.9 min, confirmed against
# this script's own numerically-sampled eclipse duration above).
inc_rad = np.radians(i)
Rsum = r1 + r2
Rdiff = r1 - r2   # signed, not abs() -- matches LC_v4.5's convention exactly

def eclipse_geometry(r_c, nu_c):
    """Analytic impact parameter, transit duration, and inclination thresholds for a
    single conjunction at true anomaly nu_c with instantaneous separation r_c (R_Sol)."""
    if r_c is None:
        return None
    b_c = r_c * np.cos(inc_rad)                                   # impact parameter, R_Sol (signed)
    sin_i = np.sin(inc_rad)
    denom_c = 1 + e * np.cos(nu_c)                                # (1+e sin omega) or (1-e sin omega)

    if sin_i > 0:
        arg = np.clip(np.sqrt(max(Rsum**2 - b_c**2, 0.0)) / (r_c * sin_i), -1.0, 1.0)
        half_angle = np.arcsin(arg)
    else:
        half_angle = 0.0
    # Local-velocity approximation (exact for e=0; treats r as ~constant across the
    # transit window, via Kepler's 2nd law angular-velocity scaling at this conjunction)
    duration = half_angle * P * (1 - e**2)**1.5 / (np.pi * denom_c**2)   # seconds

    i_min = np.degrees(np.arccos(np.clip(Rsum / r_c, -1.0, 1.0))) if r_c > 0 else np.nan
    i_grazing = np.degrees(np.arccos(np.clip(Rdiff / r_c, -1.0, 1.0))) if r_c > 0 else np.nan
    return {'b': b_c, 'duration': duration, 'i_min': i_min, 'i_grazing': i_grazing}

geo_pe = eclipse_geometry(r_c_pe, nu_c_pe)
geo_se = eclipse_geometry(r_c_se, nu_c_se)

# Minimum possible orbital period: smallest P (at this e) keeping periastron >= r1+r2,
# i.e. a_min = (r1+r2)/(1-e) -- generalizes LC_v4.5's a_min = r1+r2 (its e=0 special case)
a_min_Rsol = Rsum / (1 - e)
P_min = np.sqrt((a_min_Rsol * R_Sol / AU)**3 / (m1 + m2)) * yr   # seconds

if geo_pe:
    print(f"Primary Transit Duration: {geo_pe['duration']/60:.3f} minutes")
    print(f"Primary Impact Parameter: {geo_pe['b']:.3f} R☉,    b/r₁ = {geo_pe['b']/r1:.3f}")
    print(f"Primary Minimum Inclination for Eclipse: {geo_pe['i_min']:.2f}° < i < {180-geo_pe['i_min']:.2f}°")
    print(f"Primary Minimum Grazing Eclipse Inclination: {geo_pe['i_grazing']:.2f}° < i < {180-geo_pe['i_grazing']:.2f}°")
if geo_se:
    print(f"Secondary Transit Duration: {geo_se['duration']/60:.3f} minutes")
    print(f"Secondary Impact Parameter: {geo_se['b']:.3f} R☉,    b/r₁ = {geo_se['b']/r1:.3f}")
    print(f"Secondary Minimum Inclination for Eclipse: {geo_se['i_min']:.2f}° < i < {180-geo_se['i_min']:.2f}°")
    print(f"Secondary Minimum Grazing Eclipse Inclination: {geo_se['i_grazing']:.2f}° < i < {180-geo_se['i_grazing']:.2f}°")
print(f"Minimum Possible Orbital Period: {P_min/(24*60*60):.3f} <= P < {P/(24*60*60):.3f} days")

fig.savefig(f'{target}_{(P/86400):.3f}d_{sma/AU:.3f}AU_{e:.3f}.png', dpi=500, bbox_inches='tight')
print(f"Graph saved as '{target}_{(P/86400):.3f}d_{sma/AU:.3f}AU_{e:.3f}.png'")