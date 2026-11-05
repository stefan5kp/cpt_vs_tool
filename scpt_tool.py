# scpt_tool: shared functions for the interpretation tool and Layer-Based Regression examples.
# Import packages for use:
import numpy as np
import matplotlib.pyplot as plt
import ngl_tools.smt as smt
import pwlf
from itertools import combinations, chain, islice
from scipy.optimize import least_squares


# =============================================================================
# Utilities
# =============================================================================

def interleave(v1, v2):
    vout = np.empty(2 * len(v1), dtype=float)
    vout[0::2] = v1
    vout[1::2] = v2
    return vout


# =============================================================================
# CPT processing
# =============================================================================

def gamma_t (depth, fs, qt):
    # Replace invalid values in CPT_fs and CPT_qt with NaN
    fs1 = fs.copy()
    qt1 = qt.copy()
    fs1[fs1 <= 0] = np.nan
    qt1[qt1 <= 0] = np.nan
    # Compute gamma only for valid values
    valid_mask = ~np.isnan(fs1) & ~np.isnan(qt1)
    gamma = np.full_like(depth, np.nan)  # Initialize with NaN
    # Calculate gamma using Robertson (2010)
    gamma[valid_mask] = 9.81 * (0.27 * np.log10(fs1[valid_mask]/qt1[valid_mask] * 100) + 0.36 * np.log10(qt1[valid_mask] / 101.3) + 1.236)
     # Interpolate NaN values
    gamma = np.interp(depth, depth[valid_mask], gamma[valid_mask])

    return gamma


def get_CPT (CPT_depth, CPT_qt, CPT_fs, dGWT):

    CPT_fs = CPT_fs * 1000 #Convert to kPa
    CPT_qt = CPT_qt * 1000 #Convert to kPa
    pa = 101.3 #KPa,
    gamma = gamma_t(CPT_depth, CPT_fs, CPT_qt)

    CPT_fs[CPT_fs < 0] = 0.01
    CPT_qt[CPT_qt < 0] = 0.01

    # integrate gamma over CPT_depth to get sigmav, then compute sigmavp by subtracting the pore pressure u
    # updated 2026-02-24
    depth_interval = np.diff(CPT_depth, prepend=0)
    sigmav = np.cumsum(gamma * depth_interval)

    u = (CPT_depth - dGWT) * 9.81
    u[u<0] = 0
    sigmavp = sigmav - u

    # Compute fz
    fz = sigmavp/pa

    # Calculate Ic, Qtn and their inverse filtered values
    Ic, Qtn, Fr = smt.get_Ic_Qtn_Fr(CPT_qt, CPT_fs, sigmav, sigmavp)
    qt_inv, fs_inv, Ic_inv = smt.cpt_inverse_filter(CPT_qt, CPT_depth, fs=CPT_fs, sigmav=sigmav, sigmavp=sigmavp)
    Ic_inv, Qtn_inv, Fr_inv = smt.get_Ic_Qtn_Fr(qt_inv, fs_inv, sigmav, sigmavp)

    # Ic correction for Qtn (Zhang et al. kc)
    kc = 1 + 5 / (1 + np.exp(-ZHANG_KC[0] * Ic_inv + ZHANG_KC[1]))
    Qtncs = Qtn_inv * kc

    return fz, Ic_inv, Qtncs, CPT_depth, Qtn_inv, Ic_inv, Qtn, Ic, kc


# =============================================================================
# CPT-Vs models
# =============================================================================

# Zhang et al. model coefficients: ln(Vs) = b0 * ln(kc * Qtn) + m(Ic) * ln(sigma_v'/pa) + a
ZHANG_B = [0.224117, 1.278293, 3.385425]   # b0, and m(Ic) = 0.1 + 0.4 / (1 + exp(-b1 * Ic + b2))
ZHANG_KC = [5.482742, 12.320768]           # kc(Ic) = 1 + 5 / (1 + exp(-c1 * Ic + c2))

# Geology-specific model offset "a" for the Zhang et al. model
ZHANG_A = {'Beach and Dune': 4.402626,
           'Floodplain':     4.255830,
           'Glacial':        4.195240,
           'Lagoonal':       4.195468,
           'Valley Train':   4.238807,
           'af/Qi':          4.316504,
           'Qal1':           4.350243,
           'Qal2':           4.411610,
           'Qal3':           4.383580,
           'Qoa':            4.515429,
           'Qs':             4.420910}


def robertson_2012(Qtn, fz, Ic):
    # Robertson (2012): ln(Vs) = 1.93 + 0.5 * ln(Qtn) + 0.25 * ln(sigma_v'/pa) + 0.63 * Ic
    return np.exp(1.93 + 0.5 * np.log(Qtn) + 0.25 * np.log(fz) + 0.63 * Ic)


def zhang_lnVs(Qcncs, fz, Ic, a, b = ZHANG_B):
    # Zhang et al. model, returns ln(Vs). Qcncs is the kc-corrected Qtn from get_CPT
    m = 0.1 + 0.4 / (1 + np.exp(-b[1] * Ic + b[2]))
    return b[0] * np.log(Qcncs) + m * np.log(fz) + a


# =============================================================================
# Travel time processing and CPT-to-layer averaging
# =============================================================================

def get_TT (TT_depth, Meas_TT, source = 1.0):
    # Correction of Source Offset Distance
    tt_corrected = TT_depth / np.sqrt (TT_depth**2 + source**2) * Meas_TT
    ztop = TT_depth[:-1]
    zbot = TT_depth[1:]
    DTT = np.diff(tt_corrected)
    return ztop, zbot, DTT


def pair_layers(CPT_depth, ztop, zbot):
    # Pair each CPT point with every layer that contains it. Layers are closed intervals [ztop, zbot], so a point
    # on the boundary between two layers is used in both.
    in_layer = (CPT_depth[:, None] >= ztop[None, :]) & (CPT_depth[:, None] <= zbot[None, :])
    pair_cpt, pair_layer = np.nonzero(in_layer)
    N_layer = np.bincount(pair_layer, minlength=len(ztop))
    return pair_cpt, pair_layer, N_layer


def layer_travel_time(lnVs, CPT_depth, ztop, zbot):
    # Predicted travel time (s) through each layer, each CPT point represents an equal share of its layer
    pair_cpt, pair_layer, N_layer = pair_layers(CPT_depth, ztop, zbot)
    cpt_slowness = 1.0 / np.exp(lnVs[pair_cpt])
    sum_slowness = np.bincount(pair_layer, weights=cpt_slowness, minlength=len(ztop))
    tt_hat = np.full(len(ztop), np.nan)
    has_pts = N_layer > 0
    tt_hat[has_pts] = (zbot - ztop)[has_pts] / N_layer[has_pts] * sum_slowness[has_pts]
    return tt_hat


def layer_average_velocity(lnVs, CPT_depth, ztop, zbot):
    # Travel-time (slowness) average of CPT-based ln(Vs) over each layer (NaN for layers with no CPT points)
    return (zbot - ztop) / layer_travel_time(lnVs, CPT_depth, ztop, zbot)


def site_specific(x, CPT_depth, Qtncs, fz, Ic, ztop, zbot, DTT):
    # Residual function for fitting the site-specific "a" of the Zhang et al. model to the measured travel times.
    # Layers with no CPT points (e.g. below the end of the sounding) or a non-positive travel time are dropped.
    tt_hat = layer_travel_time(zhang_lnVs(Qtncs, fz, Ic, x), CPT_depth, ztop, zbot)
    valid = ~np.isnan(tt_hat) & (DTT > 0)
    return np.log(DTT[valid] / 1000) - np.log(tt_hat[valid])


# =============================================================================
# Layer-based regression of travel time data
# =============================================================================

def _best_breaks(tt_depth, traveltime, candidates, n_bp, chunk = 200_000, n_keep = 20):
    # Exhaustive search over all n_bp-combinations of candidate breakpoints (same as looping
    # pwlf.fit_with_breaks), but vectorized: the continuous piecewise-linear fit uses the hinge
    # basis [1, x - xmin, (x - c1)+, ...], so the normal equations for any combo are a submatrix
    # of one precomputed Gram matrix. Combos are scored in chunks, and the best few are re-checked
    # exactly with lstsq to avoid any round-off in the normal equations.
    x = np.asarray(tt_depth, dtype=float)
    y = np.asarray(traveltime, dtype=float)
    candidates = np.asarray(candidates, dtype=float)

    X = np.column_stack([np.ones_like(x), x - x.min(), np.maximum(x[:, None] - candidates[None, :], 0)])
    G = X.T @ X
    Xty = X.T @ y
    yty = y @ y

    # Number of points in [b_i, b_i+1] (inclusive) = hi[b_i+1] - lo[b_i]
    xs = np.sort(x)
    lo = np.searchsorted(xs, candidates, side='left')
    hi = np.searchsorted(xs, candidates, side='right')
    lo_min = np.searchsorted(xs, x.min(), side='left')
    hi_max = np.searchsorted(xs, x.max(), side='right')

    combo_iter = combinations(range(len(candidates)), n_bp)
    kept_combos, kept_rss = [], []
    while True:
        flat = np.fromiter(chain.from_iterable(islice(combo_iter, chunk)), dtype=np.int64)
        if flat.size == 0:
            break
        combos = flat.reshape(-1, n_bp)

        # Ensure each segment has at least two points
        seg_lo = np.column_stack([np.full(len(combos), lo_min), lo[combos]])
        seg_hi = np.column_stack([hi[combos], np.full(len(combos), hi_max)])
        combos = combos[np.all(seg_hi - seg_lo >= 2, axis=1)]
        if len(combos) == 0:
            continue

        cols = np.column_stack([np.zeros(len(combos), dtype=np.int64), np.ones(len(combos), dtype=np.int64), combos + 2])
        Gs = G[cols[:, :, None], cols[:, None, :]]
        b = Xty[cols]
        try:
            beta = np.linalg.solve(Gs, b[..., None])[..., 0]
        except np.linalg.LinAlgError:
            beta = (np.linalg.pinv(Gs) @ b[..., None])[..., 0]
        rss = yty - np.sum(b * beta, axis=1)

        k = min(n_keep, len(rss))
        idx = np.sort(np.argpartition(rss, k - 1)[:k])
        kept_combos.append(combos[idx])
        kept_rss.append(rss[idx])

    if not kept_combos:
        return None

    # Exact re-check of the shortlisted combos
    shortlist = np.concatenate(kept_combos)
    approx = np.concatenate(kept_rss)
    shortlist = shortlist[np.sort(np.argsort(approx, kind='stable')[:n_keep])]
    best_rss, best_combo = np.inf, None
    for combo in shortlist:
        Xs = X[:, np.concatenate(([0, 1], combo + 2))]
        beta = np.linalg.lstsq(Xs, y, rcond=None)[0]
        rss = np.sum((y - Xs @ beta) ** 2)
        if rss < best_rss:
            best_rss, best_combo = rss, combo

    return [x.min()] + list(candidates[best_combo]) + [x.max()]


def regression (tt_depth, Meas_TT, n, source = 1.0):

    # Correction of Source Offset Distance
    traveltime = tt_depth / np.sqrt (tt_depth**2 + source**2) * Meas_TT

    # Compute half depths from tt_depth
    half_depths = (tt_depth[:-1] + tt_depth[1:]) / 2

    # joint the half depths with tt_depth (without the first and last point) as candidate breakpoints
    candidate_breakpoints = np.sort(np.concatenate((tt_depth[1:-1], half_depths)))

    n_segments = n

    if n_segments == 1:
        piecewise = pwlf.PiecewiseLinFit(tt_depth, traveltime)
        breaks = piecewise.fit(n_segments)
        fitted_values = piecewise.predict(tt_depth)
        fitted_values_for_plot = fitted_values
        new_depths = tt_depth
    elif n_segments - 1 > len(half_depths):
        raise ValueError("Not enough candidate breakpoints for the number of segments")
    else: 
        # Search all sorted combinations of candidate breakpoints for the lowest RSS
        breaks = _best_breaks(tt_depth, traveltime, sorted(set(candidate_breakpoints)), n_segments - 1)

        if breaks is None:
            raise ValueError(f"No valid breakpoint combination found for the data with n={n_segments}")

        # Refit the best combination with pwlf for the outputs
        piecewise = pwlf.PiecewiseLinFit(tt_depth, traveltime)
        piecewise.fit_with_breaks(breaks)

        # Combine breaks and tt_depth, then sort them
        new_depths = np.array(sorted(set(breaks + list(tt_depth))))

        # Generate the fitted values
        fitted_values = piecewise.predict(tt_depth)
        fitted_values_for_plot = piecewise.predict(new_depths)
    
    # Calculate RMSE
    residuals = traveltime - fitted_values
    rmse = np.sqrt(np.mean(residuals**2))

    # Calculate slopes of each segment
    slopes = piecewise.slopes
    ztop = np.array(breaks[:-1], dtype=float)
    # Replace the first ztop with 0.01
    ztop[0] = 0.01
    zbot = breaks[1:]
    
    return {
        'breaks': breaks,
        'candidate_breakpoints': candidate_breakpoints,
        'fitted_values': fitted_values,
        'fitted_values_for_plot': fitted_values_for_plot,
        'new_depths': new_depths,
        'rmse': rmse,
        'slopes': slopes,
        'ztop': ztop,
        'zbot': zbot
    }

def regression_full_manual(tt_depth, Meas_TT, breakpoints, source = 1.0):

    # Correction of Source Offset Distance
    traveltime = tt_depth / np.sqrt (tt_depth**2 + source**2) * Meas_TT

    # Validate breakpoints
    breakpoints = sorted(set(breakpoints))  # Ensure sorted and unique
    if breakpoints[0] <= tt_depth.min() or breakpoints[-1] >= tt_depth.max():
        raise ValueError("Breakpoints must be within the range of tt_depth (excluding endpoints)")

    all_breaks = [tt_depth.min()] + breakpoints + [tt_depth.max()]

    # Ensure each segment has at least two data points
    for i in range(len(all_breaks) - 1):
        count = np.sum((tt_depth >= all_breaks[i]) & (tt_depth <= all_breaks[i + 1]))
        if count < 2:
            raise ValueError(f"Segment between {all_breaks[i]} and {all_breaks[i + 1]} has less than 2 data points.")

    # Fit with custom breakpoints
    piecewise = pwlf.PiecewiseLinFit(tt_depth, traveltime)
    rss = piecewise.fit_with_breaks(all_breaks)

    # Generate fitted values
    new_depths = np.array(sorted(set(all_breaks + list(tt_depth))))
    fitted_values = piecewise.predict(tt_depth)
    fitted_values_for_plot = piecewise.predict(new_depths)

    # RMSE
    residuals = traveltime - fitted_values
    rmse = np.sqrt(np.mean(residuals**2))

    # Slopes
    slopes = piecewise.slopes
    ztop = all_breaks[:-1]
    # Replace the first ztop with 0.01
    ztop[0] = 0.01

    zbot = all_breaks[1:]

    return {
        'breaks': all_breaks,
        'fitted_values': fitted_values,
        'fitted_values_for_plot': fitted_values_for_plot,
        'new_depths': new_depths,
        'rmse': rmse,
        'slopes': slopes,
        'ztop': ztop,
        'zbot': zbot
    }


def optimize (tt_depth, Meas_TT, return_results = False):

    all_results = []
    if len(tt_depth) < 7:
        error = np.zeros(len(tt_depth)-1)
        for i in range (1, len(tt_depth)):
            results = regression(tt_depth, Meas_TT, i)
            all_results.append(results)
            error[i-1] = results['rmse']

    else:
        error = np.zeros(6)
        for i in range (1, 7):
            results = regression(tt_depth, Meas_TT, i)
            all_results.append(results)
            error[i-1] = results['rmse']

    # Compute cost function (Updated on 2024-06-07)
    t_avg = (tt_depth.max() - tt_depth.min()) / np.arange(1, len(error)+1)
    tref = 4

    wR = 1
    wT = 1

    # Normalized RMSE Cost Function
    JR = error / error[0]
    # Thickness-Dependent Cost Function
    JT = 0.2 * (tref / t_avg) ** 3
    # Combined Cost Function
    cost = wR * JR + wT * JT

    min_cost_index = np.argmin(cost)
    Optimum = min_cost_index + 1

    if return_results:
        return Optimum, all_results[min_cost_index]
    return Optimum
    

def optimized_regression(tt_depth, Meas_TT, source = 1.0):
    n, results = optimize(tt_depth, Meas_TT, return_results = True)
    # optimize() fits with the default source = 1.0; reuse that fit instead of repeating it
    if source != 1.0:
        results = regression(tt_depth, Meas_TT, n, source = source)

    velocity = (1/np.array(results['slopes']))*1000

    return velocity, results

def LBR(tt_depth, Meas_TT, slope_break_method = 1, breakpoints = None, source = 1.0):

    tt_corrected = tt_depth / np.sqrt (tt_depth**2 + source**2) * Meas_TT

    # Full Automatic
    if slope_break_method == 0: 
        velocity, results = optimized_regression(tt_depth, Meas_TT, source = source)
    # Full Manual
    elif slope_break_method == 1:
        if breakpoints is None:
            raise ValueError("Breakpoints must be provided for manual regression.")
        results = regression_full_manual(tt_depth, Meas_TT, breakpoints, source = source)
        velocity = (1/np.array(results['slopes']))*1000

    top_depth = results['ztop']
    bottom_depth = results['zbot']

    return velocity, tt_corrected, results, top_depth, bottom_depth


# Earlier name of LBR, kept for the lbr/ example notebook
Master = LBR


# =============================================================================
# Interpretation tool: all four methods for one SCPT sounding
# =============================================================================

def predicted_DTT(lnVs, CPT_depth, ztop, zbot, DTT):
    # Predicted travel time (ms) through each layer and its half depth, for layers with CPT points and DTT > 0
    tt_hat = layer_travel_time(lnVs, CPT_depth, ztop, zbot) * 1000
    valid = ~np.isnan(tt_hat) & (DTT > 0)
    return tt_hat[valid], 0.5 * (ztop[valid] + zbot[valid])


def compute_profiles(CPT_depth, CPT_qc, CPT_qt, CPT_fs, dGWT, TT_depth, Meas_TT, geology,
                     slope_break_method = 0, breakpoints = None, source = 1.0):
    # Zhang et al. (global and site-specific) with qc, Robertson (2012) with qt, and LBR of the travel time data

    # Zhang et al. inputs: inverse-filtered Qtn and Ic, kc-corrected Qtn
    fz, Ic, Qtncs, CPT_depth_all, Qtn_inv, Ic_inv, Qtn, Ic_raw, kc = get_CPT(CPT_depth, CPT_qc, CPT_fs, dGWT)
    valid = (fz > 0) & (Qtncs > 0) & (Ic >= 1)
    CPT_depth, fz, Ic, Qtn_inv, Qtncs, kc = CPT_depth_all[valid], fz[valid], Ic[valid], Qtn_inv[valid], Qtncs[valid], kc[valid]

    # Global velocity with the geology-specific a
    a_geo = ZHANG_A[geology]
    global_velocity = np.exp(zhang_lnVs(Qtncs, fz, Ic, a_geo))

    # Robertson (2012) with Qtn and Ic that are not inverse-filtered
    fz_rob, _, _, _, _, _, Qtn_rob, Ic_rob, _ = get_CPT(CPT_depth_all, CPT_qt, CPT_fs, dGWT)
    Robertson_Vs = robertson_2012(Qtn_rob, fz_rob, Ic_rob)

    # Measured travel time layers and interval velocity
    ztop_TT, zbot_TT, DTT = get_TT(TT_depth, Meas_TT, source = source)
    TT = TT_depth / np.sqrt(TT_depth**2 + source**2) * Meas_TT
    halfdepth = 0.5 * (ztop_TT + zbot_TT)
    Vsi = (zbot_TT - ztop_TT) / DTT * 1000

    # Site-specific a fitted to the measured travel times
    result = least_squares(site_specific, x0 = 4.0, args = (CPT_depth, Qtncs, fz, Ic, ztop_TT, zbot_TT, DTT))
    a_site_specific = result.x[0]
    site_specific_velocity = np.exp(zhang_lnVs(Qtncs, fz, Ic, a_site_specific))
    tt_hat_ss, halfdepth_ss = predicted_DTT(np.log(site_specific_velocity), CPT_depth, ztop_TT, zbot_TT, DTT)
    tt_hat_global, halfdepth_global = predicted_DTT(np.log(global_velocity), CPT_depth, ztop_TT, zbot_TT, DTT)

    # LBR velocity (automatic, or manual with the selected breakpoints)
    LBR_velocity, tt_corrected, results, LBR_ztop, LBR_zbot = LBR(TT_depth, Meas_TT, slope_break_method = slope_break_method,
                                                                  breakpoints = breakpoints, source = source)

    return dict(geology=geology, a_geo=a_geo, a_site_specific=a_site_specific,
                CPT_depth=CPT_depth, sigmavp=fz, Ic=Ic, Qtn=Qtn_inv, kc=kc,
                global_velocity=global_velocity, site_specific_velocity=site_specific_velocity,
                Robertson_Vs=Robertson_Vs, CPT_depth_rob=CPT_depth_all,
                TT=TT, TT_depth=TT_depth, halfdepth=halfdepth, DTT=DTT, ztop_TT=ztop_TT, zbot_TT=zbot_TT, Vsi=Vsi,
                tt_hat_ss=tt_hat_ss, halfdepth_ss=halfdepth_ss, tt_hat_global=tt_hat_global, halfdepth_global=halfdepth_global,
                LBR_velocity=LBR_velocity, results=results, LBR_ztop=LBR_ztop, LBR_zbot=LBR_zbot)


def plot_profiles(p):
    CPT_depth, Qtn, kc, sigmavp, Ic = p['CPT_depth'], p['Qtn'], p['kc'], p['sigmavp'], p['Ic']
    TT, TT_depth, halfdepth, DTT = p['TT'], p['TT_depth'], p['halfdepth'], p['DTT']
    results = p['results']

    fig, axs = plt.subplots(1, 6, figsize=(12, 5), dpi = 300)

    axs[0].plot(Qtn, CPT_depth, label='Qtn', color = 'grey', linewidth = 1.2)
    axs[0].plot(Qtn * kc, CPT_depth, label='Qtncs', color='tab:blue', linewidth=1.2)
    axs[0].set_xlim(0, 1000)
    axs[0].set_xticks(np.arange(0, 1001, 200))
    axs[0].grid(True, which="both", ls="--")
    axs[0].set_xlabel('Qtn / Qtncs')
    axs[0].set_ylabel('Depth (m)')
    axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fontsize=8)

    axs[1].plot((0, sigmavp[0]), (0, CPT_depth[0]), color='tab:blue', linestyle='dashed', linewidth = 1.2)
    axs[1].plot(sigmavp, CPT_depth, label='$\\sigma_v^\\prime$_CPT', color='tab:blue', linewidth = 1.2)
    axs[1].set_xlabel('$\\sigma_v^\\prime$ (atm)')
    axs[1].set_xticks(np.arange(0, 11, 2))
    axs[1].set_xlim(0, 10)
    axs[1].grid(True, which="both", ls="--")
    axs[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fontsize=8)

    axs[2].plot(Ic, CPT_depth, label='$I_c$', color = 'tab:blue', linewidth = 1.2)
    axs[2].set_xlim(1, 4)
    axs[2].grid(True, which="both", ls="--")
    axs[2].set_xlabel('Ic')
    axs[2].legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fontsize=8)

    axs[3].scatter(DTT, halfdepth, label='Measured $\\Delta$tt', marker='o', edgecolor='black', color = 'white', s=10)
    axs[3].scatter(p['tt_hat_global'], p['halfdepth_global'], label='Predicted $\\Delta$tt', marker='x', color='tab:blue', s=10)
    axs[3].scatter(p['tt_hat_ss'], p['halfdepth_ss'], label='Predicted $\\Delta$tt (ss)', marker='x', color='tab:orange', s=10)
    axs[3].set_xlabel('$\\Delta$tt (ms)')
    axs[3].set_xlim(0, 10)
    axs[3].set_xticks(np.arange(0, 11, 2))
    axs[3].grid(True, which="both", ls="--")
    axs[3].legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fontsize=8)

    axs[4].scatter(TT, TT_depth, label='Travel Time Data', color='black', s=10)
    colors = plt.colormaps['tab10']

    # Plot each segment with a different color
    for i in range(len(results['breaks']) - 1):
        # Define segment depth range
        z_start = results['breaks'][i]
        z_end = results['breaks'][i + 1]

        # Get depths in this segment (new_depths ensures breaks are included)
        segment_mask = (results['new_depths'] >= z_start) & (results['new_depths'] <= z_end)
        segment_depths = results['new_depths'][segment_mask]
        segment_tt = results['fitted_values_for_plot'][segment_mask]

        # Plot segment
        axs[4].plot(segment_tt, segment_depths, color=colors(i))

    axs[4].set_xlabel('Travel Time (ms)')
    axs[4].set_xlim(0, np.max(TT) * 1.1)
    if np.max(TT) > 150:
        axs[4].set_xticks(np.arange(0, np.max(TT) * 1.1, 50))
    elif np.max(TT) > 90:
        axs[4].set_xticks(np.arange(0, np.max(TT) * 1.1, 30))
    elif np.max(TT) > 50:
        axs[4].set_xticks(np.arange(0, np.max(TT) * 1.1, 20))
    else:
        axs[4].set_xticks(np.arange(0, np.max(TT) * 1.1, 10))
    axs[4].grid(True, which="both", ls="--")
    axs[4].legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fontsize=8)

    axs[5].plot(p['global_velocity'], CPT_depth, label='Zhang et al.', color='tab:blue', linewidth = 1.2)
    axs[5].plot(p['site_specific_velocity'], CPT_depth, label='Zhang et al. (ss)', color='tab:orange', linewidth = 1.2)
    axs[5].plot(p['Robertson_Vs'], p['CPT_depth_rob'], label='Robertson (2012)', color='grey', linewidth = 1.2)
    axs[5].plot(interleave(p['LBR_velocity'], p['LBR_velocity']), interleave(p['LBR_ztop'], p['LBR_zbot']), label='LBR', color='black', linewidth = 1.2)
    axs[5].set_xlim(0, 1000)
    axs[5].set_xticks(np.arange(0, 1001, 200))
    axs[5].set_xlabel('Velocity (m/s)')
    axs[5].legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), fontsize=8)
    axs[5].grid(True, which="both", ls="--")

    # Share common y-axis limits across all subplots
    y_max = np.nanmax([np.nanmax(CPT_depth), np.nanmax(p['zbot_TT']), np.nanmax(TT_depth)])
    for ax in axs:
        ax.set_ylim(y_max + 1, 0)  # inverted depth axis, same range for all

    # Panel labels (a)-(f) in the top right corner of each subplot
    for ax, panel in zip(axs, 'abcdef'):
        ax.text(0.95, 0.97, f'({panel})', transform=ax.transAxes, ha='right', va='top', fontsize=10,
                bbox=dict(facecolor='white', edgecolor='none', pad=1.5))

    plt.subplots_adjust(bottom=0.28, wspace=0.35)
    plt.tight_layout()
    plt.show()

    return fig


# VsZ (time-averaged velocity from the surface to the bottom of each profile): VsZ = Z / sum(h_i / Vs_i)
def VsZ_CPT(depth, velocity):
    # Each CPT point represents the interval above it (first point extends to the surface)
    valid = np.isfinite(depth) & np.isfinite(velocity) & (velocity > 0)
    depth = depth[valid]
    velocity = velocity[valid]
    dz = np.diff(depth, prepend=0)
    Z = depth[-1]
    return Z / np.sum(dz / velocity), Z


def VsZ_layered(zbot, velocity):
    # Layers are contiguous from the surface (first ztop = 0.01 is treated as 0)
    zbot = np.asarray(zbot, dtype=float)
    dz = np.diff(zbot, prepend=0)
    Z = zbot[-1]
    return Z / np.sum(dz / np.asarray(velocity, dtype=float)), Z


def compute_VsZ(p):
    VsZ_results = {
        'Zhang et al.': VsZ_CPT(p['CPT_depth'], p['global_velocity']),
        'Zhang et al. (ss)': VsZ_CPT(p['CPT_depth'], p['site_specific_velocity']),
        'Robertson (2012)': VsZ_CPT(p['CPT_depth_rob'], p['Robertson_Vs']),
        'LBR': VsZ_layered(p['LBR_zbot'], p['LBR_velocity']),
    }

    # VsZp: LBR time-averaged velocity from the first travel time depth to the bottom of the profile
    top_depth_lim = np.asarray(p['LBR_ztop'], dtype=float)
    top_depth_lim[0] = np.min(p['TT_depth'])
    VsZp = (np.max(p['TT_depth']) - np.min(p['TT_depth'])) / np.sum((np.asarray(p['LBR_zbot']) - top_depth_lim) / np.asarray(p['LBR_velocity']))

    return VsZ_results, VsZp
