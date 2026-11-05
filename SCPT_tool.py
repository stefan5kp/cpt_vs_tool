# SCPT_tool: shared functions for the All-in-one tool, Ergodic model and Layer-Based Regression examples.
# Import packages for use:
import numpy as np
import ngl_tools.smt as smt
import pwlf
from itertools import product


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

def regression(tt_depth, Meas_TT, n, source = 1.0, search_radius = 2):

    # Correction of Source Offset Distance
    traveltime = tt_depth / np.sqrt (tt_depth**2 + source**2) * Meas_TT

    # Compute half depths from tt_depth
    half_depths = (tt_depth[:-1] + tt_depth[1:]) / 2

    # joint the half depths with tt_depth (without the first and last point) as candidate breakpoints
    candidate_breakpoints = np.sort(np.concatenate((tt_depth[1:-1], half_depths)))

    n_segments = n

    # One segment
    if n_segments == 1:

        piecewise = pwlf.PiecewiseLinFit(tt_depth, traveltime)
        breaks = piecewise.fit(1)
        fitted_values = piecewise.predict(tt_depth)
        fitted_values_for_plot = fitted_values
        new_depths = tt_depth

    else:
        pw = pwlf.PiecewiseLinFit(tt_depth, traveltime)
        optimized_breaks = pw.fit(n_segments)
        internal_breaks = (optimized_breaks[1:-1])

        candidate_lists = []

        for b in internal_breaks:

            # Find the index of the closest candidate breakpoint to the optimized break
            idx = np.argmin(np.abs(candidate_breakpoints - b))
            lb = max(0, idx - search_radius)
            ub = min(len(candidate_breakpoints), idx + search_radius + 1)

            nearby = (candidate_breakpoints[lb : ub])
            candidate_lists.append(nearby)

        # Only consider combinations from the candidate lists
        best_rss = np.inf
        best_breaks = None
        best_piecewise = None

        for combo in product(*candidate_lists):

            internal = np.unique(combo)

            if len(internal) != (n_segments - 1):
                continue

            all_breaks = np.array(sorted([tt_depth.min()] + list(internal) + [tt_depth.max()]))

            # minimum 2 pts/segment
            valid = True

            for i in range(len(all_breaks) - 1):
                count = np.sum((tt_depth >= all_breaks[i]) & (tt_depth <= all_breaks[i + 1]))
                if count < 2:
                    valid = False
                    break
            if not valid:
                continue

            try:
                piecewise_trial = (pwlf.PiecewiseLinFit(tt_depth, traveltime))
                rss = (piecewise_trial.fit_with_breaks(all_breaks))
                if rss < best_rss:
                    best_rss = rss
                    best_breaks = (all_breaks)
                    best_piecewise = (piecewise_trial)
            except:
                continue

        if best_breaks is None:
            raise ValueError("No valid breakpoint combination found")

        breaks = best_breaks
        piecewise = best_piecewise

        # Predictions
        fitted_values = piecewise.predict(tt_depth)
        new_depths = np.array(sorted(set(breaks.tolist() + tt_depth.tolist())))
        fitted_values_for_plot = (piecewise.predict(new_depths))

    # Calculatae RMSE
    residuals = (traveltime - fitted_values)
    rmse = np.sqrt(np.mean(residuals**2))

    # Calculate slope of each segment
    slopes = piecewise.slopes
    ztop = breaks[:-1].copy()
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


def optimize (tt_depth, Meas_TT):

    if len(tt_depth) < 7:
        error = np.zeros(len(tt_depth)-1)
        for i in range (1, len(tt_depth)):
            results = regression(tt_depth, Meas_TT, i)
            error[i-1] = results['rmse']

    else:
        error = np.zeros(6)
        for i in range (1, 7):
            results = regression(tt_depth, Meas_TT, i)
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

    return Optimum


def optimized_regression(tt_depth, Meas_TT, source = 1.0):
    n = optimize(tt_depth, Meas_TT)
    results = regression(tt_depth, Meas_TT, n, source = source)

    velocity = (1/np.array(results['slopes']))*1000

    return velocity, results


def Master (tt_depth, Meas_TT, slope_break_method = 1, breakpoints = None, source = 1.0):

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
