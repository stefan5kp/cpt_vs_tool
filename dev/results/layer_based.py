# Import packages for use:
import pwlf
import numpy as np
import pandas as pd
from itertools import combinations, chain, islice


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
    ztop = breaks[:-1]
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

def interleave(v1, v2):
    vout = np.empty(2 * len(v1), dtype=float)
    vout[0::2] = v1
    vout[1::2] = v2
    return vout 