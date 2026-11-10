"""Shear wave velocity (Vs) interpretation from seismic cone penetration test (SCPT) data.

Provides the Zhang et al. (2026) CPT-Vs model (geology-specific and site-specific),
Robertson (2012), and layer-based regression (LBR) of travel time data.
"""

from importlib.metadata import PackageNotFoundError, version

from .vs import (
    LBR,
    VsZ_CPT,
    VsZ_layered,
    a_coeff,
    b_coeff,
    compute_profiles,
    compute_VsZ,
    gamma_t,
    get_CPT,
    get_TT,
    interleave,
    layer_average_velocity,
    layer_travel_time,
    optimize,
    optimized_regression,
    pair_layers,
    plot_profiles,
    predicted_DTT,
    regression,
    regression_full_manual,
    robertson_2012,
    site_specific,
    zhang_vs,
)

try:
    __version__ = version("cpt-vs-tool")
except PackageNotFoundError:  # running from a source checkout that is not installed
    __version__ = "0+unknown"

__all__ = [
    "LBR",
    "VsZ_CPT",
    "VsZ_layered",
    "a_coeff",
    "b_coeff",
    "compute_profiles",
    "compute_VsZ",
    "gamma_t",
    "get_CPT",
    "get_TT",
    "interleave",
    "layer_average_velocity",
    "layer_travel_time",
    "optimize",
    "optimized_regression",
    "pair_layers",
    "plot_profiles",
    "predicted_DTT",
    "regression",
    "regression_full_manual",
    "robertson_2012",
    "site_specific",
    "zhang_vs",
]
