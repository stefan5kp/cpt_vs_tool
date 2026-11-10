from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import cpt_vs_tool as cvt  # noqa: E402

DATA = Path(__file__).parent / "data"


def load_avn_1():
    pd = pytest.importorskip("pandas")
    cpt = pd.read_csv(DATA / "avn_1.csv")
    tt = pd.read_csv(DATA / "avn_1_tt.csv")
    cpt = cpt.dropna(subset=["depth (m)", "qt (MPa)", "fs (MPa)"])
    return (
        cpt["depth (m)"].to_numpy(dtype=float),
        cpt["qt (MPa)"].to_numpy(dtype=float),
        cpt["fs (MPa)"].to_numpy(dtype=float),
        tt["depth (m)"].to_numpy(dtype=float),
        tt["traveltime (ms)"].to_numpy(dtype=float),
    )


def test_version():
    assert isinstance(cvt.__version__, str)


def test_interleave():
    out = cvt.interleave(np.array([1.0, 2.0]), np.array([3.0, 4.0]))
    np.testing.assert_array_equal(out, [1.0, 3.0, 2.0, 4.0])


def test_lbr_recovers_two_layer_velocities():
    # Two layers: 200 m/s down to 10 m, 400 m/s below; zero source offset
    depth = np.arange(1.0, 21.0, 1.0)
    tt_ms = np.where(depth <= 10, depth / 200, 10 / 200 + (depth - 10) / 400) * 1000
    velocity, _, results, ztop, zbot = cvt.LBR(depth, tt_ms, slope_break_method=1, breakpoints=[10.0], source=0.0)
    np.testing.assert_allclose(velocity, [200.0, 400.0], rtol=1e-6)
    assert list(zbot) == [10.0, 20.0]


def test_vsz_layered():
    vsz, z = cvt.VsZ_layered([10.0, 20.0], [200.0, 400.0])
    assert z == 20.0
    assert vsz == pytest.approx(20.0 / (10 / 200 + 10 / 400))


def test_compute_profiles_avn_1():
    cpt_depth, qt, fs, tt_depth, meas_tt = load_avn_1()
    p = cvt.compute_profiles(cpt_depth, qt, fs, 31.0, tt_depth, meas_tt, "Qal1")
    assert p["has_tt"]
    assert p["a_geo"] == cvt.a_coeff["Qal1"]
    # Same value as the saved output of interpretation_tool/scpt_tool_v1.ipynb
    assert p["a_site_specific"] == pytest.approx(4.1430, abs=1e-3)
    assert np.all(np.isfinite(p["global_velocity"]))
    assert np.all(p["LBR_velocity"] > 0)
    vsz = cvt.compute_VsZ(p)
    assert set(vsz) == {"Zhang et al.", "Zhang et al. (ss)", "Robertson (2012)", "LBR"}


def test_compute_profiles_without_travel_time():
    cpt_depth, qt, fs, _, _ = load_avn_1()
    p = cvt.compute_profiles(cpt_depth, qt, fs, 31.0, None, None, "Qal1")
    assert not p["has_tt"]
    assert np.isnan(p["a_site_specific"])
    fig = cvt.plot_profiles(p)
    assert len(fig.axes) == 4
