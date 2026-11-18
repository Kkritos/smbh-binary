"""Polytrope sampling used to seed the debris (debris_tracer.py and mutde_lab.py)."""
import numpy as np
import debris_tracer as dt
from conftest import load_module


def test_lane_emden_first_zero_known_values():
    # n = 1 has the closed form theta = sin(xi)/xi with first zero at pi; n = 3: xi_1 = 6.89685
    assert np.isclose(dt.lane_emden(1.0)[2], np.pi, rtol=1e-4)
    assert np.isclose(dt.lane_emden(3.0)[2], 6.89685, rtol=1e-4)


def test_polytrope_positions_inside_star_and_centrally_concentrated():
    rng = np.random.default_rng(1)
    dr = dt.sample_polytrope_positions(20000, 3.0, 1.0, rng)
    r = np.linalg.norm(dr, axis=1)
    assert dr.shape == (20000, 3) and r.max() <= 1.0
    # an n = 3 polytrope holds ~half its mass within ~0.24 R_star (uniform sphere: 0.79)
    assert 0.18 < np.median(r) < 0.32
    assert np.allclose(dr.mean(axis=0), 0.0, atol=0.01)        # isotropic


def test_mutde_lab_sampler_agrees_with_tracer():
    lab = load_module("mutde_lab", "analysis/debris_mutde/mutde_lab.py")
    r_lab = np.linalg.norm(lab.sample_polytrope_n3(20000, 1.0, 3.0, np.random.default_rng(2)), axis=1)
    r_dt = np.linalg.norm(dt.sample_polytrope_positions(20000, 3.0, 1.0, np.random.default_rng(3)), axis=1)
    assert abs(np.median(r_lab) - np.median(r_dt)) < 0.02
