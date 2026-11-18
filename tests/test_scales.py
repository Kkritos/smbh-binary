"""Campaign geometry conventions (analysis/bhstar_mutde/pkg/scales.py)."""
import numpy as np
import scales as sc


def test_rp_max_is_ten_tidal_radii_when_tidally_limited():
    r_t, r_foc, rp_max, b2 = sc.cell_scales(1e6, 1.0, 1.0)
    assert np.isclose(rp_max, 10.0 * r_t)
    assert np.isclose(b2, rp_max**2 + r_foc * rp_max)


def test_gravitational_radius_floor_at_high_mass():
    M = 1e10
    _, _, rp_max, _ = sc.cell_scales(M, 0.03, 1.0)
    r_g = (sc.V1 / sc.C_KMS) ** 2 * (M + sc.M_BIN)
    assert np.isclose(rp_max, sc.RG_FLOOR_FACTOR * r_g)


def test_b_max_scales_with_rp_convention():
    b10 = sc.b_max_AU(1e6, 1.0, 1.0, rp_max_over_rt=10)
    b5 = sc.b_max_AU(1e6, 1.0, 1.0, rp_max_over_rt=5)
    assert b10 > b5


def test_bb_from_beta_limits():
    # beta at the outer cut (r_p = rp_max) -> b = b_max ; deeper encounters -> b < b_max
    beta_edge = sc.BETA_TIDAL_FLOOR
    assert np.isclose(sc.bb_from_beta(beta_edge, 1e6, 1.0, 1.0), 1.0)
    bb = sc.bb_from_beta(np.array([0.5, 1.0, 5.0]), 1e6, 1.0, 1.0)
    assert np.all(np.diff(bb) < 0) and np.all((bb > 0) & (bb < 1))


def test_frac_from_vinf_snaps_to_nodes():
    a = 1.0
    v_orb = sc.V1 * np.sqrt(sc.M_BIN / a)
    for node in sc.FRAC_NODES:
        assert sc.frac_from_vinf(node * v_orb * 1.05, a) == node
