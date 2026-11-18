"""Reading SpaceHub trajectory dumps in mutde_lab.py."""
import numpy as np
from conftest import load_module

lab = load_module("mutde_lab", "analysis/debris_mutde/mutde_lab.py")


def test_parse_bss_full_units_and_body_order(tmp_path):
    # SpaceHub DefaultWriter: header, then one row per body per output time, internal
    # time unit = yr/(2 pi) and internal speed unit = 2 pi AU/yr.
    lines = ["time,id,mass,radius,px,py,pz,vx,vy,vz"]
    for k, t in enumerate([0.0, 2 * np.pi]):
        for body in (2, 0, 1):                          # deliberately unsorted ids
            lines.append(",".join(map(str, [t, body, [1e6, 10, 1][body], 0.1,
                                            body + k, 0, 0, 1.0, 0, 0])))
        lines.append("")
    p = tmp_path / "bss_test_full.txt"
    p.write_text("\n".join(lines))
    t, pos, vel = lab.parse_bss_full(str(p))
    assert np.allclose(t, [0.0, 1.0])                   # years
    assert pos.shape == (2, 3, 3)
    assert np.allclose(pos[1, :, 0], [1, 2, 3])         # sorted by id 0,1,2
    assert np.allclose(vel[..., 0], 2 * np.pi)          # AU/yr
