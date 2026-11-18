"""Smoke tests of the compiled three-body code (skipped if simulation/bin is not built)."""
import subprocess
import pytest

VALID_PREFIXES = ("FLYBY", "EXCHANGE_", "BOUND_TRIPLE", "IONIZATION", "TDE", "PLUNGE", "COLLISION")

STAR_STAR = ["--M_smbh", "1e6", "--m1", "1", "--type1", "star", "--R1", "1", "--m2", "1",
             "--type2", "star", "--R2", "1", "--a_bin", "0.0279028", "--e_bin", "0",
             "--v_inf", "25.2", "--beta", "0.5", "--mode", "classify", "--rtol", "1e-9", "--report", "0"]
BH_STAR = ["--M_smbh", "1e6", "--m1", "10", "--type1", "bh", "--R1", "-1", "--m2", "1",
           "--type2", "star", "--R2", "1", "--a_bin", "0.1", "--e_bin", "0.85", "--outer", "hyperbolic",
           "--beta", "1.5", "--v_inf", "1", "--incl", "3.06878", "--Omega", "0.445981",
           "--omega", "0.347833", "--nu", "3.60568", "--mode", "classify", "--rtol", "1e-10", "--report", "0"]


def run(binary, args):
    out = subprocess.run([binary, *args], check=True, capture_output=True, text=True, timeout=120).stdout
    rows = [l for l in out.splitlines() if l.startswith("CSV,")]
    assert len(rows) == 1, out
    return rows[0].split(",")


def header(binary):
    out = subprocess.run([binary, "--header"], check=True, capture_output=True, text=True).stdout
    return next(l for l in out.splitlines() if l.startswith("CSV,")).split(",")


@pytest.mark.parametrize("args", [STAR_STAR, BH_STAR], ids=["star-star", "bh-star"])
@pytest.mark.parametrize("which", ["bin_pn", "bin_newt"])
def test_encounter_produces_valid_csv(which, args, request):
    binary = request.getfixturevalue(which)
    cols = header(binary)
    row = run(binary, args)
    assert len(row) == len(cols)
    rec = dict(zip(cols, row))
    assert rec["outcome"].startswith(VALID_PREFIXES), rec["outcome"]
    assert float(rec["M_smbh"]) == 1e6


def test_runs_are_deterministic(bin_pn):
    a, b = run(bin_pn, BH_STAR), run(bin_pn, BH_STAR)
    cols = header(bin_pn)
    drop = cols.index("wall_s")
    assert a[:drop] + a[drop + 1:] == b[:drop] + b[drop + 1:]
