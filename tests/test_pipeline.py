"""End-to-end check of the bundle tools on a small synthetic bundle (no real data needed)."""
import csv, os, subprocess, sys
import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("scipy")
from conftest import PKG


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    d = tmp_path_factory.mktemp("bundle")
    path = d / "smbh_bhstar_bundle.h5"
    subprocess.run([sys.executable, os.path.join(PKG, "make_test_bundle.py"), str(path),
                    "--n-per-cell", "5", "--seed", "0"], check=True, capture_output=True)
    return path


def test_check_bundle_accepts_synthetic_bundle(bundle):
    out = subprocess.run([sys.executable, os.path.join(PKG, "check_bundle.py"), str(bundle)],
                         check=True, capture_output=True, text=True).stdout
    assert "OK" in out


def test_reduce_cellsummary_conserves_counts(bundle, tmp_path):
    out = tmp_path / "mutde_cellsummary.csv"
    subprocess.run([sys.executable, os.path.join(PKG, "reduce_cellsummary.py"), str(bundle),
                    "--out", str(out)], check=True, capture_output=True)
    with h5py.File(bundle) as f:
        oc = f["results/outcome"].asstr()[:]
    rows = list(csv.DictReader(open(out)))
    assert sum(int(r["N_tot"]) for r in rows) == len(oc)
    assert sum(int(r["N_mu"]) for r in rows) == int(np.sum(oc == "TDE:star2_by_bh1"))
