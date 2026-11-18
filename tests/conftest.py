"""Shared fixtures: repository paths and the compiled simulation binaries."""
import importlib.util, os, shutil, sys
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, "analysis", "bhstar_mutde", "pkg")
sys.path.insert(0, PKG)


def load_module(name, relpath):
    """Import a script that is not on the path (e.g. analysis/debris_mutde/mutde_lab.py)."""
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, relpath))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _binary(name):
    path = os.path.join(os.environ.get("SMBH_BIN", os.path.join(ROOT, "simulation", "bin")), name)
    if not os.access(path, os.X_OK):
        pytest.skip(f"{name} not built (run `make -C simulation`)")
    return path


@pytest.fixture(scope="session")
def bin_pn():
    return _binary("smbh_binary_PN")


@pytest.fixture(scope="session")
def bin_newt():
    return _binary("smbh_binary_Newt")
