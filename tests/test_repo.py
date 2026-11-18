"""Repository consistency: version strings, notebooks, data manifest, no machine-specific paths."""
import json, os, re, subprocess, sys
import pytest
from conftest import ROOT, load_module

TRACKED = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True).stdout.split()


def test_version_consistent():
    version = open(os.path.join(ROOT, "VERSION")).read().strip()
    assert re.fullmatch(r"\d+\.\d+\.\d+", version)
    assert f"version: {version}" in open(os.path.join(ROOT, "CITATION.cff")).read()
    assert f"[{version}]" in open(os.path.join(ROOT, "CHANGELOG.md")).read()


@pytest.mark.parametrize("nb", [f for f in TRACKED if f.endswith(".ipynb")])
def test_notebooks_valid_and_without_saved_errors(nb):
    data = json.load(open(os.path.join(ROOT, nb)))
    errors = [o for c in data["cells"] if c["cell_type"] == "code"
              for o in c.get("outputs", []) if o.get("output_type") == "error"]
    assert not errors


def test_data_manifest_destinations_unique():
    gd = load_module("get_data", "scripts/get_data.py")
    dests = [rel for _, files in gd.GROUPS.values() for _, rel in files]
    assert len(dests) == len(set(dests))
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "get_data.py"), "--list"], check=True,
                   capture_output=True)


def test_no_machine_specific_paths():
    pattern = re.compile("|".join(["/Us" + "ers/", "/ho" + "me/", "/gp" + "fs/"]))
    hits = []
    for f in TRACKED:
        if f.startswith("tests/") or f.endswith((".png", ".pdf", ".mp4", ".npz")):
            continue
        try:
            text = open(os.path.join(ROOT, f), encoding="utf-8").read()
        except (UnicodeDecodeError, FileNotFoundError):
            continue
        hits += [f for _ in pattern.finditer(text)][:1]
    assert not hits, f"absolute machine paths in: {sorted(set(hits))}"
