#!/usr/bin/env python3
"""Download the large data files from Zenodo into data/ (standard library only).

    python3 scripts/get_data.py --list              # what exists and what each group is for
    python3 scripts/get_data.py starstar bhstar     # download some groups
    python3 scripts/get_data.py all                 # everything (~20 GB download)

Each file is checked against the MD5 that Zenodo publishes, and `.gz` files are
decompressed in place. Files already present in data/ are skipped.
Set the record with --record or the ZENODO_RECORD environment variable.
"""
import argparse, gzip, hashlib, json, os, shutil, sys, urllib.request

ZENODO_RECORD = os.environ.get("ZENODO_RECORD", "24372716")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.environ.get("SMBH_DATA", os.path.join(ROOT, "data"))

# group -> (what it is for, [(file name on Zenodo, destination relative to data/)])
GROUPS = {
    "starstar": ("star-star collision CSVs (Newt + PN); analysis/yulai2024/comparison_plots.ipynb, ~130 MB", [
        ("smbh_starstar_collisions_Newt_results.csv", "starstar/smbh_starstar_collisions_Newt_results.csv"),
        ("smbh_starstar_collisions_PN_results.csv",   "starstar/smbh_starstar_collisions_PN_results.csv"),
    ]),
    "bhstar": ("BH-star campaign 2 (rp_max = 10 r_t) bundle + delay tables; all bhstar_mutde notebooks, ~10 GB", [
        ("smbh_bhstar_bundle.h5",  "bhstar/campaign2_rp10/smbh_bhstar_bundle.h5"),
        ("delay_per_event.csv",    "bhstar/campaign2_rp10/delay_per_event.csv"),
        ("delay_cellsummary.csv",  "bhstar/campaign2_rp10/delay_cellsummary.csv"),
    ]),
    "bhstar_rp5": ("BH-star campaign 1 (rp_max = 5 r_t) bundle; only the rp5-vs-rp10 appendix of analyze_branching, ~5 GB", [
        ("smbh_bhstar_bundle_750.h5", "bhstar/campaign1_rp5/smbh_bhstar_bundle_750.h5"),
    ]),
    "debris": ("SpaceHub trajectory of the single micro-TDE encounter; mutde_lab.py --traj (movie + final still), "
               "~14 GB uncompressed", [
        ("bss_hyperbolic_bh10_star1_i3.06878_full.txt.gz",
         "debris_mutde/s1/bss_hyperbolic_bh10_star1_i3.06878_full.txt.gz"),
    ]),
}


def zenodo_checksums(record):
    """{file name: md5} from the Zenodo records API."""
    with urllib.request.urlopen(f"https://zenodo.org/api/records/{record}") as r:
        meta = json.load(r)
    out = {}
    for f in meta.get("files", []):
        algo, _, value = f.get("checksum", "").partition(":")
        if algo == "md5":
            out[f["key"]] = value
    return out


def md5sum(path, chunk=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(url, dest):
    tmp = dest + ".part"
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as f:
        total = int(r.headers.get("Content-Length", 0)); done = 0
        while True:
            block = r.read(1 << 22)
            if not block:
                break
            f.write(block); done += len(block)
            if total:
                print(f"\r    {done/1e9:6.2f} / {total/1e9:.2f} GB", end="", flush=True)
    print()
    os.replace(tmp, dest)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("groups", nargs="*", help="groups to fetch (see --list), or 'all'")
    ap.add_argument("--list", action="store_true", help="list the groups and whether they are present")
    ap.add_argument("--record", default=ZENODO_RECORD, help="Zenodo record id")
    a = ap.parse_args()

    if a.list or not a.groups:
        for g, (desc, files) in GROUPS.items():
            print(f"{g:11s} {desc}")
            for _, rel in files:
                path = os.path.join(DATA, rel[:-3] if rel.endswith(".gz") else rel)
                print(f"{'':11s}   [{'x' if os.path.exists(path) else ' '}] data/{os.path.relpath(path, DATA)}")
        return

    groups = list(GROUPS) if a.groups == ["all"] else a.groups
    bad = [g for g in groups if g not in GROUPS]
    if bad:
        sys.exit(f"unknown group(s): {bad}; choose from {list(GROUPS)} or 'all'")
    if not a.record:
        sys.exit("No Zenodo record id: pass --record ID or set ZENODO_RECORD.")

    sums = zenodo_checksums(a.record)
    for g in groups:
        for name, rel in GROUPS[g][1]:
            dest = os.path.join(DATA, rel)
            final = dest[:-3] if dest.endswith(".gz") else dest
            if os.path.exists(final):
                print(f"[skip] data/{os.path.relpath(final, DATA)} already present"); continue
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            print(f"[get ] {name} -> data/{rel}")
            download(f"https://zenodo.org/records/{a.record}/files/{name}?download=1", dest)
            if name in sums and md5sum(dest) != sums[name]:
                os.remove(dest); sys.exit(f"MD5 mismatch for {name}; removed it, please retry.")
            if dest.endswith(".gz"):
                print(f"[gunz] {os.path.basename(final)}")
                with gzip.open(dest, "rb") as fi, open(final, "wb") as fo:
                    shutil.copyfileobj(fi, fo, 1 << 24)
                os.remove(dest)
    print("done.")


if __name__ == "__main__":
    main()
