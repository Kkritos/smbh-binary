#!/usr/bin/env python3
"""
bundle_to_h5.py
---------------
Pack the BH-star campaign into ONE compressed HDF5 file using h5py (NO PyTables):
    /results/<column>    - the full simulation table, one dataset per column
    /grid_meta/<column>  - b_max per (M_smbh, a_bin, frac) grid point
    /ecc_meta/<column>   - e_max(a) survival limit per a_bin

Column-oriented layout: each CSV column becomes its own 1-D dataset (float64, or
variable-length UTF-8 for the string columns), gzip-compressed.  The results file
is streamed in chunks so peak memory stays ~0.5 GB.  Needs only pandas (to parse
the CSV) + h5py.  The three CSVs are left untouched.

Run from the folder holding the three CSVs:
    python3 bundle_to_h5.py

Read back (h5py only, no pandas needed):
    import h5py
    f = h5py.File("smbh_bhstar_bundle.h5", "r")
    M       = f["results/M_smbh"][:]            # float array
    outcome = f["results/outcome"].asstr()[:]   # string array
    bmax    = f["grid_meta/b_max_AU"][:]
    emax    = f["ecc_meta/e_max_survival"][:]

Rebuild a DataFrame if you want one:
    import h5py, pandas as pd
    f = h5py.File("smbh_bhstar_bundle.h5", "r")
    STR = {"type1","type2","outer","outcome","exit_mode"}
    df = pd.DataFrame({c: (f[f"results/{c}"].asstr()[:] if c in STR else f[f"results/{c}"][:])
                       for c in f["results"].keys()})
"""

import sys
import pandas as pd

try:
    import h5py
except ImportError:
    sys.exit("ERROR: h5py missing.  Install with:  pip install --user h5py   (or: conda install h5py)")

RES = "smbh_bhstar_results.csv"
GRID = "grid_meta.csv"
ECC = "ecc_meta.csv"
OUT = "smbh_bhstar_bundle.h5"
CHUNK = 1_000_000

STR_COLS = {"type1", "type2", "outer", "outcome", "exit_mode"}
VSTR = h5py.string_dtype(encoding="utf-8")       # variable-length UTF-8
GZIP = dict(compression="gzip", compression_opts=4)


def write_small(h5, group, df):
    """Write a small table (metadata) whole: one dataset per column."""
    g = h5.create_group(group)
    for col in df.columns:
        if col in STR_COLS or df[col].dtype == object:
            g.create_dataset(col, data=df[col].astype(str).to_numpy(), dtype=VSTR, **GZIP)
        else:
            g.create_dataset(col, data=df[col].to_numpy("float64"), **GZIP)


def main():
    with h5py.File(OUT, "w") as h5:
        # --- metadata tables ---
        write_small(h5, "grid_meta", pd.read_csv(GRID))
        write_small(h5, "ecc_meta", pd.read_csv(ECC))
        sys.stderr.write("wrote /grid_meta and /ecc_meta\n")

        # --- results table, streamed into resizable datasets ---
        g = h5.create_group("results")
        dsets = {}
        cols = None
        n = 0
        for chunk in pd.read_csv(RES, chunksize=CHUNK):
            chunk = chunk.drop(columns=[c for c in ("CSV",) if c in chunk.columns])
            if cols is None:                         # first chunk: create the datasets
                cols = list(chunk.columns)
                for col in cols:
                    is_str = col in STR_COLS or chunk[col].dtype == object
                    dsets[col] = g.create_dataset(
                        col, shape=(0,), maxshape=(None,),
                        dtype=(VSTR if is_str else "float64"),
                        chunks=(CHUNK,), **GZIP)
            m = len(chunk)
            for col in cols:
                d = dsets[col]
                d.resize((n + m,))
                if col in STR_COLS:
                    d[n:n + m] = chunk[col].astype(str).to_numpy()
                else:
                    d[n:n + m] = chunk[col].to_numpy("float64")
            n += m
            sys.stderr.write(f"\r  results: {n:,} rows")
        sys.stderr.write("\n")
        h5["results"].attrs["n_rows"] = n
        h5.attrs["source"] = "smbh_bhstar campaign (BH-star + SMBH, e-grid)"

    # --- verify ---
    with h5py.File(OUT, "r") as h5:
        nr = h5["results"].attrs["n_rows"]
        ncols = len(h5["results"].keys())
        ng = len(h5["grid_meta/b_max_AU"]) if "grid_meta/b_max_AU" in h5 else "?"
        ne = len(h5["ecc_meta/e_max_survival"]) if "ecc_meta/e_max_survival" in h5 else "?"
    sys.stderr.write(f"DONE -> {OUT}  (results={nr:,} rows x {ncols} cols, "
                     f"grid_meta={ng}, ecc_meta={ne})\n")


if __name__ == "__main__":
    main()
