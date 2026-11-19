#!/usr/bin/env python3
"""Create, fill and publish the Zenodo record for the simulation data, from the terminal.

Steps (see scripts/zenodo/README.md):

    export ZENODO_TOKEN=...                       # personal access token (deposit:write, deposit:actions)
    python3 scripts/zenodo/zenodo_upload.py create    # new draft + reserved DOI
    python3 scripts/zenodo/zenodo_upload.py upload    # upload the data files (resumable: skips files already there)
    python3 scripts/zenodo/zenodo_upload.py metadata  # (re)apply metadata.json
    python3 scripts/zenodo/zenodo_upload.py status    # what is in the draft
    python3 scripts/zenodo/zenodo_upload.py publish --yes   # irreversible: the record becomes public

Add --sandbox to every command to rehearse on https://sandbox.zenodo.org (separate account
and token, ZENODO_SANDBOX_TOKEN). The draft id is remembered in data/zenodo_draft[_sandbox].json.
Only the standard library and curl are needed.
"""
import argparse, json, os, subprocess, sys, tempfile, urllib.error, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from get_data import GROUPS, DATA          # noqa: E402  (single source of the file list)


def api_base(sandbox):
    return "https://sandbox.zenodo.org/api" if sandbox else "https://zenodo.org/api"


def token(sandbox):
    var = "ZENODO_SANDBOX_TOKEN" if sandbox else "ZENODO_TOKEN"
    tok = os.environ.get(var)
    if not tok:
        sys.exit(f"Set {var} to a Zenodo personal access token with scopes deposit:write and deposit:actions.")
    return tok


def state_path(sandbox):
    return os.path.join(DATA, "zenodo_draft_sandbox.json" if sandbox else "zenodo_draft.json")


def call(method, url, tok, body=None):
    req = urllib.request.Request(url, method=method, headers={"Authorization": f"Bearer {tok}",
                                                              "Content-Type": "application/json"})
    data = json.dumps(body).encode() if body is not None else None
    try:
        with urllib.request.urlopen(req, data=data) as r:
            txt = r.read().decode()
            return json.loads(txt) if txt else {}
    except urllib.error.HTTPError as e:
        sys.exit(f"{method} {url} -> HTTP {e.code}: {e.read().decode()[:2000]}")


def load_draft(sandbox):
    p = state_path(sandbox)
    if not os.path.exists(p):
        sys.exit(f"No draft yet ({p} missing): run `create` first.")
    return json.load(open(p))


def metadata():
    return json.load(open(os.path.join(HERE, "metadata.json")))["metadata"]


def cmd_create(a, tok):
    if os.path.exists(state_path(a.sandbox)) and not a.force:
        sys.exit(f"A draft is already recorded in {state_path(a.sandbox)} (use --force to start another).")
    meta = metadata()
    meta["prereserve_doi"] = True
    dep = call("POST", f"{api_base(a.sandbox)}/deposit/depositions", tok, {"metadata": meta})
    doi = dep["metadata"].get("prereserve_doi", {}).get("doi")
    os.makedirs(DATA, exist_ok=True)
    json.dump({"id": dep["id"], "bucket": dep["links"]["bucket"], "html": dep["links"]["html"],
               "doi": doi}, open(state_path(a.sandbox), "w"), indent=2)
    print(f"draft {dep['id']} created: {dep['links']['html']}")
    print(f"reserved DOI: {doi}   (record id {dep['id']}; the DOI becomes active only on publish)")
    print("Put the DOI and record id into README.md, CITATION.cff and scripts/get_data.py now.")


def local_files():
    """(name on Zenodo, local path) for every file in scripts/get_data.py's GROUPS."""
    out = []
    for _, files in GROUPS.values():
        for name, rel in files:
            out.append((name, os.path.join(DATA, rel)))
    return out


def cmd_upload(a, tok):
    d = load_draft(a.sandbox)
    have = {f["filename"]: f["filesize"]
            for f in call("GET", f"{api_base(a.sandbox)}/deposit/depositions/{d['id']}/files", tok)}
    todo = local_files()
    if a.only:
        unknown = set(a.only) - {n for n, _ in todo}
        if unknown:
            sys.exit(f"--only: not in the file list: {sorted(unknown)}")
        todo = [(n, p) for n, p in todo if n in a.only]
    missing = [p for _, p in todo if not os.path.exists(p)]
    if missing:
        msg = "\n  ".join(missing)
        sys.exit(f"Missing local files:\n  {msg}\n(.gz files: create them with `gzip -k <file>` first.)")
    with tempfile.NamedTemporaryFile("w", delete=False) as h:     # keep the token off the command line
        h.write(f"Authorization: Bearer {tok}\n"); hdr = h.name
    os.chmod(hdr, 0o600)
    try:
        for name, path in todo:
            size = os.path.getsize(path)
            if have.get(name) == size:
                print(f"[skip] {name} already uploaded"); continue
            print(f"[put ] {name}  ({size/1e9:.2f} GB)")
            subprocess.run(["curl", "--fail", "--progress-bar", "-H", f"@{hdr}", "--upload-file", path,
                            f"{d['bucket']}/{name}", "-o", os.devnull], check=True)
    finally:
        os.remove(hdr)
    print("upload complete; check with `status`.")


def cmd_metadata(a, tok):
    d = load_draft(a.sandbox)
    meta = metadata()
    meta["prereserve_doi"] = True          # keep the DOI reserved at `create`
    call("PUT", f"{api_base(a.sandbox)}/deposit/depositions/{d['id']}", tok, {"metadata": meta})
    print("metadata updated from scripts/zenodo/metadata.json")


def cmd_status(a, tok):
    d = load_draft(a.sandbox)
    dep = call("GET", f"{api_base(a.sandbox)}/deposit/depositions/{d['id']}", tok)
    print(f"record {dep['id']}  state={dep['state']}  submitted={dep['submitted']}")
    print(f"title: {dep['metadata'].get('title')}")
    print(f"DOI:   {d.get('doi')}\nweb:   {dep['links']['html']}")
    total = 0
    for f in dep.get("files", []):
        total += f["filesize"]
        print(f"  {f['filesize']/1e9:8.3f} GB  md5:{f['checksum']}  {f['filename']}")
    print(f"  {total/1e9:8.3f} GB total")
    want = {n for n, _ in local_files()}
    gone = want - {f["filename"] for f in dep.get("files", [])}
    if gone:
        print(f"not yet uploaded: {sorted(gone)}")


def cmd_pull(a, tok):
    """Save the draft's current metadata (including edits made on the website) for review."""
    d = load_draft(a.sandbox)
    tag = "_sandbox" if a.sandbox else ""
    legacy = call("GET", f"{api_base(a.sandbox)}/deposit/depositions/{d['id']}", tok)
    rdm = call("GET", f"{api_base(a.sandbox)}/records/{d['id']}/draft", tok)
    for kind, obj in (("legacy", legacy), ("rdm", rdm)):
        p = os.path.join(DATA, f"zenodo_pulled{tag}_{kind}.json")
        json.dump(obj, open(p, "w"), indent=2)
        print(f"saved {p}")


def cmd_publish(a, tok):
    if not a.yes:
        sys.exit("Publishing is irreversible (files can no longer change). Re-run with --yes.")
    d = load_draft(a.sandbox)
    dep = call("POST", f"{api_base(a.sandbox)}/deposit/depositions/{d['id']}/actions/publish", tok)
    print(f"published: {dep['links'].get('record_html', dep['links']['html'])}  DOI {dep.get('doi')}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["create", "upload", "metadata", "status", "pull", "publish"])
    ap.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org (rehearsal)")
    ap.add_argument("--force", action="store_true", help="create: start a new draft even if one is recorded")
    ap.add_argument("--yes", action="store_true", help="publish: confirm")
    ap.add_argument("--only", nargs="+", metavar="NAME", help="upload: only these file names (e.g. a small sandbox test)")
    a = ap.parse_args()
    tok = token(a.sandbox)
    {"create": cmd_create, "upload": cmd_upload, "metadata": cmd_metadata,
     "status": cmd_status, "pull": cmd_pull, "publish": cmd_publish}[a.command](a, tok)


if __name__ == "__main__":
    main()
