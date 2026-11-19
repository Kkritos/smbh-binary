# Publishing the data on Zenodo (maintainers)

Two Zenodo records accompany the paper:

| record | what | how |
|---|---|---|
| **data** (`upload_type: dataset`) | the ~20 GB of simulation outputs in `data/` | `zenodo_upload.py` in this folder (terminal, via the Zenodo API) |
| **software** (`upload_type: software`) | a snapshot of this repository per release | Zenodo's GitHub integration; metadata from `/.zenodo.json` |

## 1. Data record (terminal)

**One-time setup.** On zenodo.org: *Account → Applications → Personal access tokens → New
token*, with scopes `deposit:write` and `deposit:actions`. Keep the token out of files and
shell history, for example:

```bash
read -s ZENODO_TOKEN && export ZENODO_TOKEN    # paste the token, press Enter
```

**Rehearse on the sandbox first** (a separate site with its own account and token,
https://sandbox.zenodo.org; nothing there is permanent): run every step below with
`--sandbox` and `ZENODO_SANDBOX_TOKEN`.

**Steps.**

```bash
# 0. compress the 14 GB trajectory (the downloader expects the .gz and unpacks it)
gzip -k data/debris_mutde/s1/bss_hyperbolic_bh10_star1_i3.06878_full.txt

# 1. create the draft and reserve the DOI (nothing public yet)
python3 scripts/zenodo/zenodo_upload.py create

# 2. put the reserved DOI / record id into the repository:
#    README.md (two XXXXXXX placeholders), CITATION.cff (identifiers),
#    scripts/get_data.py (ZENODO_RECORD), then commit.

# 3. upload the 8 files (~20 GB; re-run to resume, finished files are skipped)
python3 scripts/zenodo/zenodo_upload.py upload

# 4. check, and look at the draft on the website (the link printed by `status`)
python3 scripts/zenodo/zenodo_upload.py status

# 5. publish -- irreversible: files can no longer be changed (a new version can be made)
python3 scripts/zenodo/zenodo_upload.py publish --yes
```

The file list comes from `scripts/get_data.py` (`GROUPS`), so the uploader and the downloader
cannot disagree. The metadata (title, authors, ORCIDs, description, license CC BY 4.0, link to
the GitHub repository, funding, programming languages) is in `metadata.json`; edit it and run
`metadata` to update the draft. `metadata` overwrites edits made on the website: if you edit
on the website, first run `pull` (saves the draft's metadata to `data/zenodo_pulled*.json`)
and copy the changes into `metadata.json`. On the sandbox, `upload --only <names>` uploads
just a few small files to test the upload path.
After the paper is out, add it under `related_identifiers` (`"relation": "isSupplementTo"`,
the paper's DOI) and update the published record via *New version* or *Edit* on the website.

Before publishing, test the downloader against the draft's file names with the sandbox record:
`ZENODO_RECORD=<sandbox id>` will not work (the downloader uses zenodo.org), so the real check is
right after publishing: `python3 scripts/get_data.py --record <id> starstar`.

## 2. Software record (GitHub integration)

1. Make the GitHub repository public.
2. On zenodo.org: *Account → GitHub*, find the repository and switch it **On**.
3. Create a GitHub release from the tag (`gh release create v1.0.0 --title v1.0.0 --notes-file CHANGELOG.md`
   or on the website). Zenodo archives it within minutes and assigns a DOI, using `/.zenodo.json`.
4. Add that DOI to `CITATION.cff` (`identifiers`) and as a badge in `README.md`.

Each later release (`v1.0.1`, …) gets its own DOI plus a "concept DOI" that always points to
the latest version.
