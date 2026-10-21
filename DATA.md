# Data

## 1. What is public

| Path | Content | Size |
|---|---|---|
| `data/release/rlcd_input_main_nonoise.tar.gz` | main selection (7 < R < 10 kpc, \|z\| < 0.3 kpc), 12 TNG50 galaxies, no observational errors: `<ID>_e.csv` (6 RLCD inputs) + `<ID>_truth.csv` (truth) | 51 MB (113 MB unpacked) |
| `data/release/rlcd_input_main_noise.tar.gz` | same stars with observational errors (age 10%, [Fe/H] 0.02 dex, [Mg/Fe] 0.06 dex) | 61 MB (132 MB unpacked) |
| `data/tng50_derived/potentials/<ID>_{halo,disk}.coef` | AGAMA potential expansions of the 12 galaxies (`agama.Potential(file=...)`) | 0.4 MB |
| `star_condition/*/output_data/moments_*.json` | 8×8 covariance matrices (6 inputs + true birth radius + true guiding radius) of every galaxy for the wide selections; enough to rerun the wide-selection analysis | small |
| `results/` | all summary tables, logs and figures | ~60 MB |

Unpack with `bash tools/unpack_data.sh` (creates `data/rlcd_input/{nonoise,noise}/` and verifies checksums).

**Noise model**: age × (1 + 0.10·N), [Fe/H] + 0.02·N, [Mg/Fe] + 0.06·N.
- Meaning: each star's age gets a 10% relative error, [Fe/H] an absolute error of 0.02 dex and [Mg/Fe] one of 0.06 dex.
- N: an independent standard-normal random number (mean 0, standard deviation 1) for every star and quantity; the seed is 20260926 + galaxy ID, so regeneration gives identical errors.

The selection uses the error-free values, so both versions contain the same stars, row by row. Column definitions: [`docs/REPRODUCTION.md`](docs/REPRODUCTION.md), Sec. 6.5.

## 2. What is not public (~17 GB) and how it was made

These files are too large for GitHub. `data_manifest.tsv` lists every one of them (path, bytes, MD5 checksum — a fingerprint that changes if a single byte changes), so a regenerated or copied set can be verified.

| Path | Files | Size | Produced by | Regenerable? |
|---|---|---|---|---|
| `data/tng50_extract/data_obtain_all.hdf5` | 1 | 955 MB | `scripts/server/data_obtain.py` + `merge_hdf5.py` on the TNG JupyterLab | needs TNG server access and a main-progenitor trajectory file (several hours) |
| `data/tng50_derived/orbits_all.hdf5` | 1 | 271 MB | `scripts/local/compute_orbits.py` (AGAMA, ~1 h for 12 galaxies) | yes, from the file above |
| `data/tng50_derived/orbits_test.hdf5` | 1 | 0.2 MB | `compute_orbits.py --test` (old setting) | not needed |
| `PlanB/input_data/{nonoise,noise}/`, `PlanB/true_value/` | 48 | 245 MB | `PlanB/input_data/make_input.py` (copy of `data/rlcd_input/`) | yes, from the public data |
| `star_condition/R_0_30_Z_30/{input_data/*,true_value}/` | 48 | 6.3 GB | `star_condition/R_0_30_Z_30/input_data/make_input.py` (reads the two hdf5 files) | yes |
| `star_condition/R_0_30_Z_3/...` | 48 | 6.0 GB | its `make_input.py` (cuts R_0_30_Z_30 at \|z\| < 3 kpc) | yes, from R_0_30_Z_30 |
| `star_condition/R_0_30_Z_0.3/...` | 48 | 2.8 GB | its `make_input.py` (cuts R_0_30_Z_30 at \|z\| < 0.3 kpc) | yes, from R_0_30_Z_30 |
| `Causal-Galactic-Archaeology/` | — | 257 MB | git clone of the authors' repository (commit 8c0a5a4, unmodified) | `bash tools/setup_rlcd_repo.sh` |

`data/rlcd_input/` itself is listed in the manifest too, because it is what `tools/unpack_data.sh` checks.

If you have a folder with the full data in the same layout as this repository, `bash tools/link_data.sh <folder>` links it into the repository with symbolic links (nothing is copied), and `bash tools/check_setup.sh --md5` verifies it.

## 3. Rebuilding everything from the simulation
```bash
# 0) on the TNG JupyterLab (with its old AnastrisTNG); details in docs/REPRODUCTION.md Sec. 6.1
cd <working directory>
ln -s <path>/TNG50_disk_trajectory.hdf5 .            # once
nohup python3 data_obtain.py > data_obtain.log 2>&1 &
python3 merge_hdf5.py                                 # → data_obtain_all.hdf5; download into data/tng50_extract/

# 1) locally, from the repository root
$PY scripts/local/compute_orbits.py 2>&1 | tee results/compute_orbits.log   # needs AGAMA, ~1 h
$PY scripts/local/check_orbits.py
$PY scripts/local/make_rlcd_input.py
$PY PlanB/input_data/make_input.py
$PY star_condition/R_0_30_Z_30/input_data/make_input.py
$PY star_condition/R_0_30_Z_3/input_data/make_input.py
$PY star_condition/R_0_30_Z_0.3/input_data/make_input.py

# 2) compare with the original fingerprints
bash tools/check_setup.sh --md5
```
`tee` shows the output on screen while writing it to the log. Regenerated files are not necessarily byte-identical to the originals (HDF5 files store metadata, floating-point summation order can differ), in which case `--md5` reports mismatches; compare the numbers instead, e.g. rerun `scripts/local/check_orbits.py` and compare the star counts in `results/rlcd_input_selection.txt`.

The **main-progenitor trajectory file** (`TNG50_disk_trajectory.hdf5`: centre, velocity and angular-momentum direction of each galaxy's main progenitor at every snapshot) is not public; `scripts/reference/TNG50_disk_trajectory.py` shows how it was built on the TNG server.
