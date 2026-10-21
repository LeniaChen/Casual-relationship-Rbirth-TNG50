# Does causal discovery recover stellar birth radii? A test of RLCD with TNG50

This repository tests whether the **Rank-based Latent Causal Discovery (RLCD)** approach of Jin et al., *Causal Discovery of Latent Variables in Galactic Archaeology* ([arXiv:2507.00134](https://arxiv.org/abs/2507.00134); code: [ZehaoJin/Causal-Galactic-Archaeology](https://github.com/ZehaoJin/Causal-Galactic-Archaeology)), can recover the **birth radius** of disk stars from six observables — age, [Fe/H], [Mg/H], the vertical action J_z, the angular momentum J_φ and the orbital eccentricity e.

We use 12 Milky-Way-like galaxies from the **TNG50** cosmological simulation, where the true birth radius of every star is known, so the method's latent variables can be compared directly with the truth. The repository contains the full pipeline (TNG50 extraction → orbits with AGAMA → RLCD inputs → RLCD runs and identifiability analysis), all results and figures, a detailed report, and the public input data for the main star selection.

> This is an independent check, not affiliated with the authors of the original paper. The results have not been peer reviewed. Comments and corrections are welcome via GitHub issues.

## Key results

**Applying RLCD to TNG50 the way the paper does (freezing the authors' causal graph and re-estimating its edge coefficients) does not reliably recover the birth radius.**

1. **Structure learning is random.** With identical data and parameters, RLCD's output graph depends on Python's string-hash seed (`PYTHONHASHSEED`). On the authors' own simulated galaxy their published graph appears in only 7 of 30 seeds; the most common result (12/30) is a different graph.
2. **The published edge coefficients are not the maximum-likelihood solution.** The repository's parameter learning uses the absolute value of a determinant and does not require non-negative noise variances, so unphysical solutions can look better; the coefficients used in the paper are far from the optimum (likelihood-ratio statistic ≈ 6.8 × 10⁵ on 119,236 stars).
3. **With the authors' graph frozen, the edge coefficients are not identifiable.** The Jacobian of the map "15 edge coefficients → 15 observed correlations" has rank 14, so a one-parameter family of coefficient sets fits the data equally well; proper maximum-likelihood fits push one observable's noise variance to zero (a Heywood case) and have many nearly equivalent optima. As a result, statistically indistinguishable solutions give correlations between the latent L1 and the true birth radius, |r(L1, Rbirth)|, from about 0 to 0.4–0.5 in most TNG50 galaxies.
4. **Even the most favourable solution is no better than a single observable.** For the paper's selection (7 < R < 10 kpc, |z| < 0.3 kpc), L1 never beats "age alone" or "[Fe/H] alone" as a predictor of the birth radius in any of the 12 galaxies (without observational errors). Widening the selection to the whole disk raises |r(L1, Rbirth)| to ~0.8, but J_φ alone (≈ current radius) already reaches 0.8–0.9 there.
5. **Learning the graph from scratch on TNG50** (Approach B, `PlanB/`) never recovers the authors' graph (0/360 runs) and does not produce a latent that tracks the birth radius better than age.

![Range of |r(L1, Rbirth)| among statistically indistinguishable fits for each galaxy, compared with single-variable baselines](results/rlcd_A/identifiability/summary_ci_nonoise.png)

*For each galaxy (rows), the thick blue bar is the range of |r(L1, true birth radius)| over all fits statistically indistinguishable from the maximum-likelihood fit; triangles and crosses are the correlations of age alone and [Fe/H] alone with the birth radius; open circles are the linear-regression ceiling of the six inputs; the dashed line is the paper's r = 0.64.*

The full analysis, with all numbers, figures and the method explained step by step, is in **[`REPORT.md`](REPORT.md)**.

## What is in this repository

| Path | Content |
|---|---|
| [`REPORT.md`](REPORT.md) | full report: what was done, problems found, results, conclusions, method explained |
| [`STATUS.md`](STATUS.md) | current status, open questions, possible next steps |
| [`docs/REPRODUCTION.md`](docs/REPRODUCTION.md) | detailed reproduction log: every command, data format, decision, installation note and pitfall |
| [`DATA.md`](DATA.md), [`data_manifest.tsv`](data_manifest.tsv) | what data are public, what is not, and how to regenerate everything |
| `scripts/server/` | extraction on the TNG JupyterLab server (`data_obtain.py`, `merge_hdf5.py`) |
| `scripts/local/` | orbits (AGAMA), star selection, Approach A, identifiability analysis, linear ceiling, graph comparison |
| `scripts/reference/` | earlier scripts on which the extraction is based (read-only) |
| `data/release/` | **public data**: RLCD inputs and truth for the main selection, 12 galaxies, with and without observational errors |
| `data/tng50_derived/potentials/` | AGAMA potential coefficients of the 12 galaxies |
| `results/` | logs, checks, RLCD reproduction (`rlcd_g279/`), Approach A and identifiability (`rlcd_A/`), figures (`reproduce_Jin/`, `TNG/`) |
| `star_condition/` | the same analysis for wider star selections (whole disk with three |z| limits) |
| `PlanB/` | Approach B: structure learning from scratch on each galaxy with 30 seeds |
| `tools/` | helper scripts: unpack the public data, check the setup, clone the RLCD code, link the full data |
| `setup/` | requirements of the two Python environments |

Large intermediate data (~17 GB: the TNG50 particle extraction, orbit catalogues and the wide-selection inputs) are not in the repository; see [`DATA.md`](DATA.md).

## Quick start

```bash
git clone https://github.com/<user>/<repo>.git
cd <repo>

# 1. Python environments (tested on macOS arm64)
python3 -m venv .venv && .venv/bin/python -m pip install -r setup/requirements.txt              # analysis (Python 3.10+; tested 3.11)
python3.12 -m venv .venv-rlcd && .venv-rlcd/bin/python -m pip install -r setup/requirements-rlcd.txt   # RLCD (Python 3.12)
# Graphviz is needed to draw causal graphs: brew install graphviz  /  apt install graphviz

# 2. local settings: which Python to use for what
cp env.example .env            # then edit: PY=.venv/bin/python3, PY_RLCD=.venv-rlcd/bin/python3
set -a; . ./.env; set +a       # export PY and PY_RLCD into the current shell

# 3. data and the original RLCD code
bash tools/unpack_data.sh      # public main-selection data → data/rlcd_input/{nonoise,noise}/
bash tools/setup_rlcd_repo.sh  # clones Causal-Galactic-Archaeology at the commit used here (8c0a5a4)
bash tools/check_setup.sh      # checks data, RLCD code, environments, Graphviz
```

All scripts are run **from the repository root** and find their inputs by relative paths.

### Reproduce the main results

```bash
# reproduce the paper's Fig. 3 with the authors' coefficients (seconds)
$PY results/reproduce_Jin/truth_comparison/jin_L1_Rbirth_L2_Rg.py
# reproduce the paper's Fig. 1 by rerunning RLCD structure learning (≈1 min; --seed N to try other hash seeds)
$PY_RLCD results/reproduce_Jin/causal_graph/jin_causal_graph.py

# strict identifiability analysis of the 12 TNG50 galaxies (main result; 1–2 h on 2 cores, --ids 538905 for one galaxy)
$PY scripts/local/rlcd_identifiability_ci.py --version nonoise --extra 200
$PY scripts/local/linear_ceiling.py
# figures for TNG50
$PY results/TNG/causal_graph/tng_causal_graph.py
$PY results/TNG/truth_comparison/tng_L1_Rbirth_L2_Rg.py

# wide selections: runs from the committed covariance matrices, no large data needed
$PY star_condition/R_0_30_Z_0.3/output_data/identifiability_moments.py --version nonoise

# Approach B (structure learning from scratch, 30 seeds per galaxy)
$PY PlanB/input_data/make_input.py
$PY_RLCD PlanB/output_data/rlcd_B.py --version nonoise --seeds 30
$PY PlanB/output_data/rlcd_B_eval.py --version nonoise
$PY PlanB/output_data/summary_planB.py --version nonoise
```

Add `--version noise` to use the version with observational errors. The order of all steps and the expected outputs are listed in [`REPORT.md`](REPORT.md) (Sec. 3, 5.7) and [`docs/REPRODUCTION.md`](docs/REPRODUCTION.md).

### Rebuild everything from the simulation

Requires access to the IllustrisTNG JupyterLab (TNG50-1, snapshot 99) and a trajectory file of the main progenitors; see [`DATA.md`](DATA.md), Sec. 3:
`scripts/server/data_obtain.py` → `merge_hdf5.py` → `scripts/local/compute_orbits.py` (AGAMA, [`docs/INSTALL_AGAMA.md`](docs/INSTALL_AGAMA.md)) → `scripts/local/check_orbits.py` → `scripts/local/make_rlcd_input.py` → analyses above.

## Data

`data/release/rlcd_input_main_{nonoise,noise}.tar.gz` contain, for each of the 12 TNG50 galaxies (Subhalo IDs 487742, 530330, 538905, 545437, 546870, 547293, 552581, 557721, 560421, 562029, 567607, 569251 at z = 0), the stars with 7 < R < 10 kpc and |z| < 0.3 kpc (239,970 stars in total):
- `<ID>_e.csv`: the six RLCD inputs `logJZ, JPhi, ecc, age, FeH, MgH` (same format as the authors' `g2.79e12_e.csv`);
- `<ID>_truth.csv`: the same stars row by row with the true birth radius `R_bir`, guiding radius `Rg`, actions, peri-/apocentres, current position, particle ID `iord`, etc.

The `noise` version adds observational errors as in the paper (10% in age, 0.02 dex in [Fe/H], 0.06 dex in [Mg/Fe]). Column definitions are in [`docs/REPRODUCTION.md`](docs/REPRODUCTION.md), Sec. 6.5.

## Citation and credits

If you use this code or these data, please cite this repository and the works it builds on:
- Jin, Lu, Ting, Zheng & Buck, *Causal Discovery of Latent Variables in Galactic Archaeology*, ApJ, [arXiv:2507.00134](https://arxiv.org/abs/2507.00134) — the method and code being tested;
- Dong et al. (2024), *A Versatile Causal Discovery Framework to Allow Causally-Related Hidden Variables* (ICLR 2024) — the RLCD algorithm;
- the IllustrisTNG / TNG50 simulations: Nelson et al. (2019, *The IllustrisTNG simulations: public data release*), Pillepich et al. (2019) and Nelson et al. (2019) (TNG50);
- Vasiliev (2019), AGAMA.

## License

Code: MIT (see [`LICENSE`](LICENSE)). Data products in `data/` and `results/`: CC BY 4.0; they are derived from the public IllustrisTNG data, whose use should also be acknowledged as requested by the IllustrisTNG collaboration. The original RLCD code is not included and keeps its own license.
