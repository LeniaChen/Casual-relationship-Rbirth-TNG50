# Reproduction log and technical notes

This is the detailed working document of the project: every step, command, data format, decision, installation note and pitfall, written so that someone with no prior context (a person or an AI assistant) can reproduce the current state and continue. For the scientific story read [`REPORT.md`](../REPORT.md); for the current status and next steps read [`STATUS.md`](../STATUS.md).

- All paths are relative to the repository root. `$PY` is the Python of the analysis environment and `$PY_RLCD` that of the RLCD environment (see [`env.example`](../env.example) and the README).
- Only the main-selection data are public (`data/release/`, unpack with `bash tools/unpack_data.sh`); everything else must be regenerated as described in [`DATA.md`](../DATA.md).
- Last research update: 2026-09-27.

---

## ⚠ Problems found (reported neither in the paper nor in the code repository; read before using any result)
1. **Structure learning is random**: with the same α, different PYTHONHASHSEED values give different causal graphs; the authors' graph appears in 1 of 5 runs (seed = 1) and 7 of 30. Causes: iteration order of sets of variable names (depends on string-hash randomization) and an unseeded `random.sample` in `LatentGroups.py`. → Sec. 6.3, `results/rlcd_g279/`
2. **Selection of learned parameters is flawed**: the repository's `parameter_by_trek` computes the NLL with the absolute value of the determinant and does not check that implied noise variances are ≥ 0, so physically impossible solutions can get a "better" NLL. The coefficients used in the paper are the set with the second smallest NLL (1.586), not the smallest (1.498), without explanation, and are not converged (optimum ≈ −1.27). → Sec. 6.6
3. **(Most fundamental) with the authors' graph frozen, the edge coefficients cannot be determined reliably**:
   - Jacobian rank 14 < 15 (15 edges, 15 observed correlations): a family of coefficient sets fits exactly equally well (mainly L1→age trading off against L1→MgH). A property of the graph, independent of the data.
   - A proper maximum-likelihood fit (noise variances ≥ 0, strict convergence, 230 starts) almost always pushes the noise variance of one observable to 0 (22 of 24 TNG50 data sets, mostly J_φ; [Fe/H] or e for g2.79e12; Heywood case), with many local optima within ~1e-4 in NLL.
   - Consequence: in most TNG50 galaxies, solutions statistically indistinguishable from the best (2·N·ΔNLL < 3.84) give |r(L1, Rbirth)| from ~0 to 0.4–0.5; on g2.79e12, going from 30 to 230 starts changes the best solution's |r(L1, Rb)| from 0.25 to 0.62. **Even the most favourable solution essentially never beats age alone or [Fe/H] alone** (nonoise 12/12; noise 2 marginal exceptions: 0.58 vs 0.57, 0.71 vs 0.67).
   - The paper's coefficients (NLL 1.593 with the correct formula) are far from the optimum (≈ −1.272): likelihood-ratio statistic ≈ 680,000 for 119,000 stars.
   - **Correction**: the first summary (`summary_<version>.*`) used "ΔNLL < 0.01" as "equally good", which is too loose (likelihood ratio < 0.02·N); use the strict `summary_ci_<version>.*`.
   → end of Sec. 6.6, `results/rlcd_A/identifiability/`, `results/rlcd_A/README.md`, `REPORT.md` Secs. 4.3, 5.6

**Decision (2026-09-26)**: report the non-identifiability itself as the main result: for the 12 galaxies give the possible range of |r(L1, Rbirth)| and |r(L2, Rg)| among equally good fits, and flag these problems clearly in all results and documents.

## 0. Before starting (for a person or AI assistant picking up the project)
1. Read [`STATUS.md`](../STATUS.md), then Sec. 2 (progress) and Sec. 9 (pitfalls) of this file.
2. Set up the environments and data: `cp env.example .env` (set `PY`, `PY_RLCD`) → `bash tools/unpack_data.sh` (public main-selection data) → `bash tools/setup_rlcd_repo.sh` → `bash tools/check_setup.sh`. With the full data: `bash tools/link_data.sh <data folder>`.
3. Check the latest state: the logs and reports in `results/` (for `results/compute_orbits.log`, replace the progress-bar carriage returns first: `tr '\r' '\n' < results/compute_orbits.log | grep -v "orbit integration"`).
4. Follow the conventions in Sec. 1. After any progress, update Sec. 2 of this file and `STATUS.md`.

---

## 1. Conventions
- Documentation: every formula is followed by a plain-language explanation and the meaning of each symbol.
- Everything (scripts, results) lives inside this repository at the paths listed in Sec. 4; large data are never committed (see [`DATA.md`](../DATA.md) and `.gitignore`).
- Two Python environments (Sec. 5): the analysis environment (`$PY`, Python 3.11, [`setup/requirements.txt`](../setup/requirements.txt), AGAMA optional) and the RLCD environment (`$PY_RLCD`, Python 3.12, [`setup/requirements-rlcd.txt`](../setup/requirements-rlcd.txt)). Install packages into the matching environment, not into a base/system Python; write the full interpreter path in commands if the default `python3` is a different one.
- Particle data are downloaded and processed locally; the TNG server is used only for the extraction.
- Server-side extraction runs on the **TNG JupyterLab**, which has an **old** AnastrisTNG version: do not write server scripts against the 2.0.0 API.
- Long runs write their output to a log in `results/` (`2>&1 | tee results/<name>.log`).
- RLCD: always fix `PYTHONHASHSEED` and repeat with several seeds; coefficient fits: noise variances ≥ 0, converging optimizer, many starts, report all statistically indistinguishable solutions, never select using the truth; always compare with single-variable baselines.

---

## 2. Progress

### Stage 1: prepare TNG50 data (6 inputs + 2 truths) — done 2026-09-26
- [x] 1.1 read the paper and the repository, fix the input/output formats (Sec. 3)
- [x] 1.2 server extraction of 12 galaxies: `scripts/server/data_obtain.py`, merged with `merge_hdf5.py` → downloaded as `data/tng50_extract/data_obtain_all.hdf5` (0.96 GB)
- [x] 1.3 star-by-star cross-check against the earlier `001_orbit` results: consistent (`scripts/local/check_vs_001_orbit.py`, report `results/check_vs_001_orbit.txt`)
- [x] 1.4a AGAMA 1.0 installed in the analysis environment (Homebrew libomp for OpenMP, [`docs/INSTALL_AGAMA.md`](INSTALL_AGAMA.md))
- [x] 1.4b `compute_orbits.py --test` passed (5000 stars of 487742, all checks passed)
- [x] 1.4c `compute_orbits.py` done for all 12 galaxies (3.8–8.0 min per galaxy, no errors, log `results/compute_orbits.log`) → `data/tng50_derived/orbits_all.hdf5`
  - the first run failed for 530330 (ActionFinder initialization); fixed by setting the inner grid boundary of the potential to rmin = 0.2 kpc, version `v2-rmin-grid` (Sec. 6.2)
- [x] 1.5 orbit results of the 12 galaxies checked: all pass (`scripts/local/check_orbits.py` → `results/orbit_check/`, Sec. 6.4)
- [x] 1.6 star selection and RLCD input csv files (`scripts/local/make_rlcd_input.py` → `data/rlcd_input/{nonoise,noise}/<ID>_e.csv` and `<ID>_truth.csv`, Sec. 6.5)

### Stage 2: run RLCD and compare with the truth
- [x] 2.1 Causal-Galactic-Archaeology installed: separate Python 3.12 environment, dependencies `setup/requirements-rlcd.txt`, repository in `Causal-Galactic-Archaeology/` (Sec. 8.3)
- [x] 2.2 reproduce the paper with the repository's g2.79e12 data: 4 α combinations, comparison in Sec. 6.3 and `results/rlcd_g279/compare_4alpha.txt`. Installation correct; with α0 = 1e-5 the latent structure of the paper is found (all 15 of the authors' edges plus 3 more); stricter α0 gives 4 latents. None identical to the authors' graph
- [x] 2.2a cause of the difference: **RLCD results are random**. Same α, 5 PYTHONHASHSEED values: only seed = 1 is identical to the authors' graph (Sec. 6.3)
- [x] 2.2b decision (2026-09-26): **Approach A main result, Approach B supplementary**
  - A: freeze the authors' g2.79 graph and re-estimate only the edge coefficients on each TNG50 galaxy (as in paper Secs. 4.2, 4.3)
  - B: learn the graph from scratch on each TNG50 galaxy with many seeds and count how often each structure appears
- [x] **2.3A Approach A — identifiability analysis of the 12 galaxies (done 2026-09-26)**, results in `results/rlcd_A/identifiability/`:
  - round 1 (loose): `rlcd_identifiability.py` (per galaxy: Jacobian rank, 30 starts, L1→age profile) + `rlcd_identifiability_summary.py` → `summary_{nonoise,noise}.*`
  - linear ceiling: `linear_ceiling.py` → `linear_ceiling.*`
  - round 2 (strict, main result): `rlcd_identifiability_ci.py --extra 200` (230 starts, noise variances ≥ 0, strict convergence, 95% likelihood ratio) → `<tag>_ci.*`, `summary_ci_{nonoise,noise}.*`
  - run on a 2-core Linux machine; commands in `REPORT.md` Sec. 5.7
- [x] 2.4 (Approach A) L1 vs Rbirth and L2 vs Rg → above; **conclusion: not viable as done in the paper** (`REPORT.md` Sec. 6)
- [x] Report `REPORT.md`. Updated 2026-09-27: Sec. 0, Sec. 3.7 (30 seeds: authors' graph 7/30), Sec. 3.11 figure scripts, Sec. 3.12 wider selections (covariance method), Secs. 4.1, 4.3, Sec. 5.8 (why the f(L1) bands are narrow: σ_f = |r|·σ_R), Sec. 5.9 (tables and interpretation for the three wide selections), Sec. 6.1 items 4 and 6, Sec. 7, Sec. 8.4.4 (covariances suffice). **Decision: Approach B results are not part of REPORT.md and are documented only here and in `PlanB/README.md`**
- [x] Figure scripts (2026-09-27; convention: each folder contains only scripts, each script writes a png with its own name next to it):
  - `results/reproduce_Jin/causal_graph/jin_causal_graph.py` (RLCD environment): reruns RLCD structure learning in memory (restarts itself with PYTHONHASHSEED=1, `--seed` to change), compares edge by edge with the authors' graph; edge labels are the authors' published coefficients (nll1.586); reproduces paper Fig. 1. Result: identical to the authors' graph. Writes nothing into the repository's results_* folders
  - `results/reproduce_Jin/truth_comparison/jin_L1_Rbirth_L2_Rg.py` (analysis environment): authors' coefficients + hierarchical least squares for L1, L2 + linear conversion; reproduces both panels of paper Fig. 3: Rb = −1.969 L1 + 9.309, r = 0.640; Rg = −1.289 L2 + 9.518, r = 0.976
  - `results/TNG/causal_graph/tng_causal_graph.py` (analysis environment): 12 galaxies, authors' graph + each galaxy's maximum-likelihood coefficients (smallest-NLL row of `<version>_<ID>_ci.csv`), observables with noise variance ≈ 0 marked; `--version noise` writes `_noise.png`
  - (checked 2026-09-27) the TNG L1 computation is correct: differs by < 1e-7 from a notebook-style per-star Nelder–Mead; the model-optimal estimate E[L1|X] gives equally low |r|. The L1 panels are narrow vertical bands because the correlation is weak (spread of f = |r| × spread of the truth); the maximum-likelihood solution makes L1 ≈ [Fe/H], [Mg/H] (|r| ≈ 0.9), which are almost unrelated to the birth radius in TNG50
  - (2026-09-27) wider selections were previewed with `scripts/local/preview_selection.py` (output `results/selection_preview.txt`): widening to 4<R<12 without a |z| cut only raises the median spread of the birth radius from 2.21 to 2.46 kpc (authors' g2.79: 3.67), the median |r| between [Fe/H] and birth radius stays ≈ 0.15–0.17 (authors: 0.70), and the linear ceiling barely changes (0.53–0.88). Conclusion: a moderately wider selection makes f(L1) slightly wider but does not change the fact that chemistry carries almost no birth-radius information in TNG50. The selection constants are lines 47–50 of `make_rlcd_input.py` (R_MIN, R_MAX, Z_MAX, FEH_MIN); change OUT_DIR when rerunning to keep the old results
  - `results/TNG/truth_comparison/tng_L1_Rbirth_L2_Rg.py` (analysis environment): same coefficients; f(L1)–Rbirth and f(L2)–Rg for the 12 galaxies, with the |r| range of indistinguishable solutions and the age/[Fe/H] baselines in each panel
- [x] **Selection R_0_30_Z_30 (2026-09-27)**: 0<R<30 kpc, |z|<30 kpc, other cuts unchanged. All files in `star_condition/R_0_30_Z_30/` (README explains the subfolders).
  - The data (6 GB) are compressed locally by `output_data/compute_moments.py` into an 8×8 covariance matrix per galaxy; `output_data/identifiability_moments.py` runs the same analysis on the covariances only. L1, L2 are linear combinations of the inputs, so their r follows from the covariance; checked on the old data, identical to 6 decimal places. The `truth_comparison` figure is made star by star
  - Results (nonoise): |r(L1,Rbirth)| mostly 0.66–0.86 (indistinguishable range); but J_φ alone gives 0.78–0.90 (across the whole disk the current radius ≈ birth radius). The most favourable L1 only slightly beats the best single input in 4 galaxies (530330, 538905, 546870, 557721); the ranges remain wide for 546870, 560421, 569251 (lower ends 0.03–0.13). The noise version has wider ranges. Linear ceiling 0.74–0.93
  - Interpretation: the apparent "improvement" of L1 comes mainly from the trivial information "where it is now ≈ where it was born" across the whole disk; the paper uses only 7–10 kpc precisely to exclude it
  - Note: the earlier preview only widened to 4<R<12, so its conclusion "little help" does not apply to the whole-disk range
- [x] **Selection R_0_30_Z_3 (2026-09-27)**: the stars of R_0_30_Z_30 with |z|<3 kpc (`input_data/make_input.py`, only ~3% fewer stars; the noise version reuses the same errors), same pipeline. Results essentially as for R_0_30_Z_30, but the optimum changed for 546870 (best |r(L1,Rb)| 0.71 → 0.09), 569251 (0.46 → 0.05) and 557721 (0.82 → 0.58) — further evidence that the maximum-likelihood solution is unstable
- [x] **Selection R_0_30_Z_0.3 (2026-09-27)**: the stars of R_0_30_Z_30 with |z|<0.3 kpc (0.188–0.313 million per galaxy, 2.79 million in total, about half of |z|<3), same pipeline, all files in `star_condition/R_0_30_Z_0.3/` (README has the key results). nonoise maximum-likelihood |r(L1,Rb)|: 0.79–0.89 in 9 galaxies, 0.45 for 552581, 0.28 for 569251, 0.17 for 546870; no galaxy's optimum beats "J_φ alone" (0.71–0.92). Compared with |z|<3 the optimum changed again for 557721 (0.58→0.89) and 569251 (0.05→0.28); in the noise version the optimum of 545437/547293/557721/560421 drops to 0.10–0.41. Many indistinguishable solutions remain (Jacobian rank 14/15). Tip: set OMP_NUM_THREADS=1 when two processes run at once (otherwise threads compete on a 2-core machine and it is more than 2× slower)
- [x] **Approach B (2026-09-27, main selection 7<R<10, |z|<0.3, as in the paper)**: all files in `PlanB/` (input_data, true_value, output_data, causal_graph, truth_comparison; README has the run order and results). Scripts: `PlanB/input_data/make_input.py` (copies from data/rlcd_input) → `PlanB/output_data/rlcd_B.py --seeds 30` (30 PYTHONHASHSEED values per galaxy, RLCD environment, ~2 s each) → `rlcd_B_eval.py` (for structures with frequency ≥ 20%: Jacobian rank, 50-start maximum likelihood, latent weights by hierarchical least squares, also `eval_<version>_best.json`) → `summary_planB.py` (comparison with Approach A, incl. dBIC) → `PlanB/causal_graph/planB_causal_graph.py`, `PlanB/truth_comparison/planB_latent_vs_truth.py`. `rlcd_B_eval` was checked with the authors' graph (538905: NLL 0.4918, rank 14/15, consistent with Approach A). Note: `rlcd_B_eval.py` must import `rlcd_B` from its own folder (otherwise the older copy in scripts/local, which writes to old paths, would be used).
  Results: identical to the authors' graph 0/360 (nonoise and noise); in 8 galaxies all 30 seeds give the same graph; 0–4 latents, 8–17 directed edges; 545437 and 569251 have no latent and undirected edges (no coefficients); 4 galaxies violate the age prior; of the 10 fittable galaxies only 3 (546870, 562029, 567607, each with 1 latent) are identifiable. Picking with the truth the latent most similar to the birth radius (best case): median |r| 0.51 (0.24–0.77), only 530330, 552581, 562029 beat age alone; noise median 0.27. After a BIC penalty for the number of edges, the authors' graph fits better than the graph RLCD learns in every comparable galaxy (dBIC per star +0.003 to +0.42; noise: all except 487742). Conclusion: RLCD does not learn an "L1 = birth radius, L2 = guiding radius" graph on TNG50 by itself.
- [ ] Optional: Approach B with the wide selections (`rlcd_B.py --cond R_0_30_Z_0.3` etc., using covariance-matched surrogate data; output in `PlanB/output_data/<selection>/`; the plotting scripts currently read only the main selection) — not decided

### Stage 3 (later)
- [ ] If the TNG50 validation succeeds, apply to Milky Way data; decide whether to ask the authors for their data or redo it from scratch. Currently on hold (the TNG50 validation did not succeed).

---

## 3. Research goal and notes on the paper

### 3.1 Goal
Test the RLCD (Rank-based Latent Causal Discovery) method of Jin et al. with the TNG50 simulation:
- **6 observable inputs**: Age, [Fe/H], [Mg/H], J_z, J_φ, e
- **program output**: two latent variables L1, L2
- **comparison with the TNG50 truth**: L1 ↔ birth radius Rbirth, L2 ↔ guiding radius Rg

### 3.2 The paper (Jin, Lu, Ting, Zheng, Buck; [arXiv:2507.00134](https://arxiv.org/abs/2507.00134))
- Simulations: NIHAO-UHD g2.79e12 (main) and g8.26e11 (transferability). Disk stars with [Fe/H] > −1, currently at 7–10 kpc. Errors added: 10% in age, 0.02 dex in [Fe/H], 0.06 dex in [Mg/Fe]. Low-α disk defined by [Mg/Fe] < −0.01; only the low-α disk is studied.
- Orbits: AGAMA, axisymmetric multipole potential (dark matter, gas, stars).
- Observations: APOKASC-3. Actions with galpy Stäckel fudge + gala MilkyWayPotential2022, v_c(R⊙) = 229 km/s.
- Causal graph found: L1 is the root and points to [Fe/H], [Mg/H], Age, J_z and L2; L2 points to J_φ and e. L1 corresponds to Rbirth (Pearson r = 0.64 on g2.79), L2 to Rg (very tight).
- Conversion: Rb = −1.97 L1 + 9.31, Rg = −1.29 L2 + 9.52 (g2.79).
- Prior: Age has no parent among the observables.

### 3.3 Actual data format of the repository https://github.com/ZehaoJin/Causal-Galactic-Archaeology
- RLCD reads `data/g2.79e12_e.csv`, 6 columns: `logJZ, JPhi, ecc, age, FeH, MgH`
  - `logJZ` = log₁₀(J_z)
  - `JPhi` in kpc·km/s, positive in the simulation (negative in the observational file)
  - `MgH` = MgFe + FeH
- Truth in `data/g2.79e12_multi_ele_action_Jin_SN.csv`, columns: `age, iord, JR, JZ, JPhi, FeH, MgFe, ecc, Rg, R_bir, z_bir, R`
- **Differences between the data and the paper text**:
  - actual radial range R = 9.4–12.1 kpc (not 7–10 kpc)
  - actual [Fe/H] > −0.5 (not −1)
  - R_bir at most ~20 kpc, apparently an additional cut R_bir < 20
  - no visible observational errors in the csv
  - Rg ≈ JPhi / 347 (g2.79), i.e. Rg satisfies L_z = Rg · v_c(Rg)
- L1/L2 are computed in `L1L2calculation/*.ipynb` with **hard-coded coefficients**, which must be replaced with newly learned ones for new data. The notebook uses the natural log `np.log(JZ)`, equivalent to log₁₀ after standardization.
- Commands (in `rlcd/`):
  - structure learning: `python run_rlcd.py --s <scenario> --alpha0 1e-5 --alpha1 1e-15 --alpha2 1e-5 --alpha3 1e-15`
  - coefficient learning: `cd ParameterLearning; python linear_parameter_identification.py --s <scenario> --dot_path <.dot> --lr 0.005 --epochs 600 --lr_decay plateau`
  - new data need a scenario function in `rlcd/scenarios.py` (read the csv, prefix columns with `X_`)
- Dependencies: Python 3.12+, Java (Tetrad), Graphviz, `rlcd/py312_requirements.txt`

---

## 4. Repository layout

```
<repo root>/
├── README.md, STATUS.md, REPORT.md, DATA.md, data_manifest.tsv, CLAUDE.md, AGENTS.md, LICENSE, env.example
├── docs/REPRODUCTION.md                   ← this file
├── docs/INSTALL_AGAMA.md
├── tools/                                 helper scripts (unpack/link data, check setup, clone the RLCD repository)
├── setup/                                 requirements of the two Python environments
├── scripts/
│   ├── server/                            run on TNG JupyterLab
│   │   ├── data_obtain.py                 particles + birth radii of the 12 galaxies → ./data_obtain/<ID>.hdf5
│   │   └── merge_hdf5.py                  merged into data_obtain_all.hdf5
│   ├── local/                             run locally
│   │   ├── check_vs_001_orbit.py          cross-check with the earlier 001_orbit results (needs non-public legacy files)
│   │   ├── compute_orbits.py              AGAMA: potential, actions, Rg, eccentricity
│   │   ├── check_orbits.py                checks of the orbit results → results/orbit_check/
│   │   ├── make_rlcd_input.py             selection + RLCD inputs → data/rlcd_input/
│   │   ├── preview_selection.py           preview of wider selections → results/selection_preview.txt
│   │   ├── compare_rlcd_graphs.py         compare causal graphs (latent relabeling + prior check)
│   │   ├── rlcd_A.py                      Approach A: frozen authors' graph → parameter learning → L1/L2 → truth (RLCD env) ⚠ not identifiable
│   │   ├── rlcd_identifiability.py        identifiability: Jacobian rank, multistart maximum likelihood, L1->age profile (numpy/scipy)
│   │   ├── rlcd_identifiability_summary.py summary of round 1 (loose criterion ΔNLL<0.01)
│   │   ├── rlcd_identifiability_ci.py     round 2, strict: polish + noise variances ≥ 0 + 95% likelihood ratio (main result, --extra 200)
│   │   ├── linear_ceiling.py              linear ceiling and single-variable baselines (supervised, reference only)
│   │   ├── rlcd_B.py, rlcd_B_eval.py      older copies of the Approach-B scripts (use PlanB/output_data/)
│   │   └── rlcd_requirements_mac.txt      = setup/requirements-rlcd.txt
│   └── reference/                         earlier scripts (read-only): 001_orbit.py, 001_orbit_new.py, TNG50_disk_trajectory.py, obtain_data_s98.py
├── Causal-Galactic-Archaeology/           the RLCD repository (not committed; tools/setup_rlcd_repo.sh), run inside its rlcd/
├── data/
│   ├── release/                           public main-selection data (tar.gz)
│   ├── tng50_extract/data_obtain_all.hdf5 raw server extraction (not public; do not modify)
│   ├── tng50_derived/
│   │   ├── orbits_all.hdf5                output of compute_orbits.py (not public)
│   │   ├── orbits_test.hdf5               output of --test (old potential setting, can be deleted)
│   │   └── potentials/<ID>_halo.coef, <ID>_disk.coef   potential coefficients
│   └── rlcd_input/{nonoise,noise}/<ID>_e.csv, <ID>_truth.csv   RLCD inputs and truth (unpacked from data/release/)
├── results/
│   ├── check_vs_001_orbit.txt             cross-check report
│   ├── compute_orbits.log                 orbit-computation log
│   ├── rlcd_g279/                         RLCD reproduction (authors' graph, our runs, comparison reports)
│   ├── orbit_check/                       orbit checks (summary + figures)
│   ├── rlcd_input_selection.txt           stars left after each selection step
│   ├── rlcd_A/                            Approach A (g279 validation, nonoise, noise, identifiability/)
│   ├── reproduce_Jin/{causal_graph,truth_comparison}/   reproduction of paper Figs. 1, 3
│   └── TNG/{causal_graph,truth_comparison}/             TNG50 equivalents
├── star_condition/                        wider selections R_0_30_Z_30, R_0_30_Z_3, R_0_30_Z_0.3 (+ R_7_10 mirror of results/TNG)
└── PlanB/                                 Approach B
```

**Related locations outside this repository (not public):**
- the 423-galaxy sample of an earlier project: `TNG50_disk_trajectory.hdf5` (375 MB; centre, velocity, angular momentum of every galaxy at every snapshot), `001_orbit/<ID>.txt` (earlier birth-radius results) and `001_orbit_s99_3D_z/<ID>.txt` (2026-04 version, with z and 3D distance columns)
- on the TNG server: a working directory for data_obtain.py (e.g. `~/MW/casuality/`), the original trajectory file (e.g. `~/423_sample/`), simulation data in `/home/tnguser/sims.TNG/TNG50-1/output`

---

## 5. Computing environment

### 5.1 Tested local setup (macOS, Apple Silicon arm64)
- Anaconda base was the terminal default; its llvm-openmp (14.0.6) is too old to compile AGAMA.
- **Analysis environment**: a venv created from Anaconda's Python 3.11.7, containing:
  - AnastrisTNG 2.0.0 (editable install; not needed by the local scripts)
  - pynbody 1.6.0, numpy 1.26.4, scipy 1.13.1, pandas 2.2.2, h5py 3.11, matplotlib 3.9, astropy 6.1, wheel 0.48
  - **agama 1.0**, built from source with Homebrew's libomp (`/opt/homebrew/opt/libomp`)
  - no jupyter
- Homebrew and the Xcode command-line tools are installed.

### 5.2 TNG JupyterLab server
- Simulation data: `/home/tnguser/sims.TNG/TNG50-1/output`
- Old AnastrisTNG (~1.x, API different from 2.0.0) and galpos
- data_obtain.py is run from a working directory into which the trajectory file is symlinked: `ln -s <path>/TNG50_disk_trajectory.hdf5 .`

### 5.3 RLCD (stage 2)
- **A separate Python 3.12 environment, used only for RLCD.** Reason: the repository requires Python 3.12 and pins many versions (scipy 1.11.4, pandas 2.3.1, matplotlib 3.10.5, torch 2.8.0, ...) that conflict with the analysis environment; keeping them apart leaves the AGAMA environment untouched.
- `setup/requirements-rlcd.txt` = the repository's `rlcd/py312_requirements.txt` without nvidia-* and triton (Linux + GPU only) and without dataclasses and typing (backports for old Pythons).
- Graphviz (`brew install graphviz`) draws the graphs.
- Java is only needed for `--stage1_method fges` (Tetrad); the default `ges` uses causal-learn and needs no Java.
- The repository is cloned into `Causal-Galactic-Archaeology/` (~260 MB, including precomputed results and the Tetrad jar). RLCD commands run in its `rlcd/` folder and write to `rlcd/results_<scenario>/`.

---

## 6. Pipeline details

### 6.1 Server-side extraction: `scripts/server/data_obtain.py`
**Basis**: `scripts/reference/001_orbit_new.py` of the earlier project, keeping all its fixes.

**Sample (12 galaxies, Subhalo IDs at snap 99)**:
487742, 530330, 538905, 545437, 546870, 547293, 552581, 557721, 560421, 562029, 567607, 569251
(607654 from the original list is not in the 423-galaxy sample and has no trajectory; removed.)

**For each galaxy**:
1. `snap99.load_particle(ID)` loads dm/star/gas/bh. For stars additionally `GFM_StellarFormationTime, BirthPos, BirthVel, GFM_InitialMass, GFM_Metals`. Then `physical_units()`.
2. The stellar arrays BirthPos/BirthVel etc. **must be extracted before face_on**, because face_on rotates all 3D arrays in place.
3. Wind particles are removed; only stars with `0 < aform ≤ 1`, `tform > 0`, both finite, are kept.
4. Abundances from number ratios, solar values from Asplund et al. 2009:
   - [Fe/H] = log₁₀((X_Fe / A_Fe) / (X_H / A_H)) + 4.50
     Meaning: the ratio of the numbers of Fe and H atoms, relative to the same ratio in the Sun, on a log scale.
     - X_Fe, X_H: mass fractions of Fe and H (columns 8 and 0 of `GFM_Metals`)
     - A_Fe = 55.845, A_H = 1.008: atomic masses, converting mass ratios to number ratios
     - 4.50: solar 12 − logε_Fe = 12 − 7.50
   - [Mg/Fe] = log₁₀((X_Mg / A_Mg) / (X_Fe / A_Fe)) − 0.10
     - X_Mg: mass fraction of Mg (column 6 of `GFM_Metals`)
     - A_Mg = 24.305: atomic mass of Mg
     - 0.10: solar logε_Mg − logε_Fe = 7.60 − 7.50
5. Formation-time clamp: stars formed later than the end of the trajectory by less than 1e-3 Gyr (born near z=0, floating-point rounding) are clamped to the end of the trajectory.
6. **Birth radius**: galpos' `GalaxyPoseTrajectory` reads `t, cpos, cvel, ang_mom` from the trajectory file, `box_size = 35000/0.6774`. `pynbody_decorate.make_star_birth(...)` and `align_with_galaxy(orientation_align=True)` transform each star's BirthPos into "the disk frame of the main progenitor at the moment of its birth".
   - `birth_R_2D` = √(x_b² + y_b²): distance of the birth position from the centre within the disk plane at that time (cylindrical radius); x_b, y_b are the in-plane coordinates in the disk frame at birth. **This is the Rbirth truth.**
   - `birth_R_3D` = √(x_b² + y_b² + z_b²): straight-line distance of the birth position from the galaxy centre at that time; z_b is the height above the disk plane at birth.
   - Stars formed before the start of the trajectory get NaN (`has_birth = False`). They are **kept**, because their mass is needed for the potential.
7. `face_on(rmax=8, alignwith='star')`: rotate the whole subhalo into the z=0 disk frame, then save all particles.
8. One HDF5 per galaxy (float32, gzip); written to a temporary file and renamed, so a restarted run skips finished galaxies.

**Commands on the server**:
```bash
cd <working directory>
ln -s <path>/TNG50_disk_trajectory.hdf5 .        # once, before the first run
nohup python3 data_obtain.py > data_obtain.log 2>&1 &
python3 merge_hdf5.py                             # → data_obtain_all.hdf5
```
Download it from JupyterLab into `data/tng50_extract/`.
(The `script` attribute in the header of the downloaded file reads `001_orbit_causal.py`, the file's name before it was renamed; same content.)

**Structure of `data_obtain_all.hdf5`**
```
/<ID>/header   (attributes only) ID, snap=99, a≈1, h=0.6774, time_Gyr=13.80, face_on_rmax_kpc=8,
               face_on_alignwith='star', units, solar_scale, script
/<ID>/frame    center_pos(3), center_vel(3), rotation(3,3)     — return values of the z=0 face_on
/<ID>/stars    iord(u64), pos(N,3) kpc, vel(N,3) km/s, mass Msun, init_mass Msun, age Gyr, aform,
               FeH, MgFe, X_H, X_Mg, X_Fe, birth_pos(N,3) kpc [disk frame at birth],
               birth_R_2D kpc, birth_R_3D kpc, has_birth(int8, 0/1)
/<ID>/wind     pos, mass       (wind particles, only for the potential)
/<ID>/dm       pos, mass
/<ID>/gas      pos, mass
/<ID>/bh       pos, mass       (one black hole per galaxy)
```
pos and vel are in the z=0 disk frame: origin at the galaxy centre (ssc method), z axis along the angular momentum of the stars within 8 kpc.

**Particle numbers (from the cross-check report)**
| ID | stars | wind | dark matter | gas | no birth radius | median v_φ of disk stars (km/s) |
|---|---|---|---|---|---|---|
| 487742 | 608,237 | 15 | 3,011,552 | 1,869,429 | 31 | 158.8 |
| 530330 | 1,053,023 | 0 | 2,920,963 | 275,220 | 3 | 191.4 |
| 538905 | 893,618 | 0 | 2,294,758 | 507,938 | 1 | 181.7 |
| 545437 | 722,225 | 0 | 1,859,730 | 645,766 | 66 | 171.3 |
| 546870 | 739,723 | 0 | 2,323,462 | 171,851 | 24 | 158.9 |
| 547293 | 908,363 | 0 | 1,755,535 | 605,028 | 28 | 183.2 |
| 552581 | 779,508 | 1 | 1,771,153 | 419,905 | 55 | 182.4 |
| 557721 | 712,029 | 0 | 1,947,701 | 136,165 | 190 | 187.4 |
| 560421 | 526,387 | 2 | 1,296,156 | 455,951 | 7 | 166.0 |
| 562029 | 541,225 | 1 | 1,622,644 | 411,809 | 0 | 163.6 |
| 567607 | 515,043 | 27 | 1,494,362 | 388,236 | 156 | 173.2 |
| 569251 | 670,495 | 0 | 1,457,537 | 392,465 | 37 | 181.5 |

The median TNG50 [Mg/Fe] is about **+0.28 to +0.33**, much higher than in NIHAO, so the paper's low-α cut cannot be copied. Median [Fe/H] is about −0.21 to −0.27.

**Cross-check results (1.3)**
- Ages, masses and current positions agree star by star with `001_orbit` and `001_orbit_s99_3D_z` (differences < 1e-5).
- **`Current_R` and `Birth_R` in the old `001_orbit/<ID>.txt` are both 3D distances** (not cylindrical radii), consistent with `birth_R_3D`.
- Birth radii differ from `001_orbit_s99_3D_z` by ≤ 5 pc (about 3% of the stars), negligible.
- The new data give birth radii for a few stars born near z=0 that were NaN in the old data; this is the expected effect of the formation-time clamp.

### 6.2 Local orbit computation: `scripts/local/compute_orbits.py`
**Run**:
```bash
cd <repo root>
$PY scripts/local/compute_orbits.py --test      # test: 5000 random stars of the first galaxy → orbits_test.hdf5
caffeinate -i $PY scripts/local/compute_orbits.py 2>&1 | tee results/compute_orbits.log   # everything (macOS)
$PY scripts/local/compute_orbits.py --ids 487742 530330   # selected galaxies
```
- `caffeinate -i` (macOS): keep the machine awake during the run
- `2>&1 | tee ...`: show output and errors on screen and write them to the log
- Galaxies already computed with the same version are skipped; those with another version are deleted and recomputed. After an interruption simply rerun.

**Units**: `agama.setUnits(mass=1, length=1, velocity=1)`, i.e. Msun, kpc, km/s; time unit 0.978 Gyr.

**Potential** (`build_potential`):
- disk part: **CylSpline** (axisymmetric cylindrical grid, 30×30, R and z from rmin to 60 kpc, first z node 0.05 kpc): stars with r < 50 kpc plus gas with R < 30 kpc and |z| < 5 kpc
- halo part: **Multipole** (axisymmetric, lmax = 8, 40 radial nodes, rmin to 500 kpc): all dark matter, black holes, wind particles, and the remaining stars and gas
- total potential = sum of the two
- **inner grid boundary rmin = 0.2 kpc**, about the TNG50 softening at z=0. Initially AGAMA chose the grid automatically and went down to 0.01 kpc, where particle noise made the circular angular momentum non-monotonic in radius and the ActionFinder of galaxy 530330 failed. If 0.2 kpc still fails, 0.5 and 1.0 kpc are tried automatically; the value used is stored in the attribute `pot_rmin_kpc`.
- the coefficients of both parts are exported to `data/tng50_derived/potentials/<ID>_{halo,disk}.coef` and can be reloaded with `agama.Potential(file=...)`

**Quantities per star** (from `/<ID>/stars/pos` and `vel` only):
- **E** = Φ(x) + ½|v|²
  Meaning: total energy = potential energy at the position + kinetic energy (per unit mass).
  - Φ(x): potential at the star's position; |v|: speed

  Only bound stars (E < 0) are processed; all orbit quantities of the others are NaN (about 0.6%, mostly distant halo stars).
- **J_φ** = x·v_y − y·v_x
  Meaning: angular momentum about the z axis (L_z); positive means co-rotating with the disk.
  - x, y: in-plane coordinates; v_x, v_y: corresponding velocity components

  Computed for all stars.
- **J_R, J_z**: `agama.ActionFinder(pot, interp=False)`, Stäckel fudge, returns (J_R, J_z, J_φ)
- **Rg** = `pot.Rcirc(L=|J_φ|)`, the radius where |J_φ| = Rg · v_c(Rg).
  Meaning: the radius of a circular orbit with this angular momentum.
  - v_c(R) = √(R · ∂Φ/∂R): circular speed at radius R
  - ∂Φ/∂R: radial derivative of the potential, i.e. the radial force in the plane
- **r_peri, r_apo, e, z_max**: orbits integrated with `agama.orbit` for 15 × `pot.Tcirc(xv)` (period of the circular orbit with the star's energy), 1501 equally spaced points, accuracy 1e-7, 5000 orbits per batch
  - r_peri = minimum of the 3D distance r along the orbit, r_apo = maximum
  - **e = (r_apo − r_peri) / (r_apo + r_peri)**
    Meaning: the flatter the orbit, the closer e is to 1; nearly circular orbits have e close to 0. r_apo: apocentre distance, r_peri: pericentre distance.
  - z_max = maximum |z| along the orbit
  - **The paper does not say whether the eccentricity uses 3D r or cylindrical R**; 3D r is used here (to be confirmed)
- rotation curve: 60 points at R = 0.5–30 kpc, v_c = √(−R · F_R), F_R the radial force

**Structure of `orbits_all.hdf5`**:
```
/<ID>/Jr, Jz, Jphi, E, Rg, rperi, rapo, ecc, zmax    float32, one entry per star,
                                                     rows in exactly the order of /<ID>/stars in data_obtain_all.hdf5
/<ID>/rotation_curve/R, vc
/<ID> attributes: n_stars, test_mode, pot_version='v2-rmin-grid', pot_rmin_kpc, N_PERIOD=15, TRAJ_SIZE=1501, units
(test mode also stores star_index, the indices of the sampled stars)
```

**Test results (487742, 5000 stars, old automatic grid; full results with the new setting are almost identical)**:
- 4963 of 5000 stars bound
- all checks passed: J_φ equals the directly computed L_z; every current r lies within [r_peri, r_apo]; |z| never exceeds z_max; median Rg · v_c(Rg) / |L_z| = 1.000
- v_c(8 kpc) = 185 km/s, 171–188 km/s within 0.5–28 kpc, fairly flat
- disk medians (4 < R < 12, |z| < 1): J_z = 18, J_φ = 1074, e = 0.25, Rg − R = −0.7 kpc
- young stars (< 3 Gyr) have median e ≈ 0.20, "hotter" than the Milky Way thin disk (≈ 0.1); a property of the TNG50 disks, not an error

**Full run**: 487742 took 4.4 minutes (v_c(8) = 185.4 km/s with the new setting, consistent with the test); 530330, with ~1.05 million stars, about 8–10 minutes.

### 6.3 RLCD reproduction (authors' g2.79e12 data)
Commands in Sec. 8.3. Results are written to `Causal-Galactic-Archaeology/rlcd/results_causal_vertical_heating_g279_e/result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot` (overwriting the authors' committed version; their original is available with `git show HEAD:<path>`).
- **Authors' graph (15 edges)**: L1 → {L2, age, FeH, MgH, logJZ}; L2 → {JPhi, ecc}; MgH → {FeH, JPhi, logJZ, ecc}; age → {JPhi, FeH, logJZ, ecc}. This matches the coefficients hard-coded in `L1L2calculation/GetL1_L2_279sim.ipynb`.
- **Our first graph (18 edges)**: the two latent names are swapped (our L2 is the authors' L1). Same latent structure: root latent → {other latent, age, FeH, MgH, logJZ}, other latent → {JPhi, ecc}. All 15 of the authors' edges plus 3 more: **MgH → age**, age → guiding-radius latent, MgH → guiding-radius latent.
- MgH → age violates the prior of paper Sec. 3.2 (age has no observed parent). The paper grid-searches α to get several candidate graphs and filters them with this prior down to a unique graph.
- **But**: the file name of the authors' graph contains exactly this α combination (1e-5/1e-15/1e-5/1e-15), i.e. with the same parameters our graph differs from theirs. Possible reasons considered at that point: different data or code versions, different package versions, or manual editing according to the prior. (The actual reason turned out to be the randomness, see below.)
- Only 3 more α combinations were run instead of all 81. All files in `results/rlcd_g279/`.
- Comparison with `scripts/local/compare_rlcd_graphs.py` (aligns latent labels by connectivity and checks the prior; the reference can be `git:<repo dir>:<relative path>` to read the authors' original at git HEAD). Report: `results/rlcd_g279/compare_4alpha.txt`.

| α (0/1/2/3) | latents | authors' 15 edges | extra edges | prior (age has no observed parent) |
|---|---|---|---|---|
| 1e-5 / 1e-15 / 1e-5 / 1e-15 (as in the authors' file name) | 2 | 15/15 | MgH→age, MgH→L2, age→L2 | fail |
| 1e-5 / 1e-10 / 1e-5 / 1e-10 | 2 | 15/15 | age→MgH, MgH→L2, age→L2 | pass |
| 1e-10 / 1e-15 / 1e-10 / 1e-15 | 4 | 2/15 | completely different | pass |
| 1e-15 / 1e-15 / 1e-15 / 1e-15 | 4 | 3/15 | completely different | pass |

  - α0 = 1e-5: latent structure as in the paper (L1 → L2, age, FeH, MgH, logJZ; L2 → JPhi, ecc), all authors' edges present, only 3 extra: one between age and MgH (direction depends on α), and one each from age and MgH to L2.
  - α0 = 1e-10 or 1e-15: 4 latents in two pairs — one pair pointing to {age, FeH, MgH}, the other to {JPhi, ecc, logJZ} — different from the paper but also satisfying the prior.
  - So "filtering by the prior" does not lead to a unique graph as the paper states; the authors' final graph has fewer edges than any direct output.
- Stage-1 GES puts all 6 variables in one group (`result_stage1ges_sparsit2.dot`).
- Note: `grid_search.py` in the repository calls `python`; activate the RLCD environment first, otherwise the wrong interpreter is used.

**Files in `results/rlcd_g279/`**
- `author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot/.png`: the authors' original graph (from git HEAD of their repository); `author_rerendered_same_style.png` is the same graph redrawn with the local Graphviz
- `ours_seed{0..4}_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot/.png` + `ours_seed{0..4}.log`: 5 runs with the same α and different PYTHONHASHSEED
- `ours_noseed_result_alpha<four α>_ranktestscale1.dot/.png` + `ours_noseed_alpha<...>.log`: 4 unseeded runs (4 α combinations). The 1e-5/1e-15/1e-5/1e-15 file in the repository was overwritten by later runs; its dot was rewritten from the output at the time and the png redrawn
- `ours_stage1ges_sparsit2.dot/.png`: result of stage-1 GES (last run)
- `compare_4alpha.txt`, `compare_5seeds.txt`: comparison reports; `compare_author_vs_seed1.png`: authors' graph and the seed-1 graph side by side
- `compare_author_vs_seed3_diff.png` (and `diff_author_vs_seed3_*.png`): authors' graph vs seed 3 with the differing edges coloured. The two graphs share 11 edges; the 4 differing edges are authors' L1→FeH, L1→logJZ, L2→JPhi, L2→ecc vs seed-3 L1→JPhi, L1→ecc, L2→FeH, L2→logJZ: (JPhi, ecc) and (FeH, logJZ) are attached to opposite latents, so the second latent no longer means "guiding radius". The graphs look alike at first glance; check every edge
- The two overwritten files in the authors' `rlcd/results_causal_vertical_heating_g279_e/` were restored with `git checkout`; our new files were moved here. **When running RLCD inside the repository, results are written back into its results_* folders: move them out promptly to avoid overwriting the authors' files.**

**Randomness test (same α 1e-5/1e-15/1e-5/1e-15, 5 repeats)**
- Command: in `rlcd/`, `PYTHONHASHSEED=$s $PY_RLCD run_rlcd.py ...`, s = 0–4, copying the dot after each run (now `results/rlcd_g279/ours_seed$s_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot`). Report: `results/rlcd_g279/compare_5seeds.txt`.
- Results:
  - seed 1: **identical to the authors' graph** (15 edges)
  - seeds 0, 2, 3: another 15-edge structure — root latent → {other latent, JPhi, ecc, MgH, age}, other latent → {FeH, logJZ}. Satisfies the prior but has a different physical meaning
  - seed 4: the seed-0 structure plus 3 edges (X_age→X_MgH, and X_MgH, X_age → non-root latent). Satisfies the prior
  - the very first unseeded run: the authors' structure plus 3 edges, with MgH→age violating the prior
- **Cause** (from the source code):
  - (1) RLCD stores variable names in Python sets and performs tests and latent merges in set order; Python randomizes string hashing at every start-up, so the iteration order differs between runs.
  - (2) `StructureLearning/RLCD/LatentGroups.py` line 815, `findRandomLatentChild`, calls `random.sample(...)` without a seed (used only for "non-atomic covers"). `bootStrapCovariance` in `misc.py` line 211 also uses random numbers but is never called.
- **On the prior of paper Sec. 3.2**: the prior is only a filter — it removes candidate graphs where age has an observed parent, it does not modify graphs. It removes the 18-edge graph of the first run, but the graphs of seeds 0, 2, 3, 4 all pass with structures different from the authors'. "Filtering by the prior converges to a unique graph" was therefore a coincidence of the authors' runs; their graph is one possible RLCD output (1 of 5 here).
- **Frequencies over 30 seeds (2026-09-27, PYTHONHASHSEED = 0–29, same α)**: authors' graph 7 times (1, 7, 8, 11, 14, 22, 25; 23%); the seed-0 15-edge graph is most common (12 times, 40%); the other 11 runs give four different 18-edge graphs.
- **Every future RLCD run must** fix PYTHONHASHSEED and the random seeds and be repeated with several seeds, reporting structure frequencies.
- **Do the authors mention the randomness?** (checked 2026-09-26)
  - paper: no mention of randomness, seeds or repeated runs; only that the prior filter "converges to a unique graph".
  - Jin et al.'s repository: not in the README or notebooks. The code has `utils/utils.py: set_random_seed()` (sets random, numpy, torch seeds and PYTHONHASHSEED), but it is **never called**; the comment of `reorderCovers` in `misc.py` says "Reduces some randomness", so the original RLCD authors knew about order-induced randomness and handled it only partly.
  - the original RLCD repository (github.com/dongxinshuai/scm-identify, Dong et al.): the README only warns that **parameter learning** "may need multiple run to get ideal results due to random initialization"; nothing about structure learning.
  - so the 13 `*_param_by_trek_nll*.dot` files in Jin et al.'s repository are repeated parameter-learning runs for the same graph; the number in each file name is the NLL (smaller = better fit), the smallest being nll1.4979. **Parameter learning must also be repeated many times.**

### 6.4 Orbit checks (1.5): `scripts/local/check_orbits.py`
Run: `$PY scripts/local/check_orbits.py` (uses h5py, or the pure-Python pyfive reader if h5py is missing). Output in `results/orbit_check/`: `summary.txt`, `summary.csv`, `rotation_curves.png`, `ecc_disk.png`, `Rg_vs_R.png`, `logJz_disk.png`; plus `MgFe_FeH.png`, the [Mg/Fe]–[Fe/H] distribution for 4<R<12, |z|<3 (made for the selection discussion).
- All 12 galaxies pass: bound fraction 99.4–100% (lowest 99.42% for 487742); no NaN among bound stars; |J_φ − L_z| < 0.004; all r within [r_peri, r_apo]; |z| ≤ z_max; median Rg·v_c(Rg)/|L_z| = 1.000; e within [0, 1]. All galaxies use rmin = 0.2 kpc (no retries).
- v_c(8 kpc) = 185–229 km/s, maximum v_c 191–251 km/s.
- Disk medians (4<R<12, |z|<1): J_z 9.9–18.5, J_φ 1074–1403, J_R 27–99, e 0.18–0.32 (young stars < 3 Gyr: 0.09–0.25), Rg − R from −0.36 to −0.85 kpc, counter-rotating 0.4–6.5%.
- The horizontal band at Rg ≈ 0 in the Rg–R figure is non-rotating bulge/halo stars (L_z ≈ 0), as expected.
- [Mg/Fe]: most galaxies have **no clear high-/low-α sequences**; only 545437, 567607, 560421, 487742 show a clump of low-α stars around [Mg/Fe] ≈ 0.25.

### 6.5 Selection and RLCD inputs (1.6): `scripts/local/make_rlcd_input.py`
Run: `$PY scripts/local/make_rlcd_input.py` (~2 minutes). Statistics: `results/rlcd_input_selection.txt`.
**Selection (fixed 2026-09-26)**:
1. bound (E < 0), with J_R, J_z, J_φ, Rg and eccentricity defined
2. has_birth, and birth radius R_bir < 20 kpc
3. **7 < R < 10 kpc** (paper text)
4. **|z| < 0.3 kpc**
5. **[Fe/H] > −0.5, no [Mg/Fe] cut** (TNG50 has no clear α bimodality)
6. J_z > 0 (logarithm); counter-rotating stars kept (as in the authors' data)

**Observational errors (noise version)**: age multiplied by (1 + 0.10·N), [Fe/H] + 0.02·N, [Mg/Fe] + 0.06·N (N a standard normal random number), [Mg/H] from the noisy [Mg/Fe] + [Fe/H]. Random seed = 20260926 + galaxy ID. The selection uses error-free values, so both versions contain the same stars row by row.

**Outputs** (two files per galaxy):
- `data/rlcd_input/nonoise/<ID>_e.csv`, `noise/<ID>_e.csv`: 6 columns `logJZ, JPhi, ecc, age, FeH, MgH`, same format as the authors' `data/g2.79e12_e.csv`
- `<ID>_truth.csv` (row-matched to `_e.csv`): `iord, star_index` (row in data_obtain_all.hdf5 /stars)`, R, z, Rg, R_bir, R_bir_3D, z_bir, age, FeH, MgFe, MgH, JR, JZ, JPhi, ecc, rperi, rapo, zmax, E`; the noise version also has `age_obs, FeH_obs, MgFe_obs, MgH_obs`

**Final sample sizes**: 487742: 13,190; 530330: 24,154; 538905: 34,970; 545437: 16,222; 546870: 11,484; 547293: 32,525; 552581: 13,573; 557721: 25,326; 560421: 19,696; 562029: 7,299; 567607: 12,496; 569251: 29,035. Total 239,970. The R range and |z| < 0.3 are the cuts that remove most stars.
- Medians: logJZ 0.48–0.93, J_φ 1356–1753, e 0.13–0.28, age 1.3–7.6 Gyr, R_bir 5.3–9.0 kpc (5–95%: about 2.5–13.5 kpc), Rg 7.1–8.3 kpc.
- Checked: `10**logJZ` equals JZ in the truth; JPhi consistent; no NaN; in the noise version the relative age errors have standard deviation 0.100 and the [Fe/H] errors 0.020.

### 6.6 Approach A: `scripts/local/rlcd_A.py` (RLCD environment)
```bash
cd <repo root>
$PY_RLCD scripts/local/rlcd_A.py --version g279 2>&1 | tee results/rlcd_A_g279.log
$PY_RLCD scripts/local/rlcd_A.py --version nonoise --ids 538905 2>&1 | tee results/rlcd_A_test538905.log
$PY_RLCD scripts/local/rlcd_A.py --version nonoise 2>&1 | tee results/rlcd_A_nonoise.log
$PY_RLCD scripts/local/rlcd_A.py --version noise 2>&1 | tee results/rlcd_A_noise.log
```
Options: `--n-repeat` (parameter-learning repeats), `--epochs 600 --lr 0.005 --lr-decay plateau` (as in the repository README), `--dot` defaults to the authors' graph `results/rlcd_g279/author_...dot`.
1. **Data**: `data/rlcd_input/<version>/<ID>_e.csv` (g279 uses the repository's `data/g2.79e12_e.csv` with truth `g2.79e12_multi_ele_action_Jin_SN.csv`, row-matched, checked). Standardized first: subtract the mean and divide by the standard deviation of each column (as `DataModel.generate_data(normalized=True)` in the repository), columns prefixed with `X_`.
2. **Parameter learning**: the repository's `parameter_by_trek`, maximum likelihood of the linear Gaussian model from the observed covariance matrix (latent noise variances fixed at 1), from random starts. Before each repeat `torch.manual_seed(seed)` with seed = 0…N−1.
3. **Latents**: for each child c of latent L, first remove the contributions of c's other parents, r_c = c − Σ a_pc·p; then L = Σ a_Lc·r_c / Σ a_Lc². Meaning: the L value minimizing the sum of squared differences between "children predicted from L" and the actual children.
   - a_Lc: coefficient of the edge from L to c; a_pc: coefficient from c's other parent p to c
   - start from the lowest latent (L2 first, then use L2 as a child of L1)

   This is equivalent to the per-star Nelder-Mead minimization of chisq in the authors' notebook; checked with random data and the notebook's coefficients, identical to 6 decimal places.
4. **Conversion**: np.polyfit linear fit f(L1) = α·L1 + β to Rbirth, f(L2) to Rg; reports Pearson r, residual standard deviation, median absolute deviation. The sign of a latent is arbitrary, so the sign of the slope does not matter.

**Outputs**: `results/rlcd_A/<version>/<ID>/` with `nll_repeats.csv`, `best_params.dot/.png`, `latents.csv` (iord, L1, L2, f_L1, f_L2, R_bir, Rg, age), `fit.txt`; summaries `results/rlcd_A/<version>/summary.csv/.txt`; figures `L1_vs_R_bir.png`, `L2_vs_Rg.png`.

**The first trial (v1: 10 repeats, smallest NLL reported by the code) failed; cause found (2026-09-26):**
- g279 gave r(L1, Rb) = 0.08 (paper 0.64); 538905 gave r(L1, Rb) = 0.08, r(L2, Rg) = 0.25.
- **Latent computation and conversion are correct**: with the coefficients of the authors' notebook and our formulas, g279 gives exactly Rb = −1.969·L1 + 9.309 (paper −1.97, 9.31), r = 0.640; Rg = −1.289·L2 + 9.518 (paper −1.29, 9.52), r = 0.976.
- **The problem is which parameter-learning run is selected**:
  1. `parameter_by_trek` computes the NLL with `torch.slogdet(...)[1]`, the log of the **absolute value** of the determinant, and does not check that each node's implied noise variance is ≥ 0. Physically impossible solutions can therefore get a "small" NLL. The theoretical minimum NLL of g279 is −1.275 (= ½·logdet(observed correlation matrix) + 3), but repeat 5 gave −7.03, which is impossible. v1 selected exactly such invalid solutions (seed 4 for g279, seed 5 for 538905).
  2. Even among valid solutions, the smallest NLL is not what the authors used: their notebook and paper use the repository's **second smallest NLL (1.586)**, not the smallest (1.498). The 1.498 set almost disconnects L1 from the chemistry (L1→MgH = 0.012), giving |r(L1, Rb)| of only 0.01. Reason: [Fe/H] and [Mg/H] correlate at 0.997, which dominates the NLL. **The 1.586 set has the smallest RMSE (difference between model and observed correlation matrices), 0.117**; the 1.498 set has 0.260.
  3. Parameter learning is very sensitive to the random start: for the same galaxy, valid solutions give |r(L1, Rb)| anywhere from 0.01 to 0.6.
- **Changes in v2** (the current `rlcd_A.py`):
  - after each parameter-learning run, `implied_cov` computes the model correlation matrix and every node's noise variance; invalid solutions (negative noise variance or non-positive-definite matrix) are removed and the correct NLL is recomputed;
  - among valid solutions, the smallest RMSE is selected by default (`--select rmse`, which picks the authors' set); `--select nll` is also possible;
  - 50 repeats by default (about 2 s each);
  - |r(L1,Rb)| and |r(L2,Rg)| of every repeat are written to nll_repeats.csv **for diagnostics only, never for selection**; fit.txt lists |r(L1,Rb)| of the top-5 valid solutions to show the stability.
- This is the second problem not explained in the paper or repository (the first being the randomness of structure learning): parameter learning has many local optima, and the authors did not use the smallest-NLL solution.

**v2 rerun and the deeper problem: with the authors' graph frozen the coefficients are not identifiable (2026-09-26)**
- v2 rerun: g279 selects seed 31 (RMSE 0.108), r(L1, Rb) = 0.626, r(L2, Rg) = 0.86; 538905 selects seed 46 (RMSE 0.078), r(L1, Rb) = 0.056, r(L2, Rg) = 0.78. But the top-5 valid solutions by RMSE give |r(L1,Rb)| = 0.63/0.01/0.31/0.01/0.11 (g279) and 0.06/0.20/0.52/0.44/0.19 (538905): extremely unstable.
- Further checks with `scripts/local/rlcd_identifiability.py` (numpy/scipy; results in `results/rlcd_A/identifiability/`):
  1. **Jacobian rank = 14 < 15.** The authors' graph has 15 edges and the 6 observables have exactly 15 correlations; the derivative of the map "15 coefficients → 15 correlations" has only 14 independent directions. There is a direction along which the coefficients can change together while the observed correlation matrix stays **exactly the same** — mainly L1→age vs L1→MgH, compensated by the edges from age and MgH to ecc, JPhi, FeH. A structural property of the graph, the same for g279 and TNG50.
  2. **Maximum likelihood with a properly converging optimizer (L-BFGS, multiple starts)**: 14 solutions on g279 reach the theoretical optimum NLL (within 0.01), with |r(L1,Rb)| from 0.06 to 0.61; 25 such solutions on 538905, with |r(L1,Rb)| from 0.14 to 0.45.
  3. **Profile**: fixing L1→age at values in −0.95…0.95 and refitting the other 14: on g279 the NLL stays almost everywhere at the theoretical minimum −1.275 while |r(L1,Rb)| moves continuously from 0.00 to 0.63 (near L1→age ≈ 0, L1 resembles [Mg/H] most and r peaks at about 0.68, with a slightly worse NLL). On 538905 |r(L1,Rb)| moves from 0.07 to 0.43.
- **Conclusion**: with the authors' graph frozen, **the data alone cannot determine what L1 is**. The correlation of L1 with the birth radius depends on where on this curve of equivalent solutions the optimizer happens to stop. The paper's r = 0.64 corresponds to an **unconverged** solution (NLL 1.586 vs theoretical optimum −1.275), sitting at the "L1 ≈ [Mg/H]" end. Approach-A results obtained as in the paper cannot be used as evidence that the method works.
- This is the most important finding of the project (the third and most fundamental problem). How to proceed: Sec. 7, item 9.

**Results for the 12 galaxies and the strict re-analysis (2026-09-26, option (a))**
- Round 1, `rlcd_identifiability.py`: for the 13 data sets (g279 + 12 nonoise; likewise noise) the Jacobian rank is always 14/15. With "ΔNLL < 0.01", |r(L1,Rb)| ranges mostly from ~0 to 0.4–0.7 (`summary_{nonoise,noise}.*`).
- **This criterion is too loose**: the NLL is per star, so "indistinguishable" should be 2·N·ΔNLL < 3.84. Moreover the round-1 solutions were not strictly converged (NLL errors ~1e-4) and used a noise-variance floor of 0.01.
- Round 2, `rlcd_identifiability_ci.py --extra 200`: the 30 round-1 starts + 200 new starts, all re-optimized with "noise variance ≥ 0" (penalty 1e7) and strict convergence (ftol 1e-15, gtol 1e-11), keeping solutions with 2·N·ΔNLL < 3.84. Results (`summary_ci_{nonoise,noise}.*`, numbers in `REPORT.md` Sec. 5.6):
  - the optimum almost always pushes the noise variance of one observable to 0: nonoise J_φ in 10/12 (the optima of 546870 and 560421 are interior, smallest noise variances 0.04 and 0.03); noise J_φ 9, [Fe/H] 2, e 1; g279 [Fe/H] or e (Heywood case → the linear model does not match the data).
  - TNG50: in most galaxies indistinguishable solutions still give |r(L1,Rb)| from ~0 to 0.4–0.5; narrower for 546870.
  - g279: only 1 solution under the strict criterion, but it changes with the search: with the 30 round-1 starts the best solution has |r(L1,Rb)| = 0.25, |r(L2,Rg)| = 0.76; with 230 starts the best (NLL lower by 0.0009) has 0.62 and 0.55.
  - the paper's coefficients: NLL 1.593 with the correct formula (optimum −1.272), likelihood-ratio statistic ≈ 6.8e5.
  - `linear_ceiling.py`: linear-regression ceiling of the 6 quantities 0.53–0.85 (TNG50); age alone 0.23–0.69, [Fe/H] alone 0.00–0.46; g279 [Fe/H] alone 0.70. **In all 12 nonoise galaxies the most favourable L1 does not beat the better of age and [Fe/H]; noise has 2 exceptions (560421: 0.58 vs 0.57; 567607: 0.71 vs 0.67).**
- `results/rlcd_A/README.md` updated accordingly (including the correction note).

---

## 7. Open decisions

1. **Eccentricity**: 3D r (current) or cylindrical R?
2. **Additional independent checks?** Energy conservation along orbits; J_z from orbit integration compared with the Stäckel approximation.
3. **(Decided, Sec. 6.5) Selection** (1.6). Points that were considered:
   - radial range: the paper says 7–10 kpc but the repository data are at 9.4–12.1 kpc, possibly scaled with galaxy size. TNG50 galaxies have v_c ≈ 160–190 km/s, smaller than the Milky Way.
   - height cut: e.g. |z| < 1 kpc or a z_max limit
   - low-α cut: TNG50 median [Mg/Fe] ≈ +0.3; decide per galaxy from the [Mg/Fe]–[Fe/H] distribution, or no split
   - [Fe/H] lower limit (−0.5 or −1)
   - only stars with a birth radius (has_birth = True), and Rbirth < 20 kpc?
   - exclude counter-rotating (J_φ < 0) or unbound stars?
4. **(Decided) Observational errors**: both versions generated (nonoise main, noise as in the paper).
5. **(Decided, Sec. 6.5) RLCD input format**: the repository's 6 columns `logJZ, JPhi, ecc, age, FeH, MgH`, one csv per galaxy (`data/rlcd_input/<version>/<ID>_e.csv`); truth in `<ID>_truth.csv` with `iord, R, z, Rg, R_bir, ...`. J_z = 0 gives log₁₀ = −inf: such stars are removed (J_z > 0).
6. **(Decided) Python version for RLCD**: a separate 3.12 environment.
7. **(Decided, see 8) How to test**: learn a full RLCD structure on every galaxy, or, as in the paper, learn the structure on one galaxy, freeze it and re-estimate only the coefficients on others.
8. **(Decided 2.2b) Which causal graph for TNG50**:
   - A. As in paper Secs. 4.2, 4.3: freeze the authors' final g2.79 graph (15 edges), re-estimate only the coefficients on each TNG50 galaxy, compute L1, L2 (the repository's `*_with_sim_dag` scenarios do this). Most directly comparable with the paper, and independent of reproducing structure learning exactly.
   - B. Run structure learning from scratch on each TNG50 galaxy (e.g. α0 = 1e-5) and see whether two latents are found independently. A more direct test of the method, but results depend on α.
   - Decided: both, A as the main result and B as a supplement.
9. **(Decided: (a)) The coefficients are not identifiable — how to continue Approach A?** Options:
   - (a) report the non-identifiability itself as the main result: for each TNG50 galaxy give the possible range of |r(L1, Rbirth)| (best to worst) and show that the paper's 0.64 is one point on this curve;
   - (b) remove the free direction with an additional constraint (e.g. fix L1→age or L1→MgH, or require a relation between L1 and some quantity), which must be physically motivated and must not use the truth;
   - (c) add observables (e.g. another element abundance or J_R) so that there are more correlations than edges and the model becomes identifiable; this changes the graph and is closer to Approach B;
   - (d) pause Approach A and do Approach B first, to see which graph the TNG50 data produce and whether it is identifiable.

---

## 8. Installation notes (redo these if something breaks)

### 8.1 AnastrisTNG (optional)
The local scripts do not need AnastrisTNG; the server extraction uses the old version installed on TNG JupyterLab. To install 2.0.0 locally: `cd AnastrisTNG && $PY -m pip install -e .`

### 8.2 AGAMA
See [`docs/INSTALL_AGAMA.md`](INSTALL_AGAMA.md) (build from source with Homebrew's libomp; the PyPI package is only a bootstrap).

### 8.3 RLCD
```bash
brew install python@3.12 graphviz                       # macOS; use your package manager elsewhere
python3.12 -m venv .venv-rlcd                           # any location; then PY_RLCD=.venv-rlcd/bin/python3
.venv-rlcd/bin/python3 -m pip install --upgrade pip
.venv-rlcd/bin/python3 -m pip install -r setup/requirements-rlcd.txt
bash tools/setup_rlcd_repo.sh                           # clones Causal-Galactic-Archaeology at commit 8c0a5a4
```
Test: reproduce the paper's causal graph with the repository's g2.79e12 data
```bash
cd Causal-Galactic-Archaeology/rlcd
PYTHONHASHSEED=1 $PY_RLCD run_rlcd.py --s causal_vertical_heating_g279_e \
    --alpha0 1e-5 --alpha1 1e-15 --alpha2 1e-5 --alpha3 1e-15 2>&1 | tee ../../results/rlcd_g279/<log name>.log
```
(This overwrites the precomputed `results_causal_vertical_heating_g279_e/result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot/.png` in the repository; compare with `git diff` / `git status`, then move the new files to `results/rlcd_g279/` and restore the originals with `git checkout`. `results/reproduce_Jin/causal_graph/jin_causal_graph.py` does the same in memory without touching the repository.)

---

## 9. Known pitfalls

- **RLCD results are random** (set order + unseeded random numbers): the same parameters can give different graphs. Always set `PYTHONHASHSEED=<number>` and repeat with several seeds.
- **"Equally good" must be statistical**: 2·N·ΔNLL < 3.84, not a fixed ΔNLL threshold; polish (strictly converge) before comparing.
- Running RLCD inside `Causal-Galactic-Archaeology/rlcd/` overwrites the authors' committed results; move outputs out and restore with git.
- `grid_search.py` calls `python`: activate the RLCD environment first.
- Two analysis processes at once: set `OMP_NUM_THREADS=1`.
- Reading HDF5 without h5py: the scripts fall back to the pure-Python reader pyfive (`pip install pyfive`); pyfive datasets have no `len()`, use `.shape[0]`.
- `TNG50_disk_trajectory.hdf5`: at early snapshots the alignment switches to baryon/all, so the "disk plane" of the oldest stars is less reliable; `birth_R_3D` can be used for early times.
- On the server, `TNG50_disk_trajectory.hdf5` must be in the working directory (symlink), otherwise every galaxy fails with "No such file".
- The `script` attribute of the downloaded `data_obtain_all.hdf5` says `001_orbit_causal.py` (old name of `data_obtain.py`).
- Reference run times: identifiability analysis of 12 galaxies, 1–2 hours per round on 2 cores.

---

## 10. Timeline (short)
- 2026-09-26: read the paper and repository → updated local AnastrisTNG → designed the server extraction → ran data_obtain.py and merge_hdf5.py on the server, downloaded → cross-check with earlier results, consistent → AGAMA installed (after three failures) → compute_orbits.py test passed → full run failed for 530330 → inner grid boundary rmin = 0.2 kpc → rerun → RLCD installed → 4 α combinations → RLCD found to be random; seed 1 of 5 reproduces the authors' graph → orbits done → checks 1.5 passed → selection 1.6 (7–10 kpc, |z|<0.3, [Fe/H]>−0.5), RLCD inputs → **stage 1 done** → Approach A: rlcd_A v1 failed (invalid solutions selected) → v2 → non-identifiability found (Jacobian rank 14) → option (a) → round 1 for 12 galaxies nonoise/noise → ΔNLL<0.01 criterion found too loose → strict re-analysis (rlcd_identifiability_ci.py) → report written
- 2026-09-27: 30-seed structure frequencies; figure scripts (reproduction of paper Figs. 1, 3; TNG50 equivalents); wider selections R_0_30_Z_30, R_0_30_Z_3, R_0_30_Z_0.3; Approach B (30 seeds per galaxy); report updated
- 2026-10-21: packaged as a public GitHub repository (English documentation, public main-selection data, helper scripts)
