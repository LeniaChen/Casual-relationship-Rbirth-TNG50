# Testing "causal discovery of the stellar birth radius" (RLCD) with TNG50 — full report

- All paths below are relative to the root of this repository. A short summary is in [`README.md`](README.md); the current status and open questions are in [`STATUS.md`](STATUS.md).
- Dates: analysis 2026-09-26/27 (structure-learning seed frequencies, figure scripts, three wider star selections 0 < R < 30 kpc with |z| < 30 / 3 / 0.3 kpc and the corresponding conclusions added on 2026-09-27).
- Companion document: [`docs/REPRODUCTION.md`](docs/REPRODUCTION.md) (command-by-command reproduction log, installation notes, pitfalls). This report is for understanding; that one is for reproducing.
- Notation: every formula is followed by a plain-language explanation and a definition of each symbol.

---

## 0. One-page summary

1. **What was done.** For 12 Milky-Way-like galaxies from TNG50 we computed, following the paper, six "observables" (age, [Fe/H], [Mg/H], J_z, J_φ, eccentricity e) and two "truths" (birth radius Rbirth, guiding radius Rg). We installed the authors' RLCD code, reproduced their result, then applied their causal graph to TNG50 to see whether the latent variables L1, L2 it finds really correspond to Rbirth and Rg.
2. **Three problems that neither the paper nor the code repository mentions:**
   - Problem 1: RLCD structure learning is **random**. With identical parameters, a different random seed can give a different graph. The authors' graph appeared in only 7 of our 30 repeats (23%); the most common result was a different 15-edge graph (12 times).
   - Problem 2: the repository's way of picking the "best" edge coefficients is flawed (physically impossible solutions can look good). The coefficients used in the paper are **far from the maximum-likelihood solution**: on the authors' own 119,000 stars, the likelihood-ratio statistic between them and the optimum is about 680,000 (anything above 3.84 counts as "significantly worse").
   - Problem 3 (the most fundamental): **with the authors' causal graph frozen, the edge coefficients cannot be determined reliably from the data.** (i) It is mathematically non-identifiable: of the 15 coefficients only 14 directions are constrained by the data, so there is a whole family of coefficient sets with exactly the same fit. (ii) A proper maximum-likelihood fit finds many almost equally good local optima, and the best solution almost always pushes the noise variance of one observable (J_φ in most galaxies) to zero, a sign that the linear model does not describe the data. (iii) The meaning of L1 is therefore unstable: in the 12 TNG50 galaxies, solutions statistically indistinguishable from the best give |r(L1, Rbirth)| anywhere from about 0 to 0.4–0.5; on the authors' own galaxy, going from 30 to 230 random starts changed the best solution's |r(L1, Rbirth)| from 0.25 to 0.62.
3. **Conclusions.**
   - **Applying RLCD to TNG50 as done in the paper does not reliably recover the birth radius.** First, how well L1 tracks Rbirth depends on where the optimizer stops, not on the data. Second, **even the most favourable solution essentially never beats "age alone" or "[Fe/H] alone"**: never in the 12 galaxies without observational errors; with errors only 2 galaxies are marginally higher at the favourable end (0.58 vs 0.57, 0.71 vs 0.67), while the same galaxies also admit L1 solutions as low as 0.01–0.09 (Sec. 5.6).
   - **Widening the star selection does not rescue it** (Sec. 5.9). With the whole disk (0 < R < 30 kpc, |z| limits of 30, 3 and 0.3 kpc; 0.19–0.72 million stars per galaxy), |r(L1, Rbirth)| rises to about 0.8 and f(L1) spreads over 0–20 kpc like paper Fig. 3. But this is only because across the whole disk "where a star is now ≈ where it was born": **J_φ alone (≈ current radius) already reaches 0.8–0.9 in most galaxies.** In 36 galaxy × selection combinations (no errors), the maximum-likelihood L1 beats the best single observable in only 5, by at most 0.07; and the coefficients are still not unique: changing only the |z| limit moves the maximum-likelihood |r| of galaxy 546870 from 0.71 to 0.09.
   - The correlation of L2 with the guiding radius Rg is more stable (|r| = 0.6–0.95 in most galaxies), but this means little: Rg is essentially J_φ divided by a constant, so anything tied to J_φ correlates with Rg.
   - **The method is not ready to be applied to Milky Way data.** Problem 3 is a property of the causal graph and exists for Milky Way data as well, and without a true birth radius there is no way to tell which kind of L1 one has obtained. Before applying it to the Milky Way the model must be made identifiable and shown, on TNG50, to beat single observables (Secs. 6, 7).
4. **Correction.** An early version of this analysis stated that among "equally good fits" |r(L1, Rbirth)| could range from 0 to 0.6–0.7, using the criterion "NLL difference < 0.01". For tens of thousands of stars that criterion is far too loose. This report uses a strict statistical criterion (Secs. 3.9, 8.4.2): the TNG50 conclusions are essentially unchanged and the ranges slightly narrower; the authors' galaxy g2.79e12 is more subtle (point (iii) above). Both sets of numbers are kept in the results folder.

---

## 1. What this project is

### 1.1 Background: why the birth radius
After stars are born in a disk they **migrate radially**: the radius where a star is now is not the radius where it was born. The birth radius Rbirth sets the chemical environment at birth and is a key quantity for reconstructing the formation history of the Milky Way (Galactic archaeology), but it **cannot be observed** and must be inferred indirectly.

### 1.2 The work being tested
Jin, Lu, Ting, Zheng & Buck, *Causal Discovery of Latent Variables in Galactic Archaeology* (ApJ; [arXiv:2507.00134](https://arxiv.org/abs/2507.00134)); code: <https://github.com/ZehaoJin/Causal-Galactic-Archaeology>.
They use an algorithm called **RLCD** (Rank-based Latent Causal Discovery, Dong et al. 2024) which, given only six observables, searches for hidden ("latent") variables behind them. On the NIHAO simulated galaxy g2.79e12 it finds two latents: L1 correlates with the birth radius (r = 0.64) and L2 almost one-to-one with the guiding radius. They then freeze this causal graph and apply it to another simulated galaxy and to the real Milky Way (APOKASC-3).

### 1.3 Our goals
1. **Test the method on an independent simulation (TNG50)**, where the true birth radius of every star is known and can be compared directly.
2. If it works, consider applying it to Milky Way observations (the authors have done this once; one would decide whether to ask them for data or redo it from scratch).

---

## 2. Reference files and tools

### 2.1 The paper and the authors' repository
| File / folder | Content | Role |
|---|---|---|
| the paper ([arXiv:2507.00134](https://arxiv.org/abs/2507.00134)) | main text | basis for the method and conclusions (details in Sec. 8) |
| `Causal-Galactic-Archaeology/data/g2.79e12_e.csv` | the 6 input columns for g2.79e12: `logJZ, JPhi, ecc, age, FeH, MgH`, 119,236 stars | the data RLCD actually reads; our TNG50 inputs use the same format |
| `.../data/g2.79e12_multi_ele_action_Jin_SN.csv` | truth for the same stars: `age, iord, JR, JZ, JPhi, FeH, MgFe, ecc, Rg, R_bir, z_bir, R` | to check the answers |
| `.../data/g8.26e11_*.csv`, `obs_e.csv`, `APOKASC3_*.csv` | second simulated galaxy, Milky Way data | used in paper Secs. 4.2, 4.3 |
| `.../rlcd/run_rlcd.py` | main program: read data → learn the causal graph (structure learning) | writes `.dot` (text description of the graph) and `.png` |
| `.../rlcd/scenarios.py` | one function per data set (reads the csv, prefixes columns with `X_`) | new data must be registered here |
| `.../rlcd/grid_search.py` | loops over 4 significance levels α, each 1e-5/1e-10/1e-15, 81 combinations | the "grid search" of paper Sec. 3.2 |
| `.../rlcd/StructureLearning/RLCD/` | the RLCD algorithm (`RLCD_alg.py`, `LatentGroups.py`, `Cover.py`, `misc.py`, ...) | structure learning |
| `.../rlcd/utils/Chi2RankTest.py` | rank test (is the rank of a covariance block ≤ k?) | the core statistical test of structure learning |
| `.../rlcd/ParameterLearning/linear_parameter_identification.py` | `parameter_by_trek`: maximum-likelihood edge coefficients for a given graph | parameter learning |
| `.../rlcd/results_causal_vertical_heating_g279_e/` | the authors' causal graph and 13 parameter-learning results (the number in each file name is the NLL) | reference for the reproduction |
| `.../L1L2calculation/*.ipynb` | per-star L1, L2 from **hard-coded** edge coefficients, then linear conversion to Rbirth, Rg | paper Figs. 3–6 |

### 2.2 Scripts from an earlier project (by the author of this repository)
| File | Content | Role |
|---|---|---|
| `scripts/reference/001_orbit.py`, `001_orbit_new.py` | earlier birth-radius computation (13 galaxies) | basis of the new `data_obtain.py` |
| `scripts/reference/TNG50_disk_trajectory.py` | builds the main-progenitor trajectory file | defines the trajectory file |
| `scripts/reference/obtain_data_s98.py` | earlier data extraction | reference |
| `TNG50_disk_trajectory.hdf5` (not distributed; made on the TNG server by `scripts/reference/TNG50_disk_trajectory.py`) | centre, velocity and angular-momentum direction of 423 galaxies at every snapshot | transforms birth positions into the "disk frame at birth" |
| `001_orbit/<ID>.txt` (not distributed) | earlier birth-radius results | cross-check of the new results |

### 2.3 Software
| Software | Purpose | Where |
|---|---|---|
| TNG JupyterLab + old AnastrisTNG, galpos, pynbody | read TNG50 particles, compute birth radii | TNG server |
| AGAMA 1.0 | potential fitting, actions, orbit integration | analysis environment (built from source, [`docs/INSTALL_AGAMA.md`](docs/INSTALL_AGAMA.md)) |
| dependencies of the RLCD repository (Python 3.12, causal-learn, torch, statsmodels, ...) | run RLCD | RLCD environment (separate Python 3.12 environment to avoid conflicts; [`setup/requirements-rlcd.txt`](setup/requirements-rlcd.txt)) |
| Graphviz | draw causal graphs | system package |

---

## 3. What was done, step by step

Overall flow:
```
TNG50 server          local (AGAMA)                        local (RLCD)
particles ──► data_obtain.py ──► compute_orbits.py ──► make_rlcd_input.py ──► reproduce RLCD ──► Approach A ──► identifiability tests
            (age/abundances/Rbirth)  (potential/actions/Rg/e)  (selection + 6 inputs + truth)   (authors' data)   (TNG50)    (12 galaxies)
```

| Step | Script | Input | Output |
|---|---|---|---|
| 3.1 extract particles, birth radii | `scripts/server/data_obtain.py`, `merge_hdf5.py` (server) | TNG50-1 snap 99, `TNG50_disk_trajectory.hdf5` | `data/tng50_extract/data_obtain_all.hdf5` (0.96 GB) |
| 3.2 cross-check with earlier results | `scripts/local/check_vs_001_orbit.py` | the file above, earlier `001_orbit/*.txt` | `results/check_vs_001_orbit.txt` |
| 3.3 potential, actions, orbits | `scripts/local/compute_orbits.py` | `data_obtain_all.hdf5` | `data/tng50_derived/orbits_all.hdf5`, `potentials/*.coef`, `results/compute_orbits.log` |
| 3.4 check the orbit results | `scripts/local/check_orbits.py` | the two hdf5 files | `results/orbit_check/` |
| 3.5 star selection, RLCD inputs | `scripts/local/make_rlcd_input.py` | the two hdf5 files | `data/rlcd_input/{nonoise,noise}/<ID>_e.csv`, `<ID>_truth.csv`, `results/rlcd_input_selection.txt` |
| 3.6 install RLCD | `setup/requirements-rlcd.txt` | authors' repository | RLCD environment |
| 3.7 reproduce the authors' graph | repo's `run_rlcd.py` + `scripts/local/compare_rlcd_graphs.py` | authors' g2.79e12 data | `results/rlcd_g279/` |
| 3.8 Approach A: freeze the authors' graph, re-estimate coefficients | `scripts/local/rlcd_A.py` | authors' graph, TNG50 inputs | `results/rlcd_A/g279/`, `results/rlcd_A/nonoise/538905/` |
| 3.9 identifiability tests | `scripts/local/rlcd_identifiability.py`, `rlcd_identifiability_summary.py` (loose), `rlcd_identifiability_ci.py` (strict) | authors' graph, all inputs | `results/rlcd_A/identifiability/` |
| 3.10 linear ceiling (reference) | `scripts/local/linear_ceiling.py` | all inputs and truths | `results/rlcd_A/identifiability/linear_ceiling.*` |
| 3.11 figures (reproduce paper Figs. 1 and 3; TNG50 equivalents) | `results/reproduce_Jin/{causal_graph,truth_comparison}/*.py`, `results/TNG/{causal_graph,truth_comparison}/*.py` | authors' data and coefficients; results of 3.9 | png next to each script |
| 3.12 repeat 3.9 with wider selections | in `star_condition/<selection>/`: `make_input.py` → `compute_moments.py` → `identifiability_moments.py` → two plotting scripts | `data_obtain_all.hdf5`, `orbits_all.hdf5` | `star_condition/R_0_30_Z_30/`, `R_0_30_Z_3/`, `R_0_30_Z_0.3/` |

### 3.1 Server-side extraction: `scripts/server/data_obtain.py`
**Sample**: 12 TNG50-1 galaxies (Subhalo IDs at snap 99): 487742, 530330, 538905, 545437, 546870, 547293, 552581, 557721, 560421, 562029, 567607, 569251. Galaxy 607654 from the original list is not in the trajectory file of the 423-galaxy sample and was dropped.

**Why only this step runs on the server**: birth positions and abundances require reading full simulation snapshots (several TB), possible only on the server; the orbit computation needs repeated debugging and is easier locally.

For each galaxy:
1. **Read particles**: dark matter, stars, gas, black holes. For stars also the formation time `GFM_StellarFormationTime`, birth position `BirthPos`, initial mass and metal abundances `GFM_Metals`; convert to physical units (kpc, km/s, Msun).
2. **Remove wind particles**: "stars" with formation time ≤ 0 are wind particles in TNG. They are excluded from the analysis but their mass is kept for the potential.
3. **Age** = present age of the Universe (13.80 Gyr) − age of the Universe at formation. TNG stores the formation time as the scale factor a_form, converted to time with the cosmological model.
4. **Abundances** (solar values from Asplund et al. 2009):
   - [Fe/H] = log₁₀[(X_Fe/A_Fe) / (X_H/A_H)] + 4.50
     Meaning: the ratio "number of iron atoms / number of hydrogen atoms" in the star, relative to the same ratio in the Sun, on a log scale. 0 means solar, −0.3 means half the solar iron abundance.
     - X_Fe, X_H: mass fractions of iron and hydrogen (columns 8 and 0 of `GFM_Metals`)
     - A_Fe = 55.845, A_H = 1.008: atomic masses, converting a mass ratio into a number ratio
     - 4.50: 12 − log ε_Fe of the Sun (the solar iron abundance constant)
   - [Mg/Fe] = log₁₀[(X_Mg/A_Mg) / (X_Fe/A_Fe)] − 0.10: the same, for magnesium relative to iron. X_Mg is column 6 of `GFM_Metals`, A_Mg = 24.305, and 0.10 is the difference between the solar Mg and Fe abundance constants.
   - [Mg/H] = [Mg/Fe] + [Fe/H] (computed when the inputs are generated)
5. **Birth radius (truth)**: `BirthPos` is the position in the simulation box at the moment of birth, when the galaxy centre and disk orientation differed from today's. Using galpos' `GalaxyPoseTrajectory`, the centre, velocity and angular-momentum direction of the main progenitor at every snapshot are read from the trajectory file and interpolated to the birth time of each star, transforming its birth position into the "disk frame at birth" (origin at the galaxy centre at that time, z axis along the disk angular momentum at that time). Then:
   - Rbirth = `birth_R_2D` = √(x_b² + y_b²)
     Meaning: the distance of the birth position from the centre within the disk plane at that time (cylindrical radius).
     - x_b, y_b: the two in-plane coordinates in the disk frame at birth
   - also stored: `birth_R_3D` = √(x_b² + y_b² + z_b²) (straight-line distance to the centre; z_b is the height above the disk plane at birth).
   - The few old stars born before the start of the trajectory file get NaN (`has_birth = 0`); stars whose formation time is later than the end of the trajectory by < 0.001 Gyr because of floating-point rounding are clamped to the end.
6. **Rotate the present-day galaxy face-on**: `face_on(rmax=8, alignwith='star')`, i.e. the z axis along the total angular momentum of the stars within 8 kpc. All positions and velocities are then in this "z=0 disk frame".
7. **Write files**: one hdf5 per galaxy, merged by `merge_hdf5.py` into `data_obtain_all.hdf5` for download. Structure: `/<ID>/stars` (iord, pos, vel, mass, age, FeH, MgFe, birth_pos, birth_R_2D, birth_R_3D, has_birth, ...), `/<ID>/dm`, `/gas`, `/bh`, `/wind` (pos and mass only, for the potential).

### 3.2 Cross-check with earlier results: `scripts/local/check_vs_001_orbit.py`
Star by star (matched on particle ID iord), the new file is compared with the earlier `001_orbit` results: ages, masses and positions agree to < 1e-5, birth radii to ≤ 5 pc. The check also revealed that `Birth_R` in the old files is a 3D distance, corresponding to `birth_R_3D` in the new file. Report: `results/check_vs_001_orbit.txt`. **Conclusion: the new extraction is correct.**

### 3.3 Potential, actions and orbits: `scripts/local/compute_orbits.py`
Units: Msun, kpc, km/s. 3.8–8 minutes per galaxy.

**Step 1: fit the gravitational potential.** "Fitting the potential" means turning the mass distribution of millions of particles into a smooth potential Φ(R, z), so that forces can be evaluated anywhere and orbits integrated. As in the paper the potential is assumed **axisymmetric** (unchanged under rotation about the z axis). It has two parts that are added:
- **disk part**: AGAMA `CylSpline` (a 30×30 grid in the R–z plane, suited to flat disks): stars with r < 50 kpc plus gas with R < 30 kpc and |z| < 5 kpc.
- **halo part**: `Multipole` (angular expansion to order 8, 40 radial nodes, suited to nearly spherical halos): all dark matter, black holes, wind particles, and the remaining stars and gas.
- The innermost grid node is at **0.2 kpc** (about the TNG50 gravitational softening). Initially AGAMA chose the grid automatically and went down to 0.01 kpc, where there are too few particles; the noise made galaxy 530330 fail (Sec. 4.4).

**Step 2: orbit quantities of every star** (bound stars only, E < 0, about 99.4–100%):
- **Energy** E = Φ(x) + ½|v|²
  Meaning: total energy per unit mass = potential energy at the star's position + kinetic energy. E < 0 means bound to the galaxy.
  - Φ(x): potential at the star's position; |v|: speed
- **Angular momentum** J_φ = x·v_y − y·v_x
  Meaning: the angular momentum about the z axis (L_z); positive means rotating with the disk.
  - x, y: the two in-plane coordinates; v_x, v_y: the corresponding velocity components
- **Actions J_R, J_z**: AGAMA `ActionFinder` (Stäckel fudge, the same family of methods as in the paper and galpy).
  In words: J_R measures the amplitude of the radial oscillation, J_z that of the vertical oscillation through the disk; they are nearly conserved when the potential changes slowly, so they serve as an orbit's "identity". The input uses log₁₀(J_z) (`logJZ`) because J_z spans several orders of magnitude.
- **Guiding radius** Rg: the radius where |J_φ| = Rg · v_c(Rg) (AGAMA `Rcirc`).
  Meaning: the radius of the circular orbit with the star's current angular momentum.
  - v_c(R) = √(R · ∂Φ/∂R): the speed of a circular orbit at radius R
  - ∂Φ/∂R: the radial derivative of the potential, i.e. the strength of the radial force
- **Eccentricity** e: AGAMA integrates each orbit for 15 circular periods (1501 points per orbit) and finds the smallest distance r_peri and the largest r_apo:
  e = (r_apo − r_peri) / (r_apo + r_peri)
  Meaning: the flatter the orbit the closer e is to 1; a nearly circular orbit has e close to 0.
  - r_peri: pericentre distance; r_apo: apocentre distance (both 3D distances r, Sec. 4.4)
- z_max (the largest height above the plane along the orbit) and the rotation curve are also stored.

### 3.4 Checking the orbit results: `scripts/local/check_orbits.py` → `results/orbit_check/`
All 12 galaxies pass these consistency checks:
- the computed J_φ equals the directly computed L_z (difference < 0.004 kpc·km/s);
- every star's current position lies within its own [r_peri, r_apo], and |z| ≤ z_max (orbit integration is correct);
- the median of Rg·v_c(Rg)/|L_z| is 1.000 (Rg is correct);
- no NaN; all e within 0–1.
Also: v_c(8 kpc) = 185–229 km/s; the median eccentricity of disk stars (4 < R < 12 kpc, |z| < 1 kpc) is 0.18–0.32, "hotter" than the Milky Way thin disk; this is a property of the TNG50 disks, not an error.

| Figure | What it shows |
|---|---|
| `rotation_curves.png` | rotation curves of the 12 galaxies; all fairly flat and reasonable |
| `ecc_disk.png` | eccentricities of disk stars: low for young stars, high for old ones, as expected from gradual disk heating |
| `Rg_vs_R.png` | guiding radius vs current radius; most stars near Rg ≈ R; the horizontal band at Rg ≈ 0 is non-rotating bulge/halo stars, as expected |
| `logJz_disk.png` | distribution of log J_z for disk stars |
| `MgFe_FeH.png` | [Mg/Fe]–[Fe/H]: most TNG50 galaxies show **no clear high-α and low-α sequences** (the paper selects the low-α disk with [Mg/Fe] < −0.01, which cannot be copied to TNG50) |
| `summary.txt/.csv` | per-galaxy numbers |

### 3.5 Star selection and RLCD inputs: `scripts/local/make_rlcd_input.py`
**Selection**:
1. bound (E < 0), all orbit quantities defined;
2. has a birth radius, with Rbirth < 20 kpc (the authors' data appear to be cut at < 20 as well);
3. **currently at 7 < R < 10 kpc** (the range in the paper's main text);
4. **|z| < 0.3 kpc** (stars near the plane only);
5. **[Fe/H] > −0.5, no split into high/low α** (TNG50 has no clear α bimodality; the authors' repository data also use [Fe/H] > −0.5);
6. J_z > 0 (logarithm); counter-rotating stars kept (as in the authors' data).

The number of stars left after each step is in `results/rlcd_input_selection.txt`. In the end 7,299–34,970 stars per galaxy, 239,970 in total.

**Two versions**:
- `nonoise`: no observational errors (main version)
- `noise`: errors as in the paper: age × (1 + 0.10·N), [Fe/H] + 0.02·N, [Mg/Fe] + 0.06·N. N is a random number with mean 0 and standard deviation 1 (drawn independently for every star and quantity); the random seed is 20260926 + galaxy ID. The selection uses the error-free values, so both versions contain the same stars, row by row.

**Outputs**:
- `data/rlcd_input/<version>/<ID>_e.csv`: 6 columns `logJZ, JPhi, ecc, age, FeH, MgH`, exactly the format of the authors' `g2.79e12_e.csv`
- `data/rlcd_input/<version>/<ID>_truth.csv`: truth and other quantities for the same stars (`iord, R, z, Rg, R_bir, R_bir_3D, age, FeH, MgFe, JR, JZ, JPhi, ecc, rperi, rapo, zmax, E`, ...), row-matched

### 3.6 Installing RLCD
The repository requires Python 3.12 and pins many package versions that conflict with the AGAMA environment, so a separate Python 3.12 environment was created. The dependency list [`setup/requirements-rlcd.txt`](setup/requirements-rlcd.txt) is the repository's `py312_requirements.txt` without the Linux+GPU-only packages (nvidia-*, triton) and the backports for old Python versions. Graphviz is installed for drawing. Stage 1 uses causal-learn's GES by default and needs no Java.

### 3.7 Reproducing the authors' causal graph (authors' g2.79e12 data) → `results/rlcd_g279/`
1. **One run with the authors' parameters** (α = 1e-5/1e-15/1e-5/1e-15, the values in the authors' result file name): 2 latents with the same structure as the paper; all 15 of the authors' edges are present, plus 3 extra edges, one of which (MgH → age) violates the paper's prior (age may not have an observed parent).
2. **Three more α combinations** (`compare_4alpha.txt`): with α0 = 1e-5 the structure is close to the paper's; with α0 = 1e-10 or 1e-15 there are 4 latents and a completely different structure, which also satisfies the prior. So the paper's statement that "filtering by the prior converges to a unique graph" does not hold.
3. **Cause: five repeats with the same α** and different random seeds PYTHONHASHSEED = 0–4 (`compare_5seeds.txt`):
   - seed 1: **identical to the authors' graph** (the two panels of `compare_author_vs_seed1.png` are the same)
   - seeds 0, 2, 3: a different 15-edge graph. It looks similar, but (J_φ, e) and ([Fe/H], J_z) are attached to the opposite latents (`compare_author_vs_seed3_diff.png` colours the 4 differing edges), so the second latent no longer means "guiding radius"
   - seed 4: the seed-0 graph plus 3 edges
   → This is Problem 1 (Sec. 4.1). The installation itself is fine.
4. **30 seeds** (PYTHONHASHSEED = 0–29, same α): the authors' graph appears 7 times (seeds 1, 7, 8, 11, 14, 22, 25; 23%); the most common result is the seed-0 15-edge graph (12 times, 40%); the remaining 11 runs give four different 18-edge graphs.

### 3.8 Approach A: freeze the authors' graph and re-estimate only the coefficients: `scripts/local/rlcd_A.py`
This is what paper Secs. 4.2 and 4.3 do (freeze the graph learned on g2.79e12 and apply it to g8.26e11 and the Milky Way). Steps:
1. read the 6 inputs and standardize each column (subtract the mean, divide by the standard deviation), as in the repository;
2. **parameter learning**: call the repository's `parameter_by_trek`, which fits the 15 edge coefficients from random starting values; repeat 50 times;
3. remove physically impossible solutions with our own checks (Sec. 4.2), and among the rest pick the one with the smallest RMSE (root-mean-square difference between model and data correlation matrices) — a criterion that happens to select the coefficient set used by the authors;
4. **compute L1, L2 for every star** (hierarchical least squares, Sec. 8.3.6);
5. convert L1 to Rbirth and L2 to Rg with a linear fit and compute the Pearson correlation r.

**The procedure was first validated on the authors' data**: with the coefficients from the authors' notebook, our steps 4 and 5 give exactly Rbirth = −1.969·L1 + 9.309 (paper: −1.97, 9.31), r = 0.640, and Rg = −1.289·L2 + 9.518, r = 0.976. **Latent computation and conversion are exactly right.**

**But re-learning the coefficients is extremely unstable**: on g2.79e12 the selected solution gives r(L1, Rb) = 0.63, but the 5 best-fitting solutions give 0.63, 0.01, 0.31, 0.01, 0.11; on TNG50 galaxy 538905 they give 0.06, 0.20, 0.52, 0.44, 0.19 (`results/rlcd_A/*/fit.txt`). This led to Problem 3 and the tests of Sec. 3.9.

### 3.9 Identifiability tests: `rlcd_identifiability.py` → `rlcd_identifiability_ci.py`
Question: with the authors' graph frozen, are the 15 edge coefficients uniquely determined by the data? Two rounds (principles in Sec. 8.4).

**Round 1** (`scripts/local/rlcd_identifiability.py`, one run per galaxy):
1. **Count the equations**: the authors' graph has 15 edges (15 unknowns), and the 6 observables have 15 pairwise correlations (15 equations). That looks just enough, but one must check that the 15 equations are independent: compute the **rank of the Jacobian** of the map "15 coefficients → 15 correlations". At 3 random coefficient sets the rank is always **14**, not 15.
2. **Multistart fits**: maximum likelihood with a properly converging optimizer (L-BFGS-B) from 30 random starts, requiring every node's noise variance ≥ 0.01, to see whether the optimum is unique.
3. **Profile**: fix the coefficient of the single edge L1→age at 39 values between −0.95 and 0.95, refit the other 14 coefficients (8 starts per value), and see how the fit quality and the meaning of L1 change.

`rlcd_identifiability_summary.py` counts all solutions with "NLL within 0.01 of the best" as "equally good". **That criterion is too loose**: the NLL is averaged per star, and multiplied by tens of thousands of stars a difference of 0.01 is statistically distinguishable. Hence round 2.

**Round 2** (`scripts/local/rlcd_identifiability_ci.py`, strict criterion):
1. the 30 round-1 solutions plus 200 new random starts are all re-optimized with very strict convergence (polish), under the physical constraint "noise variance ≥ 0";
2. only solutions **statistically indistinguishable** from the best are kept: 2·N·ΔNLL < 3.84 (a 95% likelihood-ratio test, Sec. 8.4.2);
3. the range of |r(L1, Rbirth)| and |r(L2, Rg)| over these solutions is reported, together with the node whose noise variance is pushed to 0 at the optimum.

### 3.10 Linear ceiling (reference): `scripts/local/linear_ceiling.py`
L1 and L2 are **linear combinations** of the 6 inputs (Sec. 8.3.6), so whatever the coefficients, |r(L1, Rbirth)| cannot exceed the correlation obtained by regressing the true Rbirth directly on the 6 inputs. This regression uses the truth and is a reference only: it answers "how much birth-radius information is linearly available in these 6 quantities".

### 3.11 Figures: `results/reproduce_Jin/`, `results/TNG/`
Each folder contains only scripts; each script writes a png with its own name next to it:
| Script | What it draws | Result |
|---|---|---|
| `results/reproduce_Jin/causal_graph/jin_causal_graph.py` (RLCD environment) | reruns RLCD structure learning in memory (PYTHONHASHSEED = 1 by default, `--seed` to change), compares edge by edge with the authors' graph; edge labels are the authors' published coefficients (the nll1.586 set) | identical to paper Fig. 1. It is "identical" because seed 1 is fixed and the coefficients are the authors'; other seeds may differ (Sec. 3.7, item 4) |
| `results/reproduce_Jin/truth_comparison/jin_L1_Rbirth_L2_Rg.py` | authors' coefficients → L1, L2 per star → conversion; reproduces both panels of paper Fig. 3 | Rb = −1.969 L1 + 9.309, r = 0.640; Rg = −1.289 L2 + 9.518, r = 0.976 |
| `results/TNG/causal_graph/tng_causal_graph.py` (`--version noise`) | 12 TNG50 galaxies: authors' graph + each galaxy's maximum-likelihood coefficients; observables with noise variance ≈ 0 marked | Sec. 5.8 |
| `results/TNG/truth_comparison/tng_L1_Rbirth_L2_Rg.py` (`--version noise`) | same coefficients; f(L1)–Rbirth and f(L2)–Rg for the 12 galaxies, with the |r| range of indistinguishable solutions and the single-variable baselines in each panel | Sec. 5.8 |

(`star_condition/R_7_10/` contains copies of the two `results/TNG/` scripts and figures, so that all four selections have the same layout.)

### 3.12 Repeating the analysis with wider selections: `star_condition/`
**Why**: for 7 < R < 10 kpc the TNG50 f(L1) values form a very narrow vertical band, whereas f(L1) in paper Fig. 3 spans 0–20 kpc. Is the narrow selection the reason?

**Three selections** (all other cuts as in Sec. 3.5: bound, birth radius < 20 kpc, [Fe/H] > −0.5, J_z > 0, no [Mg/Fe] cut):
| Folder | Selection | Stars per galaxy | Total |
|---|---|---|---|
| `star_condition/R_0_30_Z_30/` | 0 < R < 30 kpc, \|z\| < 30 kpc (whole disk, almost no height cut) | 0.375–0.723 million | 6.17 million |
| `star_condition/R_0_30_Z_3/` | 0 < R < 30 kpc, \|z\| < 3 kpc | 0.371–0.699 million | 6.02 million |
| `star_condition/R_0_30_Z_0.3/` | 0 < R < 30 kpc, \|z\| < 0.3 kpc (same height cut as the paper) | 0.188–0.313 million | 2.79 million |

Each folder has five subfolders — `input_data` (RLCD inputs), `true_value` (truth), `output_data` (analysis results), `causal_graph`, `truth_comparison` — and a `README.md`. The last two selections were obtained by cutting the R_0_30_Z_30 files in |z| (the noisy version reuses the same errors).

**The method is exactly that of round 2 in Sec. 3.9** (authors' graph frozen, 230 starts, noise variances ≥ 0, strict convergence, 2·N·ΔNLL < 3.84), but the data are large (6 GB), so each galaxy is first compressed into a covariance matrix:
1. `output_data/compute_moments.py`: for each galaxy the means and the 8×8 covariance matrix of 8 quantities (6 inputs + true Rbirth + true Rg), saved as `moments_<version>.json`;
2. `output_data/identifiability_moments.py`: the full Sec. 3.9 analysis using only these covariance matrices. Why this is allowed is explained in Sec. 8.4.4: maximum likelihood uses only the covariance of the 6 inputs, and the correlation of L1 with the truth can also be computed from covariances. Checked on the 7–10 kpc data: identical to the star-by-star result to 6 decimal places;
3. the two plotting scripts (`causal_graph/causal_graph_<selection>.py`, `truth_comparison/L1_Rbirth_L2_Rg_<selection>.py`) work star by star; the r values in the figures equal those from the covariance method.

---

## 4. Problems found

### 4.1 Problem 1: structure learning is random (not mentioned in the paper or the repository)
- **Observation**: same data, same α, different PYTHONHASHSEED → different causal graphs; the authors' graph appears only 7 times in 30, and the most common result is another graph (12 times) (Sec. 3.7).
- **Cause** (confirmed in the source code):
  - RLCD stores variable names in Python sets and performs tests and latent merges in set-iteration order. Python randomizes string hashing at every start-up, so the iteration order, the order of the tests and hence the result can change.
  - `findRandomLatentChild` in `LatentGroups.py` calls `random.sample` without a seed.
  - The repository has a seeding function, `utils/utils.py: set_random_seed()`, but it is never called.
- **Consequence**: the authors' graph is one of several possible RLCD outputs, and not the most frequent one; the paper reports neither this nor the frequency.
- **Remedy**: fix PYTHONHASHSEED and the random seeds for every run; repeat with many seeds, count how often each structure appears, and report the most frequent structure and its stability.

### 4.2 Problem 2: the selection of the learned parameters is flawed (not mentioned in the paper or the repository)
- **Observation 1**: the repository's `parameter_by_trek` computes the likelihood with `torch.slogdet(...)[1]`, i.e. the log of the **absolute value** of the determinant of the covariance matrix, and never checks that each variable's "noise variance" is ≥ 0.
  In words: in the model every variable = contribution of its parents + its own noise, and the noise variance = 1 − the part explained by the parents, which must be ≥ 0. A coefficient set whose parents "explain" more than 100% of a variable's variance does not exist physically, but the code still assigns it a likelihood, possibly a "better" one. The theoretical minimum NLL on g2.79e12 is −1.275, yet −7.03 occurred in our repeats, which is impossible.
- **Observation 2**: the authors' repository contains 13 parameter-learning runs for the same graph, with the NLL in the file name (smaller is better). The smallest is 1.498, but the paper and notebook use the **second smallest, 1.586**, without explanation. Both are far from the optimum (about −1.272; theoretical lower bound −1.275), i.e. the optimizer (SGD, 600 steps) did not converge.
- **Consequence**: "taking the smallest NLL", as the repository does, can pick a physically impossible solution; our first attempt (v1) failed exactly this way (r = 0.08 on g2.79e12).
- **Remedy**: after each fit recompute the model correlation matrix and every node's noise variance, discard solutions with negative noise variance or a non-positive-definite matrix, and recompute the NLL with the correct formula (implemented in `rlcd_A.py` v2, which picks the smallest RMSE among the valid solutions); use a properly converging optimizer with many starts (`rlcd_identifiability_ci.py`). This fixes "invalid solutions and non-convergence", but not Problem 3.

### 4.3 Problem 3 (the most fundamental): with the authors' graph frozen, the edge coefficients cannot be determined reliably (not mentioned in the paper or the repository)
Three levels.

**(i) Mathematical non-identifiability (true for any data)**
- The Jacobian rank is 14 < 15: there is a direction along which a set of coefficients can change simultaneously while the predicted correlation matrix stays **exactly the same**. This direction is mainly "L1→age down, L1→MgH up", compensated by the edges from age and MgH to e, J_φ and [Fe/H]. It is a property of the graph structure, independent of the data.
- **In words**: the data only provide how strongly the 6 quantities correlate pairwise (15 numbers), and the graph has 15 coefficients, two of which (L1→age, L1→MgH) can trade off against each other without changing anything observable. The data cannot tell whether L1 is "more like age" or "more like [Mg/H]" — and the relation of L1 to the birth radius depends precisely on this balance.

**(ii) The optimum sits on the boundary "noise variance = 0", and the likelihood surface has many almost equally good optima**
- Physically, every noise variance must be ≥ 0, which cuts away part of the family of "exactly equal" solutions in (i). But a proper maximum-likelihood fit almost always pushes the noise variance of one observable to 0: in 10 of the 12 TNG50 galaxies without errors it is J_φ (in the other 2 the optimum is interior, with a smallest noise variance of 0.03–0.04); with errors J_φ in 9, [Fe/H] in 2 and e in 1; on g2.79e12 it is [Fe/H] or e (depending on which local optimum is found).
- **In words**: a noise variance of 0 means the model treats that observable as **exactly** a linear combination of its parents, with no random component of its own (e.g. "J_φ is completely determined by L2, age and [Mg/H]"). This is physically impossible and shows that the linear model of this graph can only fit the data by going to an extreme. In factor analysis this is called a Heywood case and usually signals a misspecified model.
- Near this boundary the likelihood surface has many local optima whose NLL differs by only ~1e-4, each giving L1 a very different meaning.

**(iii) Consequence: the meaning of L1 is unstable**
- **TNG50**: solutions statistically indistinguishable from the best give |r(L1, Rbirth)| from about 0 up to 0.4–0.5 in most galaxies (table in Sec. 5.6). Which one is obtained depends on the optimizer's starting values.
- **The authors' g2.79e12**: under the strict criterion only one solution is indistinguishable from the best, but the "best" itself changes with the search effort: with 30 starts the best solution gives |r(L1, Rb)| = 0.25 (|r(L2, Rg)| = 0.76); with 230 starts a solution with an NLL lower by 0.0009 gives 0.62 (0.55). A longer search might change it again. Even on the authors' own data, the correspondence between L1 and the birth radius is not a stable property of the fit.
- **The paper's coefficients are far from the optimum**: their NLL, recomputed with the correct formula, is 1.593 versus about −1.272 at the optimum; with 119,000 stars the likelihood-ratio statistic is about 680,000. The paper's r = 0.64 and 0.98 come from an unconverged solution. Fitting to convergence with the maximum-likelihood method the paper itself describes gives results that change with the search (previous point), and an L2–Rg correlation of only 0.55–0.76, far below 0.98.
- The "freeze the graph and transfer it to other galaxies and the Milky Way" results of paper Secs. 4.2 and 4.3 are affected in the same way.

**Possible remedies** (none attempted yet):
- (a) **report the unreliability itself as the result**, giving the possible range for each galaxy (what this report does);
- (b) add a physically motivated constraint that does not use the truth, removing the free direction;
- (c) add observables (e.g. J_R or more element abundances) so that there are more equations than unknowns, or allow non-linear relations.
- Widening the selection (Sec. 3.12) does **not** solve the problem: the Jacobian rank is still 14/15 in all three selections, with 1–160 indistinguishable solutions (Sec. 5.9).

### 4.4 Other (smaller) issues
| Issue | Description | Handling |
|---|---|---|
| paper vs repository data | the paper says 7–10 kpc, [Fe/H] > −1, low-α disk; the repository data actually have R = 9.4–12.1 kpc, [Fe/H] > −0.5 and Rbirth < 20 kpc; no observational errors are visible in the csv | we use 7–10 kpc as in the paper text, and [Fe/H] > −0.5 and Rbirth < 20 as in the repository |
| eccentricity definition | the paper does not say whether 3D r or cylindrical R is used; the values differ | we use 3D r (to be confirmed; small effect) |
| hard-coded coefficients in the L1/L2 notebooks | with new data the coefficients must be replaced, otherwise the results are wrong | `rlcd_A.py` reads the new coefficients automatically |
| running RLCD inside the repository overwrites the authors' result files | results are written back to `rlcd/results_*/` | move new files to `results/rlcd_g279/` after each run; restore the originals with git |
| `grid_search.py` calls `python` | uses the wrong interpreter if the RLCD environment is not activated | activate the RLCD environment first (`source <RLCD env>/bin/activate`) |
| ActionFinder failed for 530330 | the potential grid extended to 0.01 kpc, where particle noise made the circular angular momentum non-monotonic | inner grid boundary 0.2 kpc (TNG50 softening), automatic retries with 0.5 and 1.0 kpc |
| `Birth_R` in the old `001_orbit` files is a 3D distance | easily misread | the new data store both; the truth uses the 2D radius |
| AGAMA installation failed 3 times | PyPI has only a bootstrap package; Anaconda's OpenMP is too old | build from source with Homebrew's libomp ([`docs/INSTALL_AGAMA.md`](docs/INSTALL_AGAMA.md)) |

---

## 5. Results (by folder)

### 5.1 Data products
| File | Content |
|---|---|
| `data/tng50_extract/data_obtain_all.hdf5` | all particles of the 12 galaxies; stars with age, [Fe/H], [Mg/Fe], birth position and birth radius (server extraction; do not modify) — not distributed, see [`DATA.md`](DATA.md) |
| `data/tng50_derived/orbits_all.hdf5` | J_R, J_z, J_φ, E, Rg, r_peri, r_apo, e, z_max of every star, and rotation curves; rows in exactly the order of the stars in the previous file — not distributed |
| `data/tng50_derived/potentials/<ID>_{halo,disk}.coef` | potential coefficients of each galaxy, reload with `agama.Potential(file=...)` |
| `data/rlcd_input/{nonoise,noise}/<ID>_e.csv` | the 6-column RLCD inputs (distributed as `data/release/*.tar.gz`) |
| `data/rlcd_input/{nonoise,noise}/<ID>_truth.csv` | truth for the same stars (Rbirth = `R_bir`, Rg = `Rg`) and other quantities (distributed) |

### 5.2 Data and orbit checks: `results/check_vs_001_orbit.txt`, `results/orbit_check/`
The new extraction agrees star by star with the earlier results; all orbit quantities of the 12 galaxies pass the consistency checks (Secs. 3.2, 3.4). **The stage-1 results (TNG50 data preparation) are reliable.**

### 5.3 Reproducing the authors' graph: `results/rlcd_g279/`
| File | What it shows |
|---|---|
| `author_result_..._ranktestscale1.png` | the authors' final causal graph (15 edges): L1 → {L2, age, [Fe/H], [Mg/H], J_z}, L2 → {J_φ, e}, plus edges from [Mg/H] and age to other observables |
| `ours_seed1_result_....png`, `compare_author_vs_seed1.png` | our seed-1 graph, **identical** to the authors': installation and data are fine |
| `ours_seed{0,2,3,4}_result_....png`, `compare_author_vs_seed3_diff.png` | graphs from other seeds with the same parameters; the diff figure colours the 4 edges that differ: the children of the second latent change from (J_φ, e) to ([Fe/H], J_z) |
| `ours_noseed_result_alpha<four α>.png` | results for 4 α combinations: α0 = 1e-10 or 1e-15 gives 4 latents |
| `compare_4alpha.txt`, `compare_5seeds.txt` | edge-by-edge comparison reports (latent labels aligned automatically, prior checked) |

### 5.4 Single runs of Approach A (⚠ examples only, not evidence): `results/rlcd_A/g279/`, `results/rlcd_A/nonoise/538905/`
| File | Content |
|---|---|
| `fit.txt` | conversion formula and correlation of the selected solution; the last line lists |r(L1, Rb)| of the 5 best-fitting solutions |
| `nll_repeats.csv` | NLL, RMSE, validity, |r(L1,Rb)|, |r(L2,Rg)| of each of the 50 parameter-learning runs |
| `best_params.png` | the selected coefficients on the authors' graph |
| `latents.csv` | L1, L2 and converted values per star |
| `L1_vs_R_bir.png`, `L2_vs_Rg.png` | converted L1 vs true Rbirth, converted L2 vs true Rg |

g2.79e12: r(L1, Rb) = 0.63, r(L2, Rg) = 0.86; TNG50 538905: r(L1, Rb) = 0.06, r(L2, Rg) = 0.78. But equally good solutions for the same galaxy give |r(L1, Rb)| anywhere from 0.01 to 0.63, so **these numbers on their own mean nothing** (also stated in `results/rlcd_A/README.md`).

### 5.5 Identifiability tests: `results/rlcd_A/identifiability/`
`<tag>` = `g279` (authors' data) or `nonoise_<ID>`, `noise_<ID>` (TNG50).

**Round 1 (loose criterion, `rlcd_identifiability.py`)**, 4 files per data set:
| File | Content |
|---|---|
| `<tag>_summary.txt` | Jacobian rank, non-identifiable direction, conclusions of the multistart and profile (text) |
| `<tag>_multistart.csv` | the solutions from 30 random starts: NLL, validity, |r(L1,Rb)|, |r(L2,Rg)|, the 15 coefficients |
| `<tag>_profile.csv`, `<tag>_profile.png` | profile: best NLL and corresponding |r| with L1→age fixed at 39 values |

**Round 2 (strict criterion, `rlcd_identifiability_ci.py`)**, 2 files per data set:
| File | Content |
|---|---|
| `<tag>_ci.csv` | all valid solutions from 230 starts (the 30 of round 1 + 200 new) after strict convergence: NLL, likelihood-ratio statistic `lr` = 2·N·ΔNLL, the smallest noise variance and its node, |r(L1,Rb)|, |r(L2,Rg)|, the 15 coefficients |
| `<tag>_ci.txt` | text: number of solutions indistinguishable from the best, their |r| ranges, the observable whose noise variance is pushed to 0 at the optimum |

**Jacobian rank**: for all 13 data sets (authors' + 12 TNG50), at 3 random coefficient sets each, the rank is **14/15**. The non-identifiable direction is the same in all data sets (mainly L1→age −0.73, L1→MgH +0.52, the rest on the edges from age and MgH to e, J_φ, [Fe/H]), because it is a property of the graph.

**How to read the profile figures** (e.g. `g279_profile.png`, `nonoise_567607_profile.png`):
- x axis: the fixed L1→age coefficient (−0.95 to 0.95).
- left panel: the best NLL reachable with that value fixed and the other 14 refitted; the dashed line is the theoretical minimum (perfect agreement of model and data).
- right panel: the corresponding |r(L1, true Rbirth)| (blue) and |r(L2, true Rg)| (black). On g2.79e12 the blue curve goes from 0 to about 0.68: near L1→age ≈ 0, L1 resembles [Mg/H] and correlates most strongly with Rbirth ([Fe/H] and [Mg/H] correlate well with Rbirth in the authors' galaxy); near ±0.75, L1 resembles age and is almost unrelated to Rbirth. TNG50 is the other way round (567607: about 0.7 at both ends, about 0 in the middle), because in most TNG50 galaxies **age** correlates much more strongly with Rbirth (0.4–0.7) than [Fe/H] does (0–0.25).
- ⚠ The left panel looks "almost flat", but the axis scale hides real differences: the NLL is per star, and 0.001 times 119,000 stars is statistically significant. The profile shows that "the meaning of L1 is very sensitive to the coefficients", but cannot by itself prove "equally good fit". The strict judgement is round 2.

### 5.6 Summary of the 12 galaxies (main result): `summary_ci_{nonoise,noise}.*`, `linear_ceiling.*`
- "Indistinguishable solutions": after strict convergence, solutions with 2·N·ΔNLL < 3.84 relative to the best solution found (95% confidence).
- "Strict range": the minimum to maximum |r| over these solutions. Starts are finite (230 per galaxy): if a better optimum exists that was not found, the "best" and this range can change (g2.79e12 is an example), so it is the range "found so far".
- "Loose range": the range under the round-1 criterion "NLL within 0.01" (`summary_{nonoise,noise}.*`), listed for comparison.
- "Noise ≈ 0": the observable whose noise variance is pushed to about 0 at the optimum (Sec. 4.3 (ii)).
- "Linear ceiling": correlation of a linear regression of the truth on the 6 inputs (Sec. 3.10); "age alone / [Fe/H] alone": |r| of a single input with Rbirth.

**No observational errors (nonoise)**:

| Galaxy | Stars | Indist. / valid | Strict \|r(L1, Rbirth)\| | Loose \|r(L1, Rbirth)\| | Strict \|r(L2, Rg)\| | Noise ≈ 0 | Age alone | [Fe/H] alone | Linear ceiling |
|---|---|---|---|---|---|---|---|---|---|
| g2.79e12 (authors) | 119,236 | 1 / 223 | **0.62** | 0.00 – 0.63 | 0.55 | e | 0.26 | 0.70 | 0.85 |
| 487742 | 13,190 | 36 / 225 | **0.00 – 0.42** | 0.00 – 0.49 | 0.60 – 0.93 | J_φ | 0.54 | 0.15 | 0.73 |
| 530330 | 24,154 | 59 / 224 | **0.21 – 0.50** | 0.03 – 0.59 | 0.78 – 0.94 | J_φ | 0.68 | 0.09 | 0.81 |
| 538905 | 34,970 | 47 / 223 | **0.09 – 0.47** | 0.07 – 0.45 | 0.75 – 0.96 | J_φ | 0.55 | 0.07 | 0.70 |
| 545437 | 16,222 | 35 / 226 | **0.09 – 0.49** | 0.02 – 0.46 | 0.61 – 0.93 | J_φ | 0.53 | 0.00 | 0.66 |
| 546870 | 11,484 | 118 / 226 | **0.34 – 0.40** | 0.29 – 0.41 | 0.92 – 0.93 | [Fe/H] (0.039) | 0.40 | 0.20 | 0.53 |
| 547293 | 32,525 | 83 / 225 | **0.00 – 0.25** | 0.02 – 0.36 | 0.81 – 0.92 | J_φ | 0.53 | 0.25 | 0.81 |
| 552581 | 13,573 | 36 / 221 | **0.00 – 0.41** | 0.01 – 0.50 | 0.63 – 0.72 | J_φ | 0.23 | 0.46 | 0.67 |
| 557721 | 25,326 | 75 / 221 | **0.12 – 0.41** | 0.00 – 0.44 | 0.82 – 0.95 | J_φ | 0.53 | 0.11 | 0.71 |
| 560421 | 19,696 | 91 / 223 | **0.00 – 0.56** | 0.02 – 0.55 | 0.77 – 0.97 | [Fe/H] (0.026) | 0.58 | 0.25 | 0.79 |
| 562029 | 7,299 | 9 / 224 | **0.08 – 0.50** | 0.02 – 0.69 | 0.52 – 0.86 | J_φ | 0.69 | 0.08 | 0.85 |
| 567607 | 12,496 | 20 / 225 | **0.05 – 0.53** | 0.03 – 0.71 | 0.66 – 0.93 | J_φ | 0.68 | 0.14 | 0.82 |
| 569251 | 29,035 | 63 / 223 | **0.00 – 0.17** | 0.00 – 0.21 | 0.82 – 0.94 | J_φ | 0.28 | 0.22 | 0.69 |

**With observational errors (noise)**:

| Galaxy | Stars | Indist. / valid | Strict \|r(L1, Rbirth)\| | Loose \|r(L1, Rbirth)\| | Strict \|r(L2, Rg)\| | Noise ≈ 0 | Age alone | [Fe/H] alone | Linear ceiling |
|---|---|---|---|---|---|---|---|---|---|
| g2.79e12 (authors) | 119,236 | 1 / 223 | **0.62** | 0.00 – 0.63 | 0.55 | e | 0.26 | 0.70 | 0.85 |
| 487742 | 13,190 | 22 / 225 | **0.02 – 0.34** | 0.00 – 0.52 | 0.81 – 0.87 | J_φ | 0.54 | 0.15 | 0.71 |
| 530330 | 24,154 | 15 / 223 | **0.09 – 0.50** | 0.01 – 0.62 | 0.85 – 0.92 | J_φ | 0.66 | 0.09 | 0.76 |
| 538905 | 34,970 | 23 / 224 | **0.10 – 0.38** | 0.00 – 0.45 | 0.77 – 0.94 | J_φ | 0.54 | 0.07 | 0.67 |
| 545437 | 16,222 | 24 / 225 | **0.07 – 0.46** | 0.01 – 0.46 | 0.71 – 0.90 | J_φ | 0.52 | 0.00 | 0.64 |
| 546870 | 11,484 | 134 / 224 | **0.27 – 0.33** | 0.21 – 0.43 | 0.95 – 0.95 | J_φ | 0.36 | 0.20 | 0.45 |
| 547293 | 32,525 | 36 / 224 | **0.01 – 0.38** | 0.00 – 0.45 | 0.80 – 0.94 | J_φ | 0.52 | 0.25 | 0.79 |
| 552581 | 13,573 | 35 / 224 | **0.00 – 0.36** | 0.00 – 0.45 | 0.64 – 0.71 | J_φ | 0.22 | 0.45 | 0.63 |
| 557721 | 25,326 | 34 / 226 | **0.07 – 0.41** | 0.02 – 0.41 | 0.64 – 0.72 | e | 0.49 | 0.11 | 0.62 |
| 560421 | 19,696 | 35 / 225 | **0.01 – 0.58** | 0.00 – 0.58 | 0.83 – 0.98 | J_φ | 0.57 | 0.25 | 0.77 |
| 562029 | 7,299 | 19 / 224 | **0.01 – 0.67** | 0.02 – 0.69 | 0.67 – 0.83 | J_φ | 0.67 | 0.08 | 0.83 |
| 567607 | 12,496 | 17 / 226 | **0.09 – 0.71** | 0.00 – 0.71 | 0.81 – 0.89 | [Fe/H] | 0.67 | 0.14 | 0.78 |
| 569251 | 29,035 | 59 / 225 | **0.00 – 0.24** | 0.00 – 0.26 | 0.81 – 0.95 | [Fe/H] | 0.28 | 0.22 | 0.63 |

![summary_ci_nonoise](results/rlcd_A/identifiability/summary_ci_nonoise.png)

`summary_ci_nonoise.png` (likewise `summary_ci_noise.png`): one row per galaxy. Thick blue line = strict range, thin grey line = loose range, black triangle = age alone, cross = [Fe/H] alone, open circle = linear ceiling, grey dashed line = the paper's 0.64.

**How to read these results**:
1. **The meaning of L1 is unstable**: in most TNG50 galaxies, indistinguishable solutions give |r(L1, Rbirth)| from about 0 to 0.4–0.5; only a few galaxies (e.g. 546870) have a narrow range. Which one is obtained depends on the optimizer's starting values.
2. **Even the most favourable solution is essentially no better than a single input**: in `summary_ci_nonoise.png` the right end of every blue line lies left of the black triangle or the cross (546870 roughly coincides). Without errors, in none of the 12 TNG50 galaxies does L1 ever correlate with the birth radius more strongly than the better of "age alone" and "[Fe/H] alone". With errors there are 2 exceptions, only at the favourable end and only marginally: 560421 (0.58 vs 0.57 for age) and 567607 (0.71 vs 0.67 for age), and in these two galaxies L1 can also be as low as 0.01 and 0.09. On the authors' own g2.79e12, [Fe/H] alone (0.70) also beats the paper's L1 (0.64) and our best solution (0.62).
3. **Far from the linear ceiling**: the birth-radius information linearly available in the 6 quantities (open circles, 0.53–0.85) is clearly larger than what L1 captures. The structure of L1 (extracting one common factor first) loses information.
4. **L2 correlates well with Rg** (mostly 0.6–0.95), but Rg is essentially J_φ in other units (|r(J_φ, Rg)| = 0.98–1.00 except in 552581), and in most galaxies the optimum pushes the noise variance of J_φ to 0, i.e. the model treats J_φ as an exact linear combination of L2 (plus age and [Mg/H]). So "L2 resembles Rg" is not a discovery. L2 is worse in 552581 (about 0.63–0.71) because about 5% of its selected stars counter-rotate (J_φ < 0), while Rg is defined from |J_φ|, so the two are no longer linearly related.
5. **Adding observational errors (noise) does not change the conclusions.**
6. **Loose vs strict**: both criteria give similar ranges (the strict one used 200 more starts and is sometimes even wider); points 1 and 2 hold under both.

### 5.7 How these result files were produced
```bash
cd <repo root>
PY=python3      # python of the analysis environment
# round 1 (loose)
$PY scripts/local/rlcd_identifiability.py --version g279 --nstart 30
for v in nonoise noise; do
  for g in 487742 530330 538905 545437 546870 547293 552581 557721 560421 562029 567607 569251; do
    $PY scripts/local/rlcd_identifiability.py --version $v --ids $g --nstart 30
  done
  $PY scripts/local/rlcd_identifiability_summary.py --version $v
done
# linear ceiling
$PY scripts/local/linear_ceiling.py
# round 2 (strict)
$PY scripts/local/rlcd_identifiability_ci.py --version g279 --extra 200
$PY scripts/local/rlcd_identifiability_ci.py --version nonoise --extra 200
$PY scripts/local/rlcd_identifiability_ci.py --version noise --extra 200
```
(Reference run time on a 2-core Linux machine: about 1 hour per version for round 1 and 1–2 hours per version for round 2.)

### 5.8 Figures: `results/TNG/causal_graph/`, `results/TNG/truth_comparison/`
- `tng_causal_graph.png` (`_noise.png`): each galaxy's maximum-likelihood coefficients on the authors' graph. Observables marked "ω≈0" are those whose noise variance is pushed to 0 (Sec. 4.3 (ii)); n_eq in each panel title is the number of indistinguishable solutions — only the one with the smallest NLL is drawn.
- `tng_L1_Rbirth_L2_Rg.png` (`_noise.png`): f(L1) vs true Rbirth and f(L2) vs true Rg. **The L1 panels are almost all narrow vertical bands**, very different from paper Fig. 3 (which spans 0–20 kpc). The L1 computation was checked (it differs by < 10⁻⁷ from the per-star Nelder–Mead solution of the notebook); the narrow bands are an inevitable consequence of the weak correlation:

  the least-squares slope of the conversion f = α·L1 + β is α = r · σ_R / σ_L, so

  σ_f = |r| · σ_R

  Meaning: the spread of the converted f(L1) equals "the correlation of L1 with the truth" times "the spread of the truth itself". The weaker the correlation, the more f(L1) is squeezed towards the mean.
  - σ_f: standard deviation of f(L1) (how widely the points spread along the x axis); σ_R: standard deviation of the true birth radius; σ_L: standard deviation of L1
  - r: Pearson correlation of L1 with the true birth radius (Sec. 8.3.7); α, β: slope and intercept of the conversion line

  Example: in the authors' galaxy r = 0.64 and σ_R = 3.67 kpc, so σ_f ≈ 2.35 kpc and the points spread out; for the TNG50 7–10 kpc selection σ_R is only about 2.2 kpc and |r| is between 0 and 0.56 (Sec. 5.6), so σ_f is below 1.2 kpc and the points form a vertical band. This is not a plotting issue: L1 itself is only weakly related to the birth radius.

### 5.9 Wider selections: `star_condition/<selection>/output_data/summary_ci_{nonoise,noise}.*`
The `summary_ci_<version>.txt/.csv/.png` files of each selection have the same format as Sec. 5.6, with two more columns: `r1_best` (|r(L1, Rbirth)| of the maximum-likelihood solution) and `r_JPhi_Rbir` (|r| of J_φ alone with the birth radius). The Jacobian rank is 14/15 in all three selections.

**No errors (nonoise)**: each cell gives "|r(L1, Rbirth)| of the maximum-likelihood solution (range over indistinguishable solutions)", next to the single input most correlated with the birth radius and its |r|; **bold** = the maximum-likelihood solution beats the best single input.

| Galaxy | \|z\|<30: L1 | best single | \|z\|<3: L1 | best single | \|z\|<0.3: L1 | best single |
|---|---|---|---|---|---|---|
| 487742 | 0.79 (0.75–0.80) | 0.86 J_φ | 0.78 (0.74–0.80) | 0.87 J_φ | 0.79 (0.22–0.79) | 0.89 J_φ |
| 530330 | **0.80 (0.69–0.81)** | 0.78 J_φ | **0.80 (0.66–0.81)** | 0.79 J_φ | 0.81 (0.70–0.82) | 0.84 J_φ |
| 538905 | **0.80 (0.67–0.80)** | 0.79 J_φ | **0.81 (0.68–0.81)** | 0.80 J_φ | 0.84 (0.72–0.88) | 0.86 J_φ |
| 545437 | 0.78 (0.66–0.83) | 0.86 J_φ | 0.83 (0.66–0.83) | 0.87 J_φ | 0.80 (0.73–0.81) | 0.87 J_φ |
| 546870 | **0.71 (0.08–0.71)** | 0.64 [Mg/H] | 0.09 (0.09–0.09) | 0.63 [Mg/H] | 0.17 (0.06–0.63) | 0.63 [Fe/H] |
| 547293 | 0.85 (0.85–0.85) | 0.85 J_φ | 0.83 (0.76–0.83) | 0.86 J_φ | 0.82 (0.71–0.86) | 0.87 J_φ |
| 552581 | 0.49 (0.49–0.50) | 0.63 [Mg/H] | 0.54 (0.24–0.54) | 0.61 [Mg/H] | 0.45 (0.45–0.65) | 0.71 J_φ |
| 557721 | 0.82 (0.49–0.90) | 0.90 J_φ | 0.58 (0.48–0.90) | 0.90 J_φ | 0.89 (0.49–0.90) | 0.92 J_φ |
| 560421 | 0.69 (0.13–0.76) | 0.83 J_φ | 0.69 (0.13–0.76) | 0.83 J_φ | 0.79 (0.19–0.80) | 0.87 J_φ |
| 562029 | 0.77 (0.69–0.79) | 0.82 J_φ | 0.80 (0.69–0.80) | 0.82 J_φ | 0.81 (0.78–0.81) | 0.83 J_φ |
| 567607 | 0.86 (0.75–0.86) | 0.87 J_φ | 0.86 (0.81–0.87) | 0.87 J_φ | 0.89 (0.80–0.89) | 0.90 J_φ |
| 569251 | 0.46 (0.03–0.49) | 0.85 J_φ | 0.05 (0.01–0.52) | 0.85 J_φ | 0.28 (0.28–0.78) | 0.87 J_φ |

Other numbers (all three selections): standard deviation of the birth radius 2.1–5.1 kpc; linear ceiling of the 6 inputs 0.73–0.95; maximum-likelihood |r(L2, Rg)| 0.21–0.93; the node with noise variance ≈ 0 at the optimum is mostly [Fe/H] or J_φ, occasionally L2.

**With errors (noise)**: same pattern but less stable. The maximum-likelihood solution beats the best single input in only 4 cases (530330 and 546870 for |z|<30; 538905 for |z|<3; 538905 for |z|<0.3), by at most 0.06; in several galaxies the maximum-likelihood |r| drops very low (e.g. for |z|<0.3: 0.10 for 560421, 0.21 for 557721, 0.29 for 545437), while indistinguishable solutions in the same galaxy reach 0.77–0.90. Full tables in each selection's `summary_ci_noise.txt`.

**How to read these results**:
1. **L1 seems to "improve", but only because the selection lets in the current radius**: across the whole disk, where a star is now and where it was born are strongly correlated (radial migration only adds scatter on top of that trend), so J_φ alone (angular momentum ≈ current radius × circular speed) already reaches 0.8–0.9. The paper uses only stars at 7–10 kpc precisely to remove this "trivial" information and test whether the method recovers the birth radius after migration.
2. **L1 still essentially never beats the best single input**: without errors only 5 of 36 maximum-likelihood solutions do, by at most 0.07; with |z| < 0.3 (the paper's height cut) none does.
3. **The coefficients are still not unique, and the maximum-likelihood solution jumps with the selection**: for 546870 the maximum-likelihood |r| is 0.71, 0.09 and 0.17 for |z| < 30, 3 and 0.3; for 569251 it is 0.46, 0.05 and 0.28; for 557721 0.82, 0.58 and 0.89. Only the |z| limit changes (|z| < 3 has just 3% fewer stars than |z| < 30), yet the meaning of L1 changes completely.
4. **So wider selections do not change the conclusions of Sec. 6**; they rather show that how well L1 correlates with the birth radius is set mainly by how much birth-radius information the inputs already contain (i.e. by the selection), not by the method finding the birth radius.

---

## 6. Conclusions

### 6.1 Does RLCD work on TNG50?
**Not as done in the paper (frozen authors' graph, re-estimated coefficients).** In order of importance:
1. **Even in the luckiest case L1 is essentially no better than a single observable**: without errors, in none of the 12 TNG50 galaxies does the most favourable solution correlate with the birth radius better than "age alone" or "[Fe/H] alone"; with errors only 2 galaxies are marginally higher at the favourable end (0.58 vs 0.57, 0.71 vs 0.67). A method for "finding the birth radius" that is no better than using age directly has no practical value.
2. **The meaning of L1 is not determined by the data** (Problem 3): the authors' graph has one degree of freedom invisible to the data, and the likelihood surface has many almost equally good optima, so L1 can range between "age-like" and "[Mg/H]-like". Without the truth one cannot tell which one was obtained.
3. **The linear model does not fit the data**: in 22 of the 24 TNG50 data sets the optimum needs the noise variance of J_φ (or [Fe/H], e) to vanish, so the linear model of this graph does not describe TNG50 well.
4. **The graph itself is unstable** (Problem 1): with the same data and parameters the authors' graph appears in only 7 of 30 seeds, and it is not the most frequent graph.
5. **TNG50 physics differs from the authors' simulation**: for TNG50 disk stars at 7–10 kpc [Fe/H] is almost uncorrelated with the birth radius (|r| = 0–0.25 in 11 of 12 galaxies) and the birth-radius information is mostly in the age, whereas the authors' graph has chemistry as the main children of L1. A causal graph learned on one simulated galaxy does not necessarily transfer to others.
6. **Wider selections do not change the conclusion** (Sec. 5.9): with the whole disk |r(L1, Rbirth)| rises to about 0.8, but J_φ alone (≈ current radius) already gives 0.8–0.9 and L1 still essentially never beats it; the coefficients remain non-unique and the maximum-likelihood solution jumps with the |z| limit. How "good" L1 looks is set by how much birth-radius information the selection puts into the inputs, not by the method finding the birth radius.

**Positive results**: the data-preparation pipeline (birth radii, actions, Rg, eccentricities) is complete and validated and can be used to test any other method; L2 consistently tracks J_φ/Rg (but trivially).

### 6.2 Can it be applied to Milky Way observations?
**Not at present:**
1. **Problem 3 exists for Milky Way data too**: non-identifiability is a property of the graph structure, regardless of whether the data are simulated or observed. Fitting the authors' graph to Milky Way data, L1 will likewise move between different solutions depending on the starting values.
2. **The Milky Way has no ground truth, so one cannot tell which L1 one has.** The paper "validates" L1 with birth radii inferred by the method of Lu et al. (2024), but those are themselves computed from [Fe/H] and age, while L1 is also a linear combination of [Fe/H], [Mg/H], age, ...; the correlation is to a large extent a comparison of a quantity with itself, not an independent validation.
3. **Observations add further difficulties**: selection effects (APOKASC-3 is an asteroseismic sample), large age uncertainties, and actions that depend on an assumed Milky Way potential. None of these was tested in simulations.
4. **Ask the authors for data or redo it?** To reproduce paper Fig. 6 one can simply use `obs_e.csv` and `obs_with_sim_dag.csv`, already in the authors' repository; but until the problems above are solved a reproduction would not show that the method works.
5. **The paper's Milky Way results are affected by Problems 2 and 3**: the repository's coefficient-learning command for the Milky Way is the same as for g2.79e12 (SGD, 600 steps) and may also be unconverged; the frozen graph is the same and equally non-identifiable.

**Minimum requirement for it to work**: first find, on TNG50, a causal graph that is identifiable, free of Heywood cases, and stable across seeds and galaxies, and show that its L1 predicts the birth radius better than "age alone" or "[Fe/H] alone"; only then is it worth applying to the Milky Way.

---

## 7. Suggested next steps
1. **Make the model identifiable** (Problem 3, options (b), (c)): e.g. add J_R or more element abundances so that there are more correlations than edges, or propose a physically motivated constraint. After any change, first check that the Jacobian has full rank, then look at L1 vs Rbirth.
2. **Always compare with baselines**: any birth-radius estimate must be compared with "age alone", "[Fe/H] alone", "J_φ alone" (especially for wide selections, Sec. 5.9) and the "linear ceiling of the 6 quantities"; only beating the single-variable baselines is meaningful.
3. **Standard practice for coefficient estimation**: always impose "noise variance ≥ 0", use a converging optimizer with hundreds of starts, and report all solutions statistically indistinguishable from the best, not just one.
4. **Eccentricity definition**: decide between 3D r and cylindrical R (small effect, but should be consistent).
5. **Contact the authors** (optional): ask politely about the randomness of the graph, which coefficient set was used and why, and identifiability.

---

## 8. Methods

This part is for readers new to RLCD: first what the paper does (8.1), then how the code is organized (8.2), the principles of the algorithm (8.3), and the tools we used to test it (8.4).

### 8.1 What the original paper does
**Title**: Causal Discovery of Latent Variables in Galactic Archaeology (Jin, Lu, Ting, Zheng, Buck).

**Core idea**: ordinary machine learning finds "correlations"; causal discovery tries to find "what determines what". A star's birth radius is unobservable, but it affects several observables at once (stars born in the inner disk are more metal-rich, older, ...). If the correlations among several observables can all be explained by "one common hidden cause behind them", the existence of that hidden cause (a latent variable) can be inferred and its value estimated for each star.

**Steps of the paper**:
1. **Data (Sec. 2)**: NIHAO-UHD galaxies g2.79e12 (main) and g8.26e11; disk stars currently at 7–10 kpc with [Fe/H] > −1, low-α disk only ([Mg/Fe] < −0.01); observational errors added to age, [Fe/H], [Mg/Fe]. The 6 inputs: [Fe/H], [Mg/H], age, J_z, J_φ, e. Orbit quantities from AGAMA in an axisymmetric multipole potential. Milky Way data: APOKASC-3 (APOGEE abundances + Kepler asteroseismic ages + Gaia kinematics).
2. **Structure learning (Sec. 3.1)**: RLCD finds the causal graph from the covariance matrix of the 6 quantities, including how many latents there are and what they connect to.
3. **Physical prior (Sec. 3.2)**: a grid search over the significance levels α gives one candidate graph per α combination; graphs where "age has an observed parent" are removed (present abundances and orbits cannot determine age), and the paper states that the remaining ones converge to a unique graph.
4. **Parameter learning (Sec. 3.3)**: with the structure fixed, maximum likelihood gives each edge's coefficient (the numbers on the edges).
5. **Latent values (Sec. 3.3)**: for each star find L1, L2 minimizing the difference between predicted and observed observables; L2 first, then L1.
6. **Results (Sec. 4)**: on g2.79e12, 2 latents; L1 correlates with the birth radius at r = 0.64, L2 tightly with the guiding radius; freezing this graph and re-estimating only the coefficients on g8.26e11 and the Milky Way, they conclude that it "transfers".
7. **Discussion (Sec. 5)**: L1 is recovered less well than L2, which the authors attribute to intrinsic scatter in birth environments and to the "linear model" approximation.

**Not discussed in the paper**: the randomness of structure learning, the multiple optima of parameter learning and the selection criterion, whether the coefficients are identifiable once the graph is frozen, and comparison with "a single observable" baselines.

### 8.2 How the code is organized
Authors' repository `Causal-Galactic-Archaeology/`; all commands run inside `rlcd/`. A full run calls, in order:

```
run_rlcd.py --s <data set> --alpha0 .. --alpha3 ..
 ├─ scenarios.py: read the csv, prefix columns with X_, standardize (subtract mean, divide by std)
 ├─ utils/Chi2RankTest.py: prepare the rank tests (covariance of the full data)
 └─ StructureLearning/RLCD/RLCD_alg.py: RLCD()
      ├─ stage 1: GES (causal-learn) learns an ordinary causal graph among the observables, then groups variables by connectivity (getPartition)
      ├─ stage 2: for each group, LatentGroups + findClusters(): for k = 1, 2, 3 search for "rank-deficient" variable sets and introduce latents
      │     └─ structuralRankTest() → Chi2RankTest.test(): is the rank of a covariance block ≤ k?
      └─ output: adjacency matrix of observables + latents → results_<data set>/result_alpha..._ranktestscale1.dot/.png

ParameterLearning/linear_parameter_identification.py --dot_path <graph>
 ├─ get_F_bin_and_varnames_from_dotgraph(): which edges exist (F_bin, a 0/1 matrix) from the .dot
 ├─ parameter_by_trek(): random start → model covariance via the trek rule → maximum likelihood → 600 SGD iterations
 └─ add_edge_coefficients_to_graph(): writes the coefficients back into the graph, NLL in the file name

L1L2calculation/GetL1_L2_279sim.ipynb
 └─ coefficients hard-coded → per star, Nelder–Mead minimization of χ² for L2, L1 → np.polyfit linear conversion to Rbirth, Rg
```

Parameters:
- `--alpha0 … --alpha3`: significance levels for testing "rank ≤ 0, ≤ 1, ≤ 2, ≤ 3" (8.3.3).
- `--stage1_ges_sparsity 2`: penalty factor on the number of edges in the GES score; larger means sparser graphs.
- `--lr 0.005 --epochs 600 --lr_decay plateau`: learning rate, number of iterations, learning-rate schedule of parameter learning.

### 8.3 Principles of the algorithm

#### 8.3.1 Linear structural causal model
RLCD assumes each variable is a linear function of its "parents" plus its own random disturbance (paper Eq. 1):

V_i = Σ_{j ∈ Pa(V_i)} a_ji · V_j + ε_i

Meaning: the value of V_i = the sum of each parent's value times the corresponding edge coefficient, plus a noise term independent of all other variables.
- V_i: any variable (observable or latent)
- Pa(V_i): the parents of V_i (variables with an arrow pointing to V_i)
- a_ji: the coefficient of the edge from V_j to V_i (the number next to the arrow), i.e. "how much V_i changes on average when V_j changes by one unit"
- ε_i: the noise of V_i, mean 0, independent of everything else

Example: in the authors' graph the parents of [Fe/H] are L1, [Mg/H] and age, so [Fe/H] = a·L1 + b·[Mg/H] + c·age + ε_FeH.

The "causal graph" is a picture of this set of equations: one node per variable and one arrow per non-zero coefficient. Since all variables are standardized (variance = 1), the coefficients are roughly "standardized strengths of influence" between −1 and 1.

#### 8.3.2 Why the "rank" of the covariance matrix reveals hidden variables
**Simplest example**: one latent L with 4 observed children X1–X4, X_i = a_i·L + ε_i. Then the covariance of any two X's is cov(X_i, X_j) = a_i·a_j (with unit variance for L). Arrange the covariances between {X1, X2} and {X3, X4} in a 2×2 matrix:

| | X3 | X4 |
|---|---|---|
| X1 | a1·a3 | a1·a4 |
| X2 | a2·a3 | a2·a4 |

The second row is the first row times a2/a1, so the rows are proportional and **the rank is 1** (only one independent row). If the 4 quantities were directly connected to each other without a common cause, this matrix would generally have full rank (rank 2).
**In words**: the rank is "the minimum number of independent channels needed to explain the correlations between the two groups". Rank 1 means all association between {X1, X2} and {X3, X4} goes through a single channel — most likely a common hidden cause.

**General rule (trek separation, paper Sec. 3.1)**:
- **trek**: a pair of directed paths starting from a common source and ending at two variables. Two variables are correlated because there are treks between them (e.g. L → X1 and L → X3 form a trek with source L).
- For two groups of variables A and B, the rank of the covariance block Σ_A,B equals the smallest number of nodes needed to "cut" all treks between A and B.
- So if the rank of Σ_A,B is smaller than the number of variables in A and B, the association between A and B is "bottlenecked" by a few nodes; if those nodes are not among the observables, they are latents.

The whole idea of RLCD: **systematically test the ranks of covariance blocks between many groupings of variables, and infer the number and position of latents from the pattern of rank deficiencies.**

#### 8.3.3 The rank test: "is the rank ≤ k?" (`utils/Chi2RankTest.py`)
With noisy data the covariance matrix always has full rank, so one can only test statistically whether "the rank can be considered ≤ k". The method is **canonical correlation analysis** (CCA): find a linear combination of group A and one of group B with maximal correlation, the first canonical correlation ρ₁; then the next one in directions uncorrelated with the previous, ρ₂; ... giving ρ₁ ≥ ρ₂ ≥ … If the rank really is ≤ k, all canonical correlations after the k-th should be close to 0.

Test statistic (Bartlett approximation, as in the code):

T = [N − k − (p+q+1)/2 + Σ_{i≤k}(1/ρ_i² − 1)] × Σ_{i>k} [−ln(1 − ρ_i²)]

Meaning: how large the canonical correlations beyond the k-th are in total, scaled up by the sample size. The closer they are to 0, the smaller T.
- N: number of stars (sample size)
- k: the rank being tested (null hypothesis "rank ≤ k")
- p, q: numbers of variables in groups A and B
- ρ_i: the i-th canonical correlation
- ln: natural logarithm

If the rank really is ≤ k, T approximately follows a χ² distribution with (p−k)(q−k) degrees of freedom. Given a significance level α_k, if T is below the 1−α_k quantile of that distribution, **"rank ≤ k" is accepted** ("fail to reject" in the code).
- **The smaller α**, the higher the threshold and the easier a rank deficiency is accepted, so the more latents are introduced. This is why α0 = 1e-10 or 1e-15 gave 4 latents.
- **Why such extreme α (1e-5 to 1e-15)**: with 100,000 stars even tiny correlations are statistically "significant", and at the usual 0.05 almost no rank deficiency would be detected. α here effectively sets "how small a correlation can be ignored" and has to be tuned, which is why the paper uses a grid search.

#### 8.3.4 The three stages of RLCD (`StructureLearning/RLCD/RLCD_alg.py`)
1. **Stage 1: a skeleton among the observables (GES).** GES (Greedy Equivalence Search) is a score-based method: starting from an empty graph it repeatedly adds the edge that most improves the BIC score, then removes edges one at a time, ending with a causal graph among the observables only. BIC score = goodness of fit − complexity penalty (more edges and more samples → larger penalty; `sparsity = 2` doubles the penalty for sparser graphs). Tightly connected variables are then grouped and latents are searched for in each group separately. On g2.79e12 all 6 variables form one group.
2. **Stage 2: introduce latents level by level (`findClusters`).** Starting from k = 1:
   - take k+1 variables of the current set as group A and the rest as the reference group B;
   - test whether the rank of Σ_A,B is ≤ k. If so, all associations between these k+1 variables and the rest pass through k common latents, so k new latents are added as parents of A (such a group is an "atomic cover", the smallest identifiable unit);
   - replace A by the new latents and return to k = 1 (new latents can themselves become children of other latents, which is how hierarchies like L1 → L2 arise);
   - if no rank deficiency is found for some k, try k+1; stop when k exceeds 3.
   This stage also checks colliders, whether observables can be parents of latents, and orients edges.
3. **Stage 3: output.** The results of all groups are merged into one graph with observables and latents (.dot). The paper then filters candidate graphs with the prior "age has no observed parent".

**Where the randomness comes from**: the order in which stage 2 tries the candidate sets A follows Python set iteration, which depends on the random string hashing at start-up; which rank deficiency is found first affects later merges. Without a fixed PYTHONHASHSEED each run can give a different graph (Problem 1).

**What RLCD guarantees**: under "linear relations", "every latent has enough observed children" and "rank deficiencies come only from the graph structure, not from accidental cancellations of parameters" (rank faithfulness), RLCD finds the correct graph for large enough samples (up to the Markov equivalence class, i.e. some edge directions cannot be distinguished). Whether TNG50 and the Milky Way satisfy these conditions has not been tested.

#### 8.3.5 Parameter learning: coefficients for every edge (`parameter_by_trek`)
Once the graph is fixed, every set of edge coefficients implies a "model correlation matrix" Σ(θ). Parameter learning finds the coefficients for which Σ(θ) is as close as possible to the data correlation matrix S.

**Step 1: the model correlation matrix from the coefficients (Wright's path rule / trek rule).** With all variances equal to 1:

corr(V_i, V_j) = Σ_{all treks between i and j} (product of all edge coefficients on the trek)

Meaning: the correlation of two variables is the sum of the contributions of every "common-source path" connecting them, each contribution being the product of the edge coefficients along it.
- Example: in the authors' graph age and [Fe/H] are connected by the trek "age → [Fe/H]" (contribution a_age→FeH), "age ← L1 → [Fe/H]" (contribution a_L1→age · a_L1→FeH), and longer treks through [Mg/H]; their sum is the model's corr(age, [Fe/H]).

Each variable's **noise variance** follows:

ω_k = 1 − Σ_p Σ_q a_pk · a_qk · corr(V_p, V_q)

Meaning: the total variance of variable k (= 1) minus the part explained by its parents is its own noise variance.
- ω_k: noise variance of variable k, **must be ≥ 0**
- p, q: run over all parents of k
- a_pk: coefficient of the edge from parent p to k

The noise variance of latents is fixed to 1 (paper Sec. 3.3), removing the arbitrary "unit" of a latent.

**Step 2: maximum likelihood.** Assuming the data are multivariate Gaussian, the fit quality is the negative log-likelihood (NLL):

NLL(θ) = ½ ln det Σ(θ) + ½ tr[Σ(θ)⁻¹ S]

Meaning: a measure of how different the model correlation matrix Σ(θ) is from the data correlation matrix S. When Σ(θ) equals S exactly, the NLL reaches its minimum ½ ln det S + p/2.
- θ: all edge coefficients
- Σ(θ): the 6×6 model correlation matrix
- S: the 6×6 data correlation matrix
- det: determinant; tr: trace (sum of diagonal elements); Σ⁻¹: matrix inverse
- p = 6: number of observables, so the theoretical minimum is ½ ln det S + 3

**The repository's implementation**: random starts (`kaiming_normal` followed by tanh), 600 SGD iterations (momentum 0.9), learning rate halved on plateaus; a penalty when a coefficient exceeds 0.99 in absolute value. The determinant is taken in **absolute value** and ω_k ≥ 0 is not checked (Problem 2).

**Our implementation** (`rlcd_identifiability.py`, `rlcd_identifiability_ci.py`; `rlcd_A.py` v2 still uses the repository's SGD with only check (1) added): identical model (layer-by-layer recursion in topological order, equivalent to the trek rule, verified), but (1) ω_k is computed explicitly and required to be ≥ 0 (round 1 used a soft ≥ 0.01), Σ positive definite, otherwise the solution is invalid; (2) L-BFGS-B (a quasi-Newton method that really converges to a local optimum) replaces 600 SGD steps; (3) many random starts are compared.

#### 8.3.6 Latent values for each star
With the coefficients fixed, L2 and L1 are computed for every star (lower level first). For every child c of a latent L:
- first remove the contributions of c's other parents: r_c = c − Σ_{p ≠ L} a_pc · p
  Meaning: from the observed value of c subtract what its known parents (e.g. age, [Mg/H]) contribute; what remains should come mainly from L.
  - c: the observed value of the child (standardized); p: a parent of c other than L; a_pc: coefficient from p to c
- then L = Σ_c a_Lc · r_c / Σ_c a_Lc²
  Meaning: the L value that minimizes the sum of squared differences between "a_Lc · L" and each r_c (least squares).
  - a_Lc: coefficient from L to child c; the sum runs over all children of L
- Order: the children of L2 are observables only (J_φ, e), so L2 comes first; the children of L1 include L2, whose freshly computed value is used as an "observed child" of L1.

This is mathematically identical to the per-star Nelder–Mead minimization of χ² in the authors' notebook (checked with random data, identical to 6 decimal places), but much faster.
**Important consequence**: r_c is a linear combination of observables and L is a linear combination of the r_c, so **L1 and L2 are linear combinations of the 6 inputs**. Therefore |r(L1, Rbirth)| can never exceed the correlation of a linear regression of Rbirth on the 6 quantities, for any coefficients (the linear ceiling of Sec. 3.10).

#### 8.3.7 Conversion to physical units
L1 and L2 are dimensionless (standardized); a straight-line fit to the truth converts them to kpc: Rbirth ≈ α·L1 + β, Rg ≈ α'·L2 + β' (α, α′ slopes, β, β′ intercepts, by least squares). This uses the truth and only changes units; it does not affect the correlation. The strength of the relation is the Pearson correlation:

r = Σ(x − x̄)(y − ȳ) / √[Σ(x − x̄)² · Σ(y − ȳ)²]

Meaning: how consistently x and y deviate from their means in the same direction. |r| = 1 means a perfect straight line, 0 means no linear relation.
- x: L1 (or L2) of each star; y: true Rbirth (or Rg); x̄, ȳ: means; sums over all stars

The sign of a latent is arbitrary (replacing L by −L and flipping its coefficients gives the same model), so only |r| matters.

### 8.4 Tools used to test the method
#### 8.4.1 Counting equations and the Jacobian: can the coefficients be uniquely determined?
- **Counting**: the data provide only the 15 correlations in S (pairs of the 6 quantities); the authors' graph has 15 edges, i.e. 15 unknowns. Having no more unknowns than equations is necessary for uniqueness, but not sufficient: the equations may be dependent.
- **Jacobian matrix**: J_mn = ∂ρ_m / ∂θ_n
  Meaning: how much the m-th model correlation changes when the n-th coefficient changes slightly.
  - ρ_m: the m-th model correlation (m = 1…15)
  - θ_n: the n-th edge coefficient (n = 1…15)
  - ∂: partial derivative, the rate of change when only θ_n changes

  J is a 15×15 matrix. Its **rank** is the number of independent directions in coefficient space that the data can distinguish. Rank = 15: every change leaves a trace in the data, the coefficients are (locally) unique; **rank = 14**: there is a direction v (with J·v = 0) along which the coefficients can change while all 15 correlations stay the same — invisible to the data. J is computed by finite differences, the rank by singular-value decomposition (singular values below 1e-8 times the largest count as 0), and v is reported.
- Since this depends only on the graph, any random coefficient set gives the same answer; all 3 random sets gave 14.

#### 8.4.2 Multistart fits, profiles, and "statistically indistinguishable"
- **Multistart**: maximum likelihood with L-BFGS-B (a method that really converges to a local optimum) from many random coefficient sets, recording each solution's NLL and |r(L1, Rbirth)|. If the model were identifiable with a single maximum, all starts would converge to the same solution and the same L1.
- **Profile**: fix the coefficient with the largest weight in the non-identifiable direction (L1→age) at a series of values, refit the other 14, and plot "best NLL vs L1→age" and "|r(L1, Rbirth)| vs L1→age".
- **What "equally good" means — the likelihood-ratio test**: whether the difference in fit between two coefficient sets "really matters" is judged by

  Λ = 2 · N · (NLL₂ − NLL₁)

  Meaning: twice the difference in log-likelihood of the two solutions (the NLL is per star, hence the factor N). If Λ < 3.84, solution 1 cannot be said to be better than solution 2 at 95% confidence; the two are "statistically indistinguishable".
  - N: number of stars; NLL₁: NLL of the best solution; NLL₂: NLL of the other solution
  - 3.84: the 95% quantile of the χ² distribution with 1 degree of freedom (the standard threshold for models differing by one degree of freedom)

  The round-1 criterion "ΔNLL < 0.01" corresponds to Λ < 0.02·N, i.e. Λ < 200 for 10,000 stars, far looser than 3.84; round 2 uses the strict criterion above.
- **Noise variance ≥ 0**: physical requirement (8.3.5). Round 2 imposes it as a hard constraint (via a large penalty). An optimum lying exactly on the boundary "some noise variance = 0" is the Heywood case of Sec. 4.3.
- **Strict convergence (polish)**: when an optimizer stops, its NLL can still be off by ~10⁻⁴, enough to change Λ for tens of thousands of stars. Round 2 therefore re-optimizes every solution with very strict stopping criteria to make this error negligible.
- The |r| values are used only to show ranges; **they are never used to select a solution** (selection never looks at the truth).
- **Limitation**: the number of starts is finite (230 per galaxy). If a better optimum exists that was not found, the "best" changes, and so does the set of solutions indistinguishable from it (g2.79e12 is an example), so the strict range is "what has been found so far", not an exact confidence interval.

#### 8.4.3 Linear ceiling and single-variable baselines
- **Linear ceiling**: the correlation between true Rbirth and its least-squares linear regression on the 6 standardized inputs (the multiple correlation coefficient). By the consequence in 8.3.6 it bounds any L1.
- **Single-variable baselines**: |r| of age alone or [Fe/H] alone with Rbirth (and J_φ alone for wide selections). A method that cannot beat these baselines cannot be said to "find the birth radius".

#### 8.4.4 Covariance matrices are enough (how the large data of Sec. 3.12 were handled)
- **Maximum likelihood uses only the covariance**: the NLL of the linear Gaussian model depends only on the sample correlation matrix S of the 6 inputs and the number of stars N (8.3.5), not on the individual stars. Compressing each galaxy into a covariance matrix therefore gives exactly the same fit as using all stars.
- **So does the correlation of L1 with the truth**: by 8.3.6, L1 is a linear combination of the standardized inputs, L1 = Σ_i w_i · z_i, so

  r(L1, R) = (w · c) / √[(w · C · w) · var(R)]

  Meaning: the covariance of L1 with the true birth radius divided by the product of their standard deviations (the definition of the Pearson correlation, written with matrices).
  - z_i: the i-th input minus its mean, divided by its standard deviation (i = 1…6)
  - w: the vector of the 6 weights w_i, set by the edge coefficients and the hierarchical least squares (the same for all stars)
  - C: the 6×6 correlation matrix of the standardized inputs; c: the vector of covariances of the 6 standardized inputs with the true birth radius R
  - var(R): variance of the true birth radius; "·": vector/matrix products

  So each galaxy needs only the 8×8 covariance matrix of 8 quantities (6 inputs + Rbirth + Rg) and N. Checked on the 7–10 kpc data: identical to the star-by-star computation to 6 decimal places.
