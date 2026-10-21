# Project status

- Last research update: **2026-09-27**. Packaged as a public repository: 2026-10-21.
- One-line summary: **RLCD, applied to TNG50 as in the paper (frozen authors' graph, re-estimated coefficients), does not reliably recover the stellar birth radius**; learning the graph from scratch (Approach B) does not produce a birth-radius latent either. All planned analyses are complete; the project is at a decision point on what to do next (Sec. 4).
- Details: [`REPORT.md`](REPORT.md) (results and method), [`docs/REPRODUCTION.md`](docs/REPRODUCTION.md) (commands, data formats, decisions, pitfalls).
- When the project advances, update this file and Sec. 2 of `docs/REPRODUCTION.md`.

## 1. Done

### Stage 1 — TNG50 data (complete, 2026-09-26)
- [x] server extraction of particles, ages, abundances and birth radii of 12 galaxies → `data/tng50_extract/data_obtain_all.hdf5` (not public)
- [x] star-by-star cross-check with earlier birth-radius results → `results/check_vs_001_orbit.txt`
- [x] potentials, actions, guiding radii, eccentricities with AGAMA → `data/tng50_derived/orbits_all.hdf5` (not public); all consistency checks pass → `results/orbit_check/`
- [x] main selection (7 < R < 10 kpc, |z| < 0.3 kpc, [Fe/H] > −0.5) and RLCD inputs → `data/rlcd_input/{nonoise,noise}/` (public, `data/release/`)

### Stage 2 — RLCD vs truth (all planned analyses complete, 2026-09-27)
- [x] RLCD installed; the authors' causal graph reproduced exactly with `PYTHONHASHSEED=1` → `results/rlcd_g279/`, `results/reproduce_Jin/`
- [x] **Approach A** (authors' graph frozen, coefficients re-estimated on 12 galaxies) + strict identifiability analysis → `results/rlcd_A/identifiability/`, `results/TNG/`
- [x] Approach A for three wider selections → `star_condition/R_0_30_Z_30/`, `R_0_30_Z_3/`, `R_0_30_Z_0.3/`
- [x] **Approach B** (structure learning from scratch, 30 seeds per galaxy) → `PlanB/`
- [x] report → `REPORT.md`

### Stage 3 — Milky Way data (on hold)
- [ ] Only meaningful if the method works on simulations, which it currently does not (Sec. 3).

## 2. Problems found (not reported in the paper or the code repository)
1. **Structure learning is random**: the output graph depends on `PYTHONHASHSEED`; the authors' graph appears in 7 of 30 seeds, another graph in 12. Causes: iteration order of Python sets of variable names and an unseeded `random.sample` in `LatentGroups.py`.
2. **Parameter selection is flawed**: the NLL is computed with the absolute value of a determinant and negative noise variances are allowed, so unphysical solutions can look better; the paper's coefficients are far from the maximum-likelihood solution.
3. **(Most fundamental) the edge coefficients of the frozen graph are not identifiable**: Jacobian rank 14 of 15; the optimum sits on a "noise variance = 0" boundary (Heywood case) with many nearly equivalent optima. Solutions statistically indistinguishable from the best give |r(L1, Rbirth)| from ~0 to 0.4–0.5 in most galaxies, and even the most favourable one never beats "age alone" or "[Fe/H] alone" (12/12 galaxies without errors).
   - "Statistically indistinguishable" means 2·N·ΔNLL < 3.84: the difference in fit quality is too small to prefer one solution at 95% confidence. N is the number of stars, ΔNLL the per-star negative log-likelihood above the best solution, and 3.84 the 95% quantile of the χ² distribution with one degree of freedom.
4. **Wider selections do not help**: over the whole disk |r(L1, Rbirth)| rises to ~0.8, but J_φ alone (≈ current radius) gives 0.8–0.9; the maximum-likelihood solution still jumps with the selection.
5. **Approach B**: the authors' graph is never learned (0/360); only 3 galaxies have an identifiable structure; the best-matching latent (chosen with the truth) has median |r| = 0.51 and rarely beats age alone.

## 3. Conclusion so far
The method, as published, is not validated by TNG50, and the non-identifiability (problem 3) is a property of the causal graph that would equally affect Milky Way data, where no ground truth is available to detect it.

## 4. Open questions and possible next steps
1. **Make the model identifiable**: add observables (J_R, more element abundances) so that there are more correlations than edges; or a physically motivated constraint that does not use the truth; or non-linear relations. Check full Jacobian rank first, then compare with baselines.
2. **Always report baselines**: age alone, [Fe/H] alone, J_φ alone, and the linear-regression ceiling of the inputs (`scripts/local/linear_ceiling.py`).
3. **Optional**: Approach B for the wide selections (supported by `PlanB/output_data/rlcd_B.py --cond ...`, not run).
4. **Small open points**: eccentricity from 3D r (current) or cylindrical R; independent checks such as energy conservation along orbits and J_z from orbit integration vs the Stäckel approximation.
5. **Contact the original authors** about the randomness, the coefficient selection and the identifiability issue.
6. Milky Way application (stage 3): only after the points above.
