# results/rlcd_A — Approach A (authors' causal graph frozen)

Full explanation in [`REPORT.md`](../../REPORT.md) (Secs. 4 and 5).

## ⚠ Read this before using any of these results

**Problem 3 (the most fundamental): with the authors' causal graph frozen, the edge coefficients cannot be determined reliably from the data.**
1. Mathematically non-identifiable: the authors' graph has 15 edges and the 6 observables have 15 pairwise correlations, but the Jacobian of "coefficients → correlations" has rank 14. Along one direction (mainly L1→age trading off against L1→MgH) the coefficients can change while the observed correlation matrix stays exactly the same. This is a property of the graph, independent of the data.
2. A proper maximum-likelihood fit almost always pushes the noise variance of one observable to 0 (22 of the 24 TNG50 data sets, mostly J_φ; [Fe/H] or e for g2.79e12), and the likelihood surface has many local optima whose NLL differs by only ~1e-4.
3. Consequence: the meaning of L1 is unstable. In most TNG50 galaxies, solutions statistically indistinguishable from the best (2·N·ΔNLL < 3.84) give |r(L1, Rbirth)| from about 0 to 0.4–0.5; on g2.79e12, going from 30 to 230 starts changes the best solution's |r(L1, Rb)| from 0.25 to 0.62. **Even the most favourable solution essentially never correlates with Rbirth better than "age alone" or "[Fe/H] alone"** (no exception in the 12 nonoise galaxies; 2 noise galaxies are marginally higher at the favourable end: 0.58 vs 0.57 and 0.71 vs 0.67).

→ The r(L1, Rbirth) values of the single coefficient sets selected by `rlcd_A.py` in `g279/` and `nonoise/` are examples only, not evidence that the method works or fails.
→ The main result is `identifiability/summary_ci_<version>.txt/.png` (strict criterion).

**Correction (2026-09-26)**: the first summary (`identifiability/summary_<version>.*`) counted all solutions within 0.01 in NLL of the best as "equally good". The NLL is averaged per star, and for tens of thousands of stars this criterion is far too loose (it corresponds to a likelihood-ratio statistic < 0.02·N). Those files are kept for comparison; the conclusions are based on `summary_ci_*`.

**Problem 2: selection of the learned parameters.** The repository's `parameter_by_trek` computes the NLL with the absolute value of the determinant and does not check that the implied noise variances are ≥ 0, so physically impossible solutions can get a "better" NLL (−7.03 occurred on g2.79e12, below the theoretical bound −1.275). The coefficients used in the paper are the repository's set with the **second smallest** NLL (1.586), not the smallest (1.498), without explanation; recomputed with the correct formula its NLL is 1.593 versus about −1.272 at the optimum, a likelihood-ratio statistic of about 680,000 for 119,000 stars — far from the maximum-likelihood solution.

**Problem 1: structure learning is random.** With the same α, different PYTHONHASHSEED values give different causal graphs; the authors' graph appeared once in 5 runs (seed = 1) and 7 times in 30. See `results/rlcd_g279/`.

None of these three problems is reported in the paper or the repository.

## Contents
- `g279/`: `rlcd_A.py` v2 on the authors' data (g2.79e12) (the solution selected by RMSE gives r(L1,Rb) = 0.63, but the top 5 range from 0.01 to 0.63)
- `nonoise/538905/`: trial run of `rlcd_A.py` v2 on TNG50 galaxy 538905 (r(L1,Rb) = 0.06)
- `identifiability/`: **the main result**
  - round 1 (loose): `scripts/local/rlcd_identifiability.py` → `<tag>_summary.txt` (Jacobian rank, non-identifiable direction), `<tag>_multistart.csv` (solutions from 30 starts), `<tag>_profile.csv/.png` (profile with L1→age fixed); `rlcd_identifiability_summary.py` → `summary_<version>.csv/.txt/.png`
  - round 2 (strict): `scripts/local/rlcd_identifiability_ci.py` → `<tag>_ci.csv/.txt` (all solutions from 230 starts after strict convergence), `summary_ci_<version>.csv/.txt/.png`
  - reference: `scripts/local/linear_ceiling.py` → `linear_ceiling.csv/.txt` (linear-regression ceiling of the 6 quantities; correlations of age / [Fe/H] alone)
- The `|r(...)|` values are used only for diagnostics and ranges; **they were never used to select a solution**.
