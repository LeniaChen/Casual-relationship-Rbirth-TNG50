# Approach B: learn the causal graph from scratch on each TNG50 galaxy (2026-09-27)

**Question**: Approach A takes the graph the authors learned on g2.79e12, freezes it and transfers it to TNG50. Approach B asks the opposite: if RLCD learns from scratch on each TNG50 galaxy, which graph does it find, and are its latents the birth radius and the guiding radius?

**Selection**: same as the paper and the main result of Approach A: 7 < R < 10 kpc, |z| < 0.3 kpc (7.3–35 thousand stars per galaxy, 240 thousand in total). The scripts also support `--cond R_0_30_Z_30 / R_0_30_Z_3 / R_0_30_Z_0.3` (output to `output_data/<selection>/`, not run yet).

| Folder | Content | Produced by |
|---|---|---|
| `input_data/` | `nonoise/`, `noise/<ID>_e.csv` (6-column RLCD input); `selection.txt` | `input_data/make_input.py` (copied unchanged from `data/rlcd_input`) |
| `true_value/` | `nonoise/`, `noise/<ID>_truth.csv` (true birth radius R_bir, true guiding radius Rg, ...) | same |
| `output_data/` | `structures_<version>.csv/.txt`: graphs learned with 30 seeds per galaxy; `eval_<version>.csv/.txt/_best.json`: Jacobian rank, maximum-likelihood coefficients, |r| of each latent with the truth for every main structure; `summary_planB_<version>.csv/.txt`: one-row-per-galaxy summary compared with Approach A | `rlcd_B.py` → `rlcd_B_eval.py` → `summary_planB.py` |
| `causal_graph/` | most frequent graph of each galaxy + maximum-likelihood coefficients | `planB_causal_graph.py` (`--version noise`; `--rank 2` for the second most frequent graph) |
| `truth_comparison/` | best-matching latent f(L) vs true birth radius / true guiding radius per galaxy | `planB_latent_vs_truth.py` |

Order (from the repo root; `$PY` = analysis environment, `$PY_RLCD` = RLCD environment):
1. `$PY PlanB/input_data/make_input.py` (after `bash tools/unpack_data.sh`)
2. `$PY_RLCD PlanB/output_data/rlcd_B.py --version nonoise --seeds 30` (and `noise`; about 2 s per RLCD run, about 12 minutes for 12 galaxies × 30 seeds; finished runs are skipped)
3. `$PY PlanB/output_data/rlcd_B_eval.py --version nonoise` (and `noise`; about 10 minutes; set `OMP_NUM_THREADS=1` when running both at once)
4. `summary_planB.py` and the two plotting scripts (add `--version noise` for the noisy version)

Details:
- Each seed runs in its own Python process with `PYTHONHASHSEED` set (the RLCD result depends on set iteration order); RLCD settings as the authors' (GES, sparsity 2, α 1e-5/1e-15/1e-5/1e-15). Latent labels are unified by the lexicographically smallest relabeling, so the names L1, L2 are arbitrary.
- Evaluation: for each galaxy, every structure with frequency ≥ 20%: 50 random starts for maximum likelihood (noise variances ≥ 0), keeping all solutions statistically indistinguishable from the best (2·N·ΔNLL < 3.84; N = number of stars, ΔNLL = per-star NLL above the best). Latents are computed from their children by the hierarchical least squares of the authors' notebook.
- Check: with the authors' graph, `rlcd_B_eval.py` on 538905 gives NLL 0.4918, Jacobian rank 14/15, and L1/L2 |r| ranges within those of Approach A (consistent with Approach A).
- Rerunning `summary_planB.py` on two different machines (macOS, Linux) gives `summary_planB_nonoise.csv` files differing only in the last floating-point digit (the .txt files are identical); all other files have identical checksums.

## Results (nonoise)
- **Identical to the authors' graph: 0 of 360 runs** (also 0/360 for noise). In 8 galaxies all 30 seeds give the same graph, in 4 there are 2 kinds.
- The learned graphs vary widely: 0–4 latents, 8–17 directed edges. 545437 and 569251 have no latent and undirected edges, so no coefficients can be fitted. 4 galaxies violate the paper's prior (age has an observed parent).
- Of the 10 galaxies where coefficients can be fitted, only 3 have an identifiable main structure (Jacobian rank = number of edges): 546870, 562029, 567607, each with a single latent; the other 7 are still not identifiable.
- Picking, with the truth, the latent most similar to the birth radius (best case): median |r| 0.51 (0.24–0.77); only 530330 (0.71 vs 0.68), 552581 (0.53 vs 0.23) and 562029 (0.77 vs 0.69) beat "age alone". The guiding radius is usually matched better (mostly 0.6–0.9).
- After a BIC penalty for the different number of edges (dBIC = (NLL_B − NLL_A) + (edges_B − 15)·ln N/(2N); positive means the authors' graph is better), **the authors' graph fits better than the graph RLCD learns on that galaxy in every comparable galaxy** (dBIC per star +0.003 to +0.42) — RLCD's rank-test-based structure learning does not find the best graph in the maximum-likelihood sense.
  - Meaning of the formula: the difference in per-star negative log-likelihood between the learned graph (B) and the authors' graph (A), plus a penalty for every edge beyond the authors' 15.
  - NLL_B, NLL_A: per-star negative log-likelihoods (smaller = better fit); edges_B: number of edges of the learned graph; N: number of stars; ln: natural logarithm.
- noise: simpler structures (mostly 1–2 latents), median |r(birth radius)| of the best-matching latent only 0.27; dBIC positive everywhere except 487742 (−0.007).

**Conclusion**: on TNG50, RLCD does not learn a graph with "L1 = birth radius, L2 = guiding radius" by itself; the truth is needed to tell which latent looks like the birth radius, and even then it is mostly no better than age alone.

(Approach B is a supplementary analysis and is not discussed in `REPORT.md`.)
