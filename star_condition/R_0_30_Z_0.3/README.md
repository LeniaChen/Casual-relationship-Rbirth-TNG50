# Selection R_0_30_Z_0.3 (2026-09-27)

0 < R < 30 kpc, |z| < 0.3 kpc; all other cuts as in the main selection (bound with all orbit quantities, birth radius < 20 kpc, [Fe/H] > −0.5, J_z > 0).
0.188–0.313 million stars per galaxy, 2.79 million in total (the stars of R_0_30_Z_30 with |z| < 0.3 kpc, see `input_data/make_input.py`).
The |z| limit is the paper's (the paper uses 7 < R < 10 kpc, |z| < 0.3 kpc); only R is widened to the whole disk.

| Folder | Content | Produced by |
|---|---|---|
| `input_data/` | `nonoise/`, `noise/<ID>_e.csv` (6-column RLCD input; not distributed, see [`DATA.md`](../../DATA.md)); `selection.txt` (selection statistics) | `input_data/make_input.py` (cut from R_0_30_Z_30) |
| `true_value/` | `nonoise/`, `noise/<ID>_truth.csv` (truth: birth radius R_bir, guiding radius Rg, ...; not distributed) | same |
| `output_data/` | `moments_<version>.json` (8×8 covariance matrix of 8 quantities per galaxy); `<version>_<ID>_ci.csv/_summary.txt` (maximum-likelihood and indistinguishable solutions); `summary_ci_<version>.csv/.txt/.png` (summary of the 12 galaxies) | `compute_moments.py` → `identifiability_moments.py` |
| `causal_graph/` | authors' causal graph + maximum-likelihood coefficients of each galaxy | `causal_graph_R_0_30_Z_0.3.py` |
| `truth_comparison/` | f(L1) vs birth radius, f(L2) vs guiding radius | `L1_Rbirth_L2_Rg_R_0_30_Z_0.3.py` |

Order (from the repo root, with the analysis environment `$PY`): `make_input.py` → `compute_moments.py` → `identifiability_moments.py --version nonoise` (and `noise`; set `OMP_NUM_THREADS=1` when running both at once, it is more than twice as fast) → the two plotting scripts (add `--version noise` for `_noise.png`).
`identifiability_moments.py` needs only the committed `moments_<version>.json`, so it can be rerun without the large data.

## Key results (nonoise)
- |r(L1, birth radius)| of the maximum-likelihood solution: 0.79–0.89 in 9 galaxies, 0.45 for 552581, 0.28 for 569251, 0.17 for 546870.
- In no galaxy does the maximum-likelihood solution beat "J_φ alone" (J_φ alone vs birth radius: 0.71–0.92); only 538905 and 546870 exceed it at the upper end of the indistinguishable solutions.
- There are still many statistically indistinguishable solutions (2–79), and the |r| range of L1 is wide in 487742 (0.22–0.79), 560421 (0.19–0.80), 546870 (0.06–0.63), ... → the coefficients are still not unique (Jacobian rank 14/15).
- Compared with |z|<3, the optimum changed for 557721 (0.58→0.89), 569251 (0.05→0.28), 560421 (0.69→0.79): the maximum-likelihood solution is very sensitive to the selection.
- The noise version is even less stable: the best |r| drops to 0.10–0.41 for 545437, 547293, 557721, 560421, while indistinguishable solutions reach 0.77–0.90.

⚠ This selection covers the whole disk, where the current angular momentum J_φ (≈ current radius) is already strongly correlated with the birth radius, so L1 must be compared with "J_φ alone".
