# Selection R_0_30_Z_3 (2026-09-27)

0 < R < 30 kpc, |z| < 3 kpc; all other cuts as in the main selection (bound with all orbit quantities, birth radius < 20 kpc, [Fe/H] > −0.5, J_z > 0).
0.371–0.699 million stars per galaxy, 6.02 million in total (the stars of R_0_30_Z_30 with |z| < 3 kpc, see `input_data/make_input.py`).

| Folder | Content | Produced by |
|---|---|---|
| `input_data/` | `nonoise/`, `noise/<ID>_e.csv` (6-column RLCD input; not distributed, see [`DATA.md`](../../DATA.md)); `selection.txt` (selection statistics) | `input_data/make_input.py` (cut from R_0_30_Z_30) |
| `true_value/` | `nonoise/`, `noise/<ID>_truth.csv` (truth: birth radius R_bir, guiding radius Rg, ...; not distributed) | same |
| `output_data/` | `moments_<version>.json` (8×8 covariance matrix of 8 quantities per galaxy); `<version>_<ID>_ci.csv/_summary.txt` (maximum-likelihood and indistinguishable solutions); `summary_ci_<version>.csv/.txt/.png` (summary of the 12 galaxies) | `compute_moments.py` → `identifiability_moments.py` |
| `causal_graph/` | authors' causal graph + maximum-likelihood coefficients of each galaxy | `causal_graph_R_0_30_Z_3.py` |
| `truth_comparison/` | f(L1) vs birth radius, f(L2) vs guiding radius | `L1_Rbirth_L2_Rg_R_0_30_Z_3.py` |

Order (from the repo root, with the analysis environment `$PY`): `make_input.py` → `compute_moments.py` → `identifiability_moments.py --version nonoise` (and `noise`) → the two plotting scripts (add `--version noise` for `_noise.png`).
`identifiability_moments.py` needs only the committed `moments_<version>.json`, so it can be rerun without the large data.

⚠ This selection covers the whole disk, where the current angular momentum J_φ (≈ current radius) is already strongly correlated with the birth radius (0.78–0.90 in most galaxies), so L1 must be compared with "J_φ alone". The paper uses only stars at 7–10 kpc precisely to remove this information.
