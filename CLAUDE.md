# Notes for AI coding assistants (and new contributors)

This repository tests whether the RLCD causal-discovery method of Jin et al. (arXiv:2507.00134) recovers stellar birth radii, using TNG50 galaxies where the truth is known. `AGENTS.md` points here.

## Orient yourself first
1. Read `STATUS.md` (one page: status, problems found, open questions).
2. For details read `docs/REPRODUCTION.md`: Sec. 0 (start-up checklist), Sec. 2 (progress), Sec. 6 (every pipeline step and data format), Sec. 9 (pitfalls). For the scientific argument read `REPORT.md`.
3. Check the setup: `bash tools/check_setup.sh` (public data: `bash tools/unpack_data.sh`; full data, if available: `bash tools/link_data.sh <folder>`).
4. When you finish a piece of work, tell the user where the project stands and what comes next.

## Conventions
- Reply in the language the user writes in. Repository documentation is in English.
- Documentation style: every formula is followed by a plain-language explanation of what it means and of each symbol.
- Report progress during long tasks instead of staying silent until the end.
- After any research progress, update `STATUS.md` and Sec. 2 of `docs/REPRODUCTION.md`, detailed enough that someone without context can reproduce and continue.
- Approach B results live in `PlanB/README.md` and `docs/REPRODUCTION.md`; `REPORT.md` covers Approach A.
- Two Python environments: the analysis environment (`$PY`, `setup/requirements.txt`, AGAMA optional) and the RLCD environment (`$PY_RLCD`, Python 3.12, `setup/requirements-rlcd.txt`). They are defined in `.env` (copy of `env.example`); load with `set -a; . ./.env; set +a`. Install new packages into the matching environment, never into a base/system Python. Prefer giving the user complete, copyable commands for anything that must run on their machine, with output teed to a log in `results/`.
- The TNG JupyterLab server has an **old** AnastrisTNG; server scripts must not assume the 2.0.0 API. Particle data are downloaded and processed locally.

## Where things are
- All scripts locate files relative to the repository root (`ROOT = Path(__file__).parents[n]`); run them from the root.
- Large data (~17 GB) are not in Git (`.gitignore`, `DATA.md`, `data_manifest.tsv`). Never commit files > 50 MB; if a new large output appears, add it to `.gitignore`, `tools/link_data.sh` and `data_manifest.tsv`.
- The authors' code `Causal-Galactic-Archaeology/` is not in Git (`tools/setup_rlcd_repo.sh`, commit 8c0a5a4). **Running RLCD inside its `rlcd/` folder overwrites the authors' committed results**: move new outputs to `results/` and restore originals with `git checkout`.
- Runnable without the large data: `star_condition/*/output_data/identifiability_moments.py` (reads the committed `moments_*.json`), everything using `data/rlcd_input/` (public), `results/reproduce_Jin/` (needs the authors' repository).

## Technical rules that matter
- RLCD structure learning is random: always set `PYTHONHASHSEED=<n>`, repeat with several seeds and report structure frequencies.
- Edge coefficients: require noise variances ≥ 0, use a converging optimizer (L-BFGS-B) with hundreds of starts, report every solution statistically indistinguishable from the best (2·N·ΔNLL < 3.84), and never use the truth to select a solution.
- Always compare |r| with baselines: age alone, [Fe/H] alone, J_φ alone, and the linear-regression ceiling of the six inputs.
- Set `OMP_NUM_THREADS=1` when running two analysis processes at once.
- Without h5py the scripts fall back to the pure-Python reader `pyfive`.

## Document index
| File | Purpose |
|---|---|
| `README.md` | overview, key results, quick start, data, citation |
| `STATUS.md` | status, problems, next steps (read first) |
| `REPORT.md` | full report: steps, problems, results, conclusions, method |
| `docs/REPRODUCTION.md` | reproduction log: commands, data structures, decisions, installation, pitfalls |
| `DATA.md`, `data_manifest.tsv` | public vs non-public data, sizes, checksums, regeneration |
| `PlanB/README.md`, `star_condition/*/README.md`, `results/rlcd_A/README.md` | sub-project notes and key results |
