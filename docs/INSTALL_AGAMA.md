# Installing AGAMA (only needed to recompute orbits)

AGAMA is used only by `scripts/local/compute_orbits.py` (potential fitting, actions, guiding radius, eccentricity). The public data in `data/release/` already contains these quantities for the main selection, so most users never need AGAMA.

The PyPI package `agama` is only a 1 kB bootstrap and fails with "output.json not found"; build from source instead. Recipe that worked on macOS (Apple Silicon) inside a virtual environment:

```bash
xcode-select -p                         # make sure the compiler is installed
python -m pip install wheel             # needed by the build
brew install libomp                     # a recent OpenMP runtime
git clone https://github.com/GalacticDynamics-Oxford/Agama.git && cd Agama
rm -f Makefile.local && make clean      # always remove a stale Makefile.local before rebuilding
# If you use Anaconda: hide it during the build so the old OpenMP in Anaconda is not picked up
env -u CONDA_EXE -u CONDA_PREFIX -u CONDA_DEFAULT_ENV \
  python -m pip install . --no-build-isolation --config-settings --build-option=--yes
grep -n omp Makefile.local              # should point to /opt/homebrew/opt/libomp
cd ~ && python -c "import agama; print(agama.Potential(type='Plummer').potential([1,0,0]))"   # ≈ -0.707
```

Pitfalls we hit:
- Anaconda's `libiomp5.dylib` (llvm-openmp 14) lacks `__kmpc_dispatch_deinit`, which recent clang needs → link error. Hiding the conda variables makes the installer use Homebrew's libomp.
- `--yes` silently reuses an existing `Makefile.local`; delete it before rebuilding.

More details: `docs/REPRODUCTION.md`, Sec. 8.
