# -*- coding: utf-8 -*-
"""
compute_orbits.py — fit the gravitational potential of each galaxy with AGAMA and compute for all stars:
    actions (J_R, J_z, J_phi), energy E, guiding radius Rg, peri-/apocentre and eccentricity (orbit integration), z_max

Run (from the repo root, with the analysis environment):
    cd <repo root>
    $PY scripts/local/compute_orbits.py --test     # test: 1 galaxy, 5000 random stars
    $PY scripts/local/compute_orbits.py            # all 12 galaxies, all stars
    $PY scripts/local/compute_orbits.py --ids 487742 530330   # selected galaxies only

Input:   data/tng50_extract/data_obtain_all.hdf5
Outputs: data/tng50_derived/orbits_all.hdf5          (test mode: orbits_test.hdf5)
           /<ID>/Jr, Jz, Jphi, E, Rg, rperi, rapo, ecc, zmax, n_period  — each array has one entry per star,
                                                                          in exactly the order of /<ID>/stars in data_obtain_all.hdf5
           /<ID>/rotation_curve/R, vc                  — circular-velocity curve of the fitted potential (for checks)
         data/tng50_derived/potentials/<ID>_halo.coef, <ID>_disk.coef  — potential expansion coefficients, reload with agama.Potential(file=...)
Galaxies already done are skipped; after an interruption simply rerun.

Method (as in Jin et al.: AGAMA + axisymmetric multipole expansion including dark matter, gas and stars):
  - frame: the z=0 disk frame saved by data_obtain (centre at the origin, stellar angular momentum along +z)
  - potential = sum of two parts, both forced axisymmetric; inner grid boundary rmin = 0.2 kpc (about the TNG50
    gravitational softening), smoothly extrapolated by AGAMA inside (if ActionFinder still fails, 0.5 and 1.0 kpc are tried):
      "halo" part: Multipole expansion (lmax=8): dark matter, black holes, wind particles, gas and stars outside the disk
      "disk" part: CylSpline expansion (cylindrical grid): stars with r<50 kpc and gas with R<30 kpc and |z|<5 kpc
    (the standard AGAMA recipe for simulated galaxies: a flat disk is described more accurately on a cylindrical grid)
  - actions: AGAMA ActionFinder (Staeckel fudge); J_phi is the z-component of the angular momentum, L_z
  - Rg: the radius where L_z = Rg * v_c(Rg), AGAMA pot.Rcirc(L=|L_z|)
  - eccentricity: each orbit is integrated for N_PERIOD circular-orbit periods in this potential; min/max of the 3D distance
        e = (r_apo - r_peri) / (r_apo + r_peri)
  - unbound stars (E >= 0): all orbit quantities are NaN
Units: length kpc, velocity km/s, mass Msun, actions kpc*km/s, energy (km/s)^2
"""
import os
import sys
import time
import argparse
import numpy as np
import h5py
import agama

# ==================== parameters ====================
IN_FILE = 'data/tng50_extract/data_obtain_all.hdf5'
OUT_DIR = 'data/tng50_derived'
POT_DIR = os.path.join(OUT_DIR, 'potentials')

N_PERIOD = 15        # integration time = N_PERIOD * circular-orbit period Tcirc(E)
TRAJ_SIZE = 1501     # points stored per orbit (about 100 per circular period)
CHUNK = 5000         # orbits per integration batch (memory: 5000*1501*6*4 bytes ≈ 180 MB)
ACCURACY = 1e-7      # relative accuracy of the orbit integration

# particles assigned to the "disk" (CylSpline) part
DISK_STAR_RMAX = 50.0    # stars: 3D distance r < 50 kpc
DISK_GAS_RMAX = 30.0     # gas: cylindrical radius R < 30 kpc
DISK_GAS_ZMAX = 5.0      #      and |z| < 5 kpc

N_TEST = 5000

# Inner boundary of the potential grids (kpc). The TNG50 softening at z=0 is ~0.3 kpc; the mass distribution is unreliable below that.
# The automatic grid extends to ~0.01 kpc, where particle noise makes the circular angular momentum non-monotonic and ActionFinder fails.
# These inner boundaries are tried in turn until ActionFinder initializes.
RMIN_TRY = [0.2, 0.5, 1.0]
RMAX_HALO = 500.0        # outer boundary of the Multipole (kpc)
RMAX_DISK = 60.0         # outer boundary in R and z of the CylSpline (kpc)
ZMIN_DISK = 0.05         # first CylSpline grid node in z (kpc)
POT_VERSION = 'v2-rmin-grid'   # bump when the potential settings change; galaxies stored with another version are recomputed

agama.setUnits(mass=1, length=1, velocity=1)   # Msun, kpc, km/s (time unit = 0.978 Gyr)
TIME_UNIT_GYR = 0.9778


def build_potential(g, gid, rmin):
    """Build the axisymmetric potential from dm / gas / stars / wind / bh particles (inner grid boundary rmin kpc)"""
    rd = lambda grp: (g[grp]['pos'][:].astype(np.float64), g[grp]['mass'][:].astype(np.float64))
    s_pos, s_m = rd('stars')
    g_pos, g_m = rd('gas')
    d_pos, d_m = rd('dm')
    w_pos, w_m = rd('wind')
    b_pos, b_m = rd('bh')

    s_disk = np.linalg.norm(s_pos, axis=1) < DISK_STAR_RMAX
    g_disk = (np.hypot(g_pos[:, 0], g_pos[:, 1]) < DISK_GAS_RMAX) & (np.abs(g_pos[:, 2]) < DISK_GAS_ZMAX)

    disk_pos = np.vstack([s_pos[s_disk], g_pos[g_disk]])
    disk_m = np.concatenate([s_m[s_disk], g_m[g_disk]])
    halo_pos = np.vstack([d_pos, b_pos, w_pos, s_pos[~s_disk], g_pos[~g_disk]])
    halo_m = np.concatenate([d_m, b_m, w_m, s_m[~s_disk], g_m[~g_disk]])

    print(f"    potential (rmin = {rmin} kpc): disk part {len(disk_m):,} particles ({disk_m.sum():.3e} Msun), "
          f"halo part {len(halo_m):,} particles ({halo_m.sum():.3e} Msun)")

    pot_halo = agama.Potential(type='Multipole', particles=(halo_pos, halo_m),
                               symmetry='axisymmetric', lmax=8, gridSizeR=40,
                               rmin=rmin, rmax=RMAX_HALO)
    pot_disk = agama.Potential(type='CylSpline', particles=(disk_pos, disk_m),
                               symmetry='axisymmetric', mmax=0, gridSizeR=30, gridSizez=30,
                               rmin=rmin, rmax=RMAX_DISK, zmin=ZMIN_DISK, zmax=RMAX_DISK)
    os.makedirs(POT_DIR, exist_ok=True)
    pot_halo.export(os.path.join(POT_DIR, f'{gid}_halo.coef'))
    pot_disk.export(os.path.join(POT_DIR, f'{gid}_disk.coef'))
    return agama.Potential(pot_halo, pot_disk)


def peri_apo(pot, xv, T):
    """Integrate orbits in batches; return r_peri, r_apo, z_max"""
    n = len(xv)
    rperi = np.full(n, np.nan)
    rapo = np.full(n, np.nan)
    zmax = np.full(n, np.nan)
    t0 = time.time()
    for i0 in range(0, n, CHUNK):
        i1 = min(i0 + CHUNK, n)
        _, traj = agama.orbit(potential=pot, ic=xv[i0:i1], time=T[i0:i1],
                              trajsize=TRAJ_SIZE, accuracy=ACCURACY,
                              separateTime=True, verbose=False)
        traj = np.asarray(traj, dtype=np.float32)          # (m, TRAJ_SIZE, 6)
        r = np.sqrt(np.sum(traj[:, :, :3].astype(np.float64)**2, axis=2))
        rperi[i0:i1] = r.min(axis=1)
        rapo[i0:i1] = r.max(axis=1)
        zmax[i0:i1] = np.abs(traj[:, :, 2]).max(axis=1)
        done = i1
        el = time.time() - t0
        print(f"\r    orbit integration {done:,}/{n:,}  elapsed {el/60:.1f} min, "
              f"remaining ~{el/done*(n-done)/60:.1f} min   ", end='', flush=True)
    print()
    return rperi, rapo, zmax


def process_galaxy(fin, fout, gid, test=False, rng=None):
    g = fin[gid]
    t_start = time.time()
    print(f"\n==== galaxy {gid} ====")

    pot = af = None
    for rmin in RMIN_TRY:
        pot = build_potential(g, gid, rmin)
        try:
            af = agama.ActionFinder(pot, interp=False)
            break
        except RuntimeError as err:
            print(f"    ActionFinder initialization failed with rmin = {rmin} kpc: {err}")
            print("    retrying with a larger rmin")
    if af is None:
        raise RuntimeError(f"galaxy {gid}: all rmin values failed {RMIN_TRY}")

    pos = g['stars']['pos'][:].astype(np.float64)
    vel = g['stars']['vel'][:].astype(np.float64)
    n_all = len(pos)
    if test:
        idx = np.sort(rng.choice(n_all, size=min(N_TEST, n_all), replace=False))
    else:
        idx = np.arange(n_all)
    xv = np.hstack([pos[idx], vel[idx]])
    n = len(idx)
    print(f"    stars: {n:,}" + (f" (test: random subset of {n_all:,})" if test else ""))

    # ---- energy, bound or not ----
    E = pot.potential(xv[:, :3]) + 0.5 * np.sum(xv[:, 3:]**2, axis=1)
    bound = np.isfinite(E) & (E < 0)
    print(f"    bound stars: {bound.sum():,} / {n:,}")

    out = {k: np.full(n, np.nan) for k in ['Jr', 'Jz', 'Jphi', 'Rg', 'rperi', 'rapo', 'ecc', 'zmax']}
    out['E'] = E
    out['Jphi'] = xv[:, 0] * xv[:, 4] - xv[:, 1] * xv[:, 3]     # L_z, defined for all stars

    # ---- actions ----
    t0 = time.time()
    J = af(xv[bound])
    out['Jr'][bound] = J[:, 0]
    out['Jz'][bound] = J[:, 1]
    print(f"    actions done ({time.time()-t0:.1f} s)")

    # ---- guiding radius ----
    out['Rg'][bound] = pot.Rcirc(L=np.abs(out['Jphi'][bound]))

    # ---- orbit integration: peri-, apocentre, eccentricity ----
    T = N_PERIOD * pot.Tcirc(xv[bound])
    ok = np.isfinite(T) & (T > 0)
    b_idx = np.where(bound)[0][ok]
    rp, ra, zm = peri_apo(pot, xv[b_idx], T[ok])
    out['rperi'][b_idx] = rp
    out['rapo'][b_idx] = ra
    out['zmax'][b_idx] = zm
    out['ecc'][b_idx] = (ra - rp) / (ra + rp)

    # ---- circular-velocity curve ----
    Rgrid = np.linspace(0.5, 30, 60)
    pts = np.column_stack([Rgrid, 0 * Rgrid, 0 * Rgrid])
    vc = np.sqrt(-Rgrid * pot.force(pts)[:, 0])

    # ---- write ----
    grp = fout.create_group(gid)
    grp.attrs['n_stars'] = n
    grp.attrs['test_mode'] = bool(test)
    grp.attrs['pot_version'] = POT_VERSION
    grp.attrs['pot_rmin_kpc'] = rmin
    grp.attrs['N_PERIOD'] = N_PERIOD
    grp.attrs['TRAJ_SIZE'] = TRAJ_SIZE
    grp.attrs['units'] = 'kpc, km/s, Msun; actions kpc*km/s; E (km/s)^2'
    if test:
        grp.create_dataset('star_index', data=idx)
    for k, v in out.items():
        grp.create_dataset(k, data=v.astype(np.float32), compression='gzip',
                           compression_opts=4, shuffle=True)
    rc = grp.create_group('rotation_curve')
    rc.create_dataset('R', data=Rgrid)
    rc.create_dataset('vc', data=vc)

    # ---- quick diagnostics ----
    R = np.hypot(xv[:, 0], xv[:, 1])
    disk = bound & (R > 4) & (R < 12) & (np.abs(xv[:, 2]) < 1)
    def med(a): return np.nanmedian(a[disk]) if disk.sum() else np.nan
    print(f"    v_c(8 kpc) = {np.interp(8, Rgrid, vc):.1f} km/s,  max v_c = {vc.max():.1f} km/s")
    print(f"    disk stars (4<R<12, |z|<1): {disk.sum():,}; medians: "
          f"J_z = {med(out['Jz']):.1f}, J_phi = {med(out['Jphi']):.0f}, e = {med(out['ecc']):.3f}, "
          f"Rg - R = {med(out['Rg'] - R):+.2f} kpc")
    nan_e = np.sum(bound & ~np.isfinite(out['ecc']))
    if nan_e:
        print(f"    note: no eccentricity for {nan_e} bound stars")
    el = time.time() - t_start
    print(f"    galaxy {gid} done in {el/60:.1f} min")
    return el, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--test', action='store_true', help='only 5000 random stars of the first (or --ids) galaxy')
    ap.add_argument('--ids', nargs='*', help='only these galaxies')
    args = ap.parse_args()

    if not os.path.exists(IN_FILE):
        sys.exit(f"{IN_FILE} not found; run from the repo root")
    os.makedirs(OUT_DIR, exist_ok=True)
    out_file = os.path.join(OUT_DIR, 'orbits_test.hdf5' if args.test else 'orbits_all.hdf5')
    rng = np.random.default_rng(42)

    with h5py.File(IN_FILE, 'r') as fin:
        ids = sorted(fin.keys())
        if args.ids:
            ids = [i for i in ids if i in args.ids]
        if args.test:
            ids = ids[:1]
            if os.path.exists(out_file):
                os.remove(out_file)
        print(f"output: {out_file}\ngalaxies: {ids}")

        for gid in ids:
            with h5py.File(out_file, 'a') as fout:
                if gid in fout:
                    if fout[gid].attrs.get('pot_version', '') == POT_VERSION:
                        print(f"\ngalaxy {gid} already done, skipped")
                        continue
                    print(f"\ngalaxy {gid} was computed with an older potential setting; deleting and recomputing")
                    del fout[gid]
            # the galaxy group is written only after everything is computed; on error it is deleted, leaving no partial results
            with h5py.File(out_file, 'a') as fout:
                try:
                    el, n = process_galaxy(fin, fout, gid, test=args.test, rng=rng)
                except BaseException:
                    if gid in fout:
                        del fout[gid]
                    raise
            if args.test:
                print(f"\n[test] {n} stars took {el:.0f} s. At this speed all stars (~0.6–1 million per galaxy) "
                      f"take ~{el/n*8e5/60:.0f} min per galaxy (rough estimate; the potential is built once).")
    print("\nall done.")


if __name__ == '__main__':
    main()
