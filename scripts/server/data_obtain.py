# -*- coding: utf-8 -*-
"""
data_obtain.py — TNG50 data extraction for the causal-archaeology (RLCD) project
Based on 001_orbit_new.py (final version); runs on the TNG JupyterLab server.

Differences from 001_orbit_new.py
------------------------------------------------------------
  kept: birth radius (galpos + TNG50_disk_trajectory.hdf5), abundance formulas [I][J],
        wind-particle removal, clamping of tform at the upper boundary [B], alignment assertion [A], young-star self-check [L],
        BirthPos must be extracted before face_on [F].
  new:
    [N] save stellar pos + vel in the z=0 disk frame (for actions Jz/Jphi, eccentricity e, Rg)
    [O] save pos + mass of dark matter / gas / black holes / wind particles (for the local AGAMA potential fit)
    [P] keep stars born "before the start of the trajectory" (birth_R = NaN, has_birth = False):
        they have no birth radius, but their mass is part of the potential and must not be dropped
    [Q] also store the raw mass fractions X_H / X_Mg / X_Fe, so a different solar scale needs no rerun
    [R] output one HDF5 per galaxy (float32 + gzip), smaller and faster than txt;
        written to a temporary file and renamed, so an interrupted run skips finished galaxies when restarted

Output: ./data_obtain/<ID>.hdf5
  /header (attrs)  ID, snap, a, h, time_Gyr, face_on parameters, units
  /frame           center_pos, center_vel, rotation   (return values of the z=0 face_on)
  /stars           iord, pos(N,3), vel(N,3), mass, init_mass, age, aform,
                   FeH, MgFe, X_H, X_Mg, X_Fe,
                   birth_pos(N,3), birth_R_2D, birth_R_3D, has_birth
  /wind            pos, mass                 (wind particles, only for the potential)
  /dm              pos, mass
  /gas             pos, mass
  /bh              pos, mass
Units: pos kpc (physical), vel km/s, mass Msol, age Gyr.
Frame: at z=0, origin at the galaxy centre and z axis along the stellar angular momentum (r<8 kpc) — identical to
Current_R_2D / Current_z of 001_orbit_new.py.
birth_pos is in the "disk frame at birth" given by galpos, which differs from the z=0 frame.
"""

from AnastrisTNG import TNGsimulation
import numpy as np
import os
import h5py
from galpos import GalaxyPoseTrajectory
from galpos import pynbody_decorate

# ==================== settings ====================
run = 'TNG50'
BasePath = '/home/tnguser/sims.TNG/' + run + '-1/output'
snap99 = TNGsimulation.Snapshot(BasePath, 99)

# Coordinates/Velocities/Masses/ParticleIDs of dm, star, gas, bh are loaded by default;
# extra fields are added for stars only
snap99.load_particle_para['star_fields'] = [
    'GFM_StellarFormationTime', 'BirthPos', 'BirthVel', 'GFM_InitialMass',
    'GFM_Metals'
]

orbit_file_path = 'TNG50_disk_trajectory.hdf5'   # relative to the working directory; may be an absolute path
if not os.path.exists(orbit_file_path):
    raise FileNotFoundError(
        f"trajectory file {os.path.abspath(orbit_file_path)} not found. Symlink TNG50_disk_trajectory.hdf5 "
        f"into the working directory (ln -s <original path> .) or edit orbit_file_path.")

# 607654 is not in the 423-galaxy sample (trajectory file) and was removed
halo_ids = np.array([487742, 530330, 538905, 545437, 546870,
                     547293, 552581, 557721, 560421, 562029,
                     567607, 569251])

FACE_ON_RMAX = 8          # kpc, as in 001_orbit_new.py
FACE_ON_ALIGN = 'star'

output_dir = './data_obtain/'
os.makedirs(output_dir, exist_ok=True)

print(f"processing {len(halo_ids)} galaxies: {halo_ids}")

# ==================== [I] abundance constants ====================
A_H, A_Mg, A_Fe = 1.008, 24.305, 55.845
LOGEPS_MG_SUN = 7.60      # Asplund et al. 2009
LOGEPS_FE_SUN = 7.50
SUN_FEH_TERM  = LOGEPS_FE_SUN - 12.0            # = -4.50
SUN_MGFE_TERM = LOGEPS_MG_SUN - LOGEPS_FE_SUN   # = +0.10

T_EDGE_TOL = 1e-3         # [B] Gyr

processing_errors = []
ledger_all = {}


def load_galaxy_trajectory(orbit_file, galaxy_id):
    with h5py.File(orbit_file, "r") as file:
        if str(galaxy_id) not in file:
            return None
        g = file[str(galaxy_id)]
        snap_all = np.sort(np.int16(list(g.keys())))
        return {
            'snap_all': snap_all,
            'ang_mom': np.array([g[f"{i}"]['ang_mom'][...] for i in snap_all]),
            'cpos':    np.array([g[f"{i}"]['cpos'][...]    for i in snap_all]),
            'cvel':    np.array([g[f"{i}"]['cvel'][...]    for i in snap_all]),
            't':       np.array([g[f"{i}"]['t'][...]       for i in snap_all]),
        }


def arr(x, unit=None, dtype=np.float32):
    """SimArray -> plain numpy array (optional unit conversion)"""
    if unit is not None:
        x = x.in_units(unit)
    return np.asarray(x.view(np.ndarray), dtype=dtype).copy()


def write_group(h5, name, data):
    grp = h5.create_group(name)
    for k, v in data.items():
        v = np.asarray(v)
        if v.ndim >= 1 and v.size > 1000:
            grp.create_dataset(k, data=v, compression='gzip',
                               compression_opts=4, shuffle=True)
        else:
            grp.create_dataset(k, data=v)


# ==================== main loop ====================
for ID in halo_ids:
    ID = int(ID)
    out_file = os.path.join(output_dir, f'{ID}.hdf5')
    if os.path.exists(out_file):
        print(f"Subhalo {ID}: output exists, skipped")
        continue
    print(f"Processing Subhalo ID: {ID}")

    try:
        # ---------- trajectory data ----------
        trajectory_data = load_galaxy_trajectory(orbit_file_path, ID)
        if trajectory_data is None:
            raise RuntimeError(f"Galaxy {ID}: not in the trajectory file")

        t_min_orbit = trajectory_data['t'].min()
        t_max_orbit = trajectory_data['t'].max()
        print(f"    trajectory time range: {t_min_orbit:.3f} - {t_max_orbit:.3f} Gyr")

        orbits = GalaxyPoseTrajectory(
            trajectory_data['t'],
            trajectory_data['cpos'],
            trajectory_data['cvel'],
            angular_momentum=trajectory_data['ang_mom'],
            box_size=35000/0.6774
        )

        # ---------- particle data (dm, star, gas, bh all loaded) ----------
        testsub = snap99.load_particle(ID)
        testsub.physical_units()

        n_total = len(testsub.s)
        if n_total == 0:
            print(f"    No stellar particles for {ID}, skip")
            continue

        # ==========================================================
        # [F] the following stellar arrays must be extracted before face_on
        # ==========================================================
        birth_pos  = testsub.s['BirthPos'].in_units('a kpc').view(np.ndarray).copy()
        birth_vel  = testsub.s['BirthVel'].in_units('kpc Gyr**-1 a**1/2').view(np.ndarray).copy()
        mass       = testsub.s['mass'].in_units('Msol').view(np.ndarray).copy()
        init_mass  = testsub.s['GFM_InitialMass'].in_units('Msol').view(np.ndarray).copy()
        age   = testsub.s['age'].in_units('Gyr').view(np.ndarray).copy()
        aform = testsub.s['aform'].view(np.ndarray).copy()
        tform = testsub.s['tform'].in_units('Gyr').view(np.ndarray).copy()
        iord  = testsub.s['iord'].view(np.ndarray).copy()
        metals_all = testsub.s['GFM_Metals'].view(np.ndarray).copy()

        X_H_all  = metals_all[:, 0]
        X_Mg_all = metals_all[:, 6]
        X_Fe_all = metals_all[:, 8]

        # ---------- remove wind particles ----------
        valid_mask = (
            (aform > 0) & (aform <= 1) &
            (tform > 0) & np.isfinite(tform) & np.isfinite(aform)
        )
        n_valid = int(np.sum(valid_mask))
        n_wind  = n_total - n_valid
        print(f"    stellar-family particles: {n_total}, real stars: {n_valid}, wind particles: {n_wind}")
        if n_valid == 0:
            print(f"    warning: {ID} has no valid stars, skipped")
            continue

        mass_valid      = mass[valid_mask]
        age_valid       = age[valid_mask]
        aform_valid     = aform[valid_mask]
        tform_valid     = tform[valid_mask]
        X_H_v  = X_H_all[valid_mask]
        X_Mg_v = X_Mg_all[valid_mask]
        X_Fe_v = X_Fe_all[valid_mask]

        # ---------- [I][J] abundances ----------
        feh_ok  = (X_H_v > 0)  & (X_Fe_v > 0)
        mgfe_ok = (X_Fe_v > 0) & (X_Mg_v > 0)
        FeH  = np.full(n_valid, np.nan)
        MgFe = np.full(n_valid, np.nan)
        FeH[feh_ok]   = np.log10((X_Fe_v[feh_ok] / A_Fe) / (X_H_v[feh_ok] / A_H)) - SUN_FEH_TERM
        MgFe[mgfe_ok] = np.log10((X_Mg_v[mgfe_ok] / A_Mg) / (X_Fe_v[mgfe_ok] / A_Fe)) - SUN_MGFE_TERM

        # ---------- [B] clamp tform at the upper boundary ----------
        edge_mask = (tform_valid > t_max_orbit) & (tform_valid < t_max_orbit + T_EDGE_TOL)
        n_edge = int(np.sum(edge_mask))
        tform_for_orbit = tform_valid.copy()
        tform_for_orbit[edge_mask] = t_max_orbit
        n_early = int(np.sum(tform_valid < t_min_orbit))
        print(f"    clamped at upper boundary: {n_edge}, born before trajectory start: {n_early}")

        # ---------- birth-coordinate transformation (galpos) ----------
        birth = pynbody_decorate.make_star_birth(
            orbits,
            birth_pos=birth_pos[valid_mask],
            mass=mass_valid,
            birth_velocity=birth_vel[valid_mask],
            birth_time=tform_for_orbit,
            birth_pos_units="a kpc",
            birth_velocity_units='kpc Gyr**-1 a**1/2',
            scale_factor=aform_valid
        )
        birth.align_with_galaxy(orientation_align=True)

        # [A] alignment assertions
        assert len(birth.s) == n_valid, f"{ID}: make_star_birth changed the number of particles"
        assert np.allclose(np.asarray(mass_valid, dtype=np.float64),
                           np.asarray(birth.s['mass'], dtype=np.float64),
                           rtol=1e-6), f"{ID}: make_star_birth may have reordered the particles"

        bpos = np.asarray(birth.s['pos'], dtype=np.float64)
        birth_r_2d = np.sqrt(bpos[:, 0]**2 + bpos[:, 1]**2)
        birth_r_3d = np.sqrt(np.sum(bpos**2, axis=1))
        has_birth = np.isfinite(birth_r_2d) & np.isfinite(birth_r_3d)   # [P]
        n_nan = int(np.sum(~has_birth))
        if n_nan != n_early:
            print(f"    *** warning *** number of NaN ({n_nan}) != number born before trajectory start ({n_early})")

        # ---------- z=0 face_on (the whole subhalo rotated together: dm/star/gas/bh) ----------
        center_pos, center_vel, rot = testsub.face_on(rmax=FACE_ON_RMAX,
                                                      alignwith=FACE_ON_ALIGN)

        s_pos = arr(testsub.s['pos'], 'kpc', np.float64)
        s_vel = arr(testsub.s['vel'], 'km s**-1', np.float64)

        # [L] young-star self-check: birth radius vs current radius
        cur_r_2d = np.sqrt(s_pos[valid_mask, 0]**2 + s_pos[valid_mask, 1]**2)
        for lo, hi in [(0, 0.02), (0.02, 0.05), (0.05, 0.1)]:
            m = (age_valid >= lo) & (age_valid < hi) & has_birth
            if np.sum(m) > 10:
                d = birth_r_2d[m] - cur_r_2d[m]
                print(f"    [self-check] age {lo:.2f}-{hi:.2f} Gyr ({np.sum(m):>5d} stars): "
                      f"median {np.median(d):+.3f} kpc, |dr|95% = "
                      f"{np.percentile(np.abs(d), 95):.3f} kpc")

        # ---------- collect all particle types ----------
        stars = {
            'iord':       iord[valid_mask].astype(np.uint64),
            'pos':        s_pos[valid_mask].astype(np.float32),
            'vel':        s_vel[valid_mask].astype(np.float32),
            'mass':       mass_valid.astype(np.float32),
            'init_mass':  init_mass[valid_mask].astype(np.float32),
            'age':        age_valid.astype(np.float32),
            'aform':      aform_valid.astype(np.float32),
            'FeH':        FeH.astype(np.float32),
            'MgFe':       MgFe.astype(np.float32),
            'X_H':        X_H_v.astype(np.float32),
            'X_Mg':       X_Mg_v.astype(np.float32),
            'X_Fe':       X_Fe_v.astype(np.float32),
            'birth_pos':  bpos.astype(np.float32),
            'birth_R_2D': birth_r_2d.astype(np.float32),
            'birth_R_3D': birth_r_3d.astype(np.float32),
            'has_birth':  has_birth,
        }
        wind = {
            'pos':  s_pos[~valid_mask].astype(np.float32).reshape(-1, 3),
            'mass': mass[~valid_mask].astype(np.float32),
        }

        others = {}
        for fam_name, fam_attr in [('dm', 'dm'), ('gas', 'g'), ('bh', 'bh')]:
            try:
                sub = getattr(testsub, fam_attr)
                n = len(sub)
            except Exception:
                n = 0
            if n > 0:
                others[fam_name] = {
                    'pos':  arr(sub['pos'], 'kpc'),
                    'mass': arr(sub['mass'], 'Msol'),
                }
            else:
                others[fam_name] = {
                    'pos':  np.zeros((0, 3), np.float32),
                    'mass': np.zeros(0, np.float32),
                }

        # ---------- write HDF5 (temporary file, then rename) ----------
        tmp_file = out_file + '.tmp'
        with h5py.File(tmp_file, 'w') as h5:
            hdr = h5.create_group('header')
            hdr.attrs['ID'] = ID
            hdr.attrs['snap'] = 99
            for key, getter in [('a', lambda: float(testsub.properties['a'])),
                                ('h', lambda: float(testsub.properties['h'])),
                                ('time_Gyr', lambda: float(testsub.properties['t'].in_units('Gyr')))]:
                try:
                    hdr.attrs[key] = getter()
                except Exception:
                    hdr.attrs[key] = np.nan
            hdr.attrs['face_on_rmax_kpc'] = FACE_ON_RMAX
            hdr.attrs['face_on_alignwith'] = FACE_ON_ALIGN
            hdr.attrs['units'] = 'pos: kpc (physical); vel: km/s; mass: Msol; age: Gyr'
            hdr.attrs['solar_scale'] = 'Asplund 2009: logeps_Fe=7.50, logeps_Mg=7.60'
            hdr.attrs['script'] = 'data_obtain.py'

            write_group(h5, 'frame', {
                'center_pos': arr(center_pos, None, np.float64),
                'center_vel': arr(center_vel, None, np.float64),
                'rotation':   np.asarray(rot, dtype=np.float64),
            })
            h5['frame'].attrs['center_pos_units'] = str(center_pos.units)
            h5['frame'].attrs['center_vel_units'] = str(center_vel.units)

            write_group(h5, 'stars', stars)
            write_group(h5, 'wind', wind)
            for fam_name, d in others.items():
                write_group(h5, fam_name, d)
        os.replace(tmp_file, out_file)

        size_mb = os.path.getsize(out_file) / 1e6
        ledger_all[ID] = {
            'star_family': n_total, 'wind': n_wind, 'valid': n_valid,
            'edge_clipped': n_edge, 'no_birth': n_nan,
            'dm': len(others['dm']['mass']), 'gas': len(others['gas']['mass']),
            'bh': len(others['bh']['mass']), 'MB': size_mb,
        }
        print(f"    saved {out_file} ({size_mb:.1f} MB)")
        print(f"    dm {ledger_all[ID]['dm']}, gas {ledger_all[ID]['gas']}, "
              f"bh {ledger_all[ID]['bh']}, stars {n_valid} (without birth radius {n_nan})")
        print()

        del testsub, birth

    except Exception as e:
        print(f"*** ERROR *** Subhalo {ID}: {e}")
        processing_errors.append((ID, str(e)))
        tmp_file = out_file + '.tmp'
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
        continue

# =================== final report ===================
print("\n" + "=" * 90)
if processing_errors:
    print(f"*** {len(processing_errors)} galaxies failed:")
    for gid, err in processing_errors:
        print(f"  - {gid}: {err}")
print(f"{'ID':>8} {'star fam.':>9} {'wind':>6} {'stars':>9} {'noRbir':>7} "
      f"{'DM':>9} {'gas':>9} {'bh':>4} {'MB':>7}")
for gid, L in ledger_all.items():
    print(f"{gid:>8} {L['star_family']:>9} {L['wind']:>6} {L['valid']:>9} "
          f"{L['no_birth']:>7} {L['dm']:>9} {L['gas']:>9} {L['bh']:>4} {L['MB']:>7.1f}")
print(f"\nsucceeded: {len(ledger_all)} / {len(halo_ids)}, output directory: {output_dir}")
print("=" * 90)
