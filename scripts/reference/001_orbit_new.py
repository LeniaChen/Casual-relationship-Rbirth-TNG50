# -*- coding: utf-8 -*-
"""
001_orbit data extraction script — final version (class-1 + class-2 changes + Current_z column)

Class-1 changes (logic / alignment, verified):
  [A] alignment assertion after make_star_birth (verified on 487742: length and order preserved)
  [B] floating-point clamp of tform at the upper boundary (stars born at z=0 were misjudged because of rounding)
  [C] missing trajectory data raises directly; fallback branch removed (trajectories confirmed complete for the 13 galaxies)
  [D] output the real particle ID (iord)
  [E] empty-galaxy check moved after load_particle
  [F] BirthPos/BirthVel must be extracted before face_on (order locked)
  [G] removed unused gas_fields / reading of trajectory align, a, re
  [H] per-galaxy particle ledger

Class-2 changes (numerical correctness):
  [I] abundance formulas: number-density ratios + Asplund 2009 logε solar reference
      [Fe/H]  = log10( (X_Fe/A_Fe) / (X_H /A_H ) ) + 4.50
      [Mg/Fe] = log10( (X_Mg/A_Mg) / (X_Fe/A_Fe) ) - 0.10
  [J] FeH / MgFe masks decoupled
  [K] InitMass column added to the output
  [L] young-star self-check (binned by age): end-to-end validation of the coordinate chain

  [M] Current_z column added (vertical coordinate in the z=0 face_on disk frame, physical kpc).
      All valid stars are written, without a |z| cut — disk-star selection (e.g. |z|<0.3 kpc)
      is left to downstream scripts, so changing the threshold does not require rerunning the expensive extraction.

Verified (from the code that generated the trajectory file):
  cpos = comoving kpc without h (in_units("a kpc") after physical_units), consistent with box_size=35000/0.6774
  t has the same origin as tform in the main script (same pynbody cosmology conversion)
  ang_mom is a direction vector in the box frame, the same frame as BirthPos
Notes (do not affect this script's output):
  cvel is stored as "a kpc Gyr**-1", unlike the "kpc Gyr**-1 a**1/2" convention of BirthVel;
  only radii/coordinates are written, no velocities. If birth velocities are used in the future,
  first check which units GalaxyPoseTrajectory expects for cvel.
  At early snapshots of the trajectory, alignwith switches to baryon/all, so the "disk plane" of the oldest stars
  is defined differently from later times; Birth_R_3D is recommended for early-time analyses.
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
# [G] this script does not use gas; gas_fields are not loaded
snap99.load_particle_para['star_fields'] = [
    'GFM_StellarFormationTime', 'BirthPos', 'BirthVel', 'GFM_InitialMass',
    'GFM_Metals'
]

orbit_file_path = 'TNG50_disk_trajectory.hdf5'

halo_ids = np.array([487742, 530330, 538905, 545437, 546870,
                     547293, 552581, 557721, 560421, 562029,
                     567607, 569251, 607654])

print(f"processing {len(halo_ids)} galaxies: {halo_ids}")

output_dir = './001_orbit_new/'
os.makedirs(output_dir, exist_ok=True)

# ==================== [I] abundance constants ====================
# atomic masses
A_H, A_Mg, A_Fe = 1.008, 24.305, 55.845
# Asplund et al. 2009 photospheric abundances logε (scale logε_H = 12.00)
LOGEPS_MG_SUN = 7.60
LOGEPS_FE_SUN = 7.50
SUN_FEH_TERM  = LOGEPS_FE_SUN - 12.0            # = -4.50 → +4.50 in the formula
SUN_MGFE_TERM = LOGEPS_MG_SUN - LOGEPS_FE_SUN   # = +0.10 → -0.10 in the formula

# [B] tolerance of the tform upper-boundary clamp (Gyr): much larger than rounding (~1e-6), much smaller than physical time scales
T_EDGE_TOL = 1e-3

processing_errors = []
ledger_all = {}   # [H] particle ledger

def load_galaxy_trajectory(orbit_file, galaxy_id):
    """[G] read only the fields actually used: t, cpos, cvel, ang_mom"""
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

# ==================== main loop ====================
for ID in halo_ids:
    ID = int(ID)
    print(f"Processing Subhalo ID: {ID}")

    try:
        # ---------- trajectory data ----------
        trajectory_data = load_galaxy_trajectory(orbit_file_path, ID)
        # [C] missing trajectory = hard error
        if trajectory_data is None:
            raise RuntimeError(
                f"Galaxy {ID}: not in the trajectory file. Trajectories of the 13 target galaxies were confirmed complete; "
                f"this error means the trajectory file or the ID list was changed, please check.")

        t_min_orbit = trajectory_data['t'].min()
        t_max_orbit = trajectory_data['t'].max()
        print(f"    trajectory time range: {t_min_orbit:.3f} - {t_max_orbit:.3f} Gyr")

        orbits = GalaxyPoseTrajectory(
            trajectory_data['t'],
            trajectory_data['cpos'],
            trajectory_data['cvel'],
            angular_momentum=trajectory_data['ang_mom'],
            box_size=35000/0.6774   # verified: comoving kpc without h, same as cpos in the trajectory file
        )

        # ---------- particle data ----------
        testsub = snap99.load_particle(ID)
        testsub.physical_units()

        # [E] check for empty galaxies right away
        n_total = len(testsub.s['mass'])
        if n_total == 0:
            print(f"    No stellar particles found for Subhalo ID: {ID}, skip")
            continue

        # ==========================================================
        # [F] !!! extraction order locked: all arrays below must be extracted before face_on !!!
        # pynbody's face_on rotates all 3-vector arrays of the SimSnap in place
        # (including BirthPos/BirthVel). The coordinate transformation of BirthPos/BirthVel is done by
        # galpos (using the centre and disk plane at each star's birth time) and must never be
        # contaminated first by face_on with the z=0 centre/disk plane.
        # ==========================================================
        birth_pos  = testsub.s['BirthPos'].in_units('a kpc').view(np.ndarray).copy()
        birth_vel  = testsub.s['BirthVel'].in_units('kpc Gyr**-1 a**1/2').view(np.ndarray).copy()
        mass       = testsub.s['mass'].in_units('Msol').view(np.ndarray).copy()
        init_mass  = testsub.s['GFM_InitialMass'].in_units('Msol').view(np.ndarray).copy()  # [K]
        age   = testsub.s['age'].in_units('Gyr').view(np.ndarray).copy()
        aform = testsub.s['aform'].view(np.ndarray).copy()
        tform = testsub.s['tform'].in_units('Gyr').view(np.ndarray).copy()
        iord  = testsub.s['iord'].view(np.ndarray).copy()      # [D] real particle IDs
        metals_all = testsub.s['GFM_Metals'].view(np.ndarray).copy()  # (N, 10)

        # order of TNG GFM_Metals: [H, He, C, N, O, Ne, Mg, Si, Fe, other]
        X_H_all  = metals_all[:, 0]
        X_Mg_all = metals_all[:, 6]
        X_Fe_all = metals_all[:, 8]

        # ---------- selection: remove wind particles (aform<0 are wind particles, TNG convention) ----------
        valid_mask = (
            (aform > 0) & (aform <= 1) &
            (tform > 0) & np.isfinite(tform) & np.isfinite(aform)
        )
        n_valid = int(np.sum(valid_mask))
        n_wind  = n_total - n_valid
        print(f"    star particles before filtering: {n_total}")
        print(f"    star particles after filtering: {n_valid} (removed wind etc. {n_wind})")

        if n_valid == 0:
            print(f"    warning: Galaxy ID {ID} has no valid star particles, skipped")
            continue

        birth_pos_valid = birth_pos[valid_mask]
        birth_vel_valid = birth_vel[valid_mask]
        mass_valid      = mass[valid_mask]
        init_mass_valid = init_mass[valid_mask]     # [K]
        age_valid   = age[valid_mask]
        aform_valid = aform[valid_mask]
        tform_valid = tform[valid_mask]
        iord_valid  = iord[valid_mask]

        X_H_v  = X_H_all[valid_mask]
        X_Mg_v = X_Mg_all[valid_mask]
        X_Fe_v = X_Fe_all[valid_mask]

        # ---------- [I][J] abundances (number-density ratios + Asplund 2009) ----------
        feh_ok  = (X_H_v > 0)  & (X_Fe_v > 0)
        mgfe_ok = (X_Fe_v > 0) & (X_Mg_v > 0)

        FeH  = np.full(n_valid, np.nan)
        MgFe = np.full(n_valid, np.nan)
        FeH[feh_ok]   = np.log10((X_Fe_v[feh_ok]  / A_Fe) / (X_H_v[feh_ok]   / A_H)) \
                        - SUN_FEH_TERM                     # i.e. + 4.50
        MgFe[mgfe_ok] = np.log10((X_Mg_v[mgfe_ok] / A_Mg) / (X_Fe_v[mgfe_ok] / A_Fe)) \
                        - SUN_MGFE_TERM                    # i.e. - 0.10

        # ---------- [B] floating-point clamp of tform at the upper boundary ----------
        # stars born at z=0 (aform=1) have a tform physically equal to the end of the trajectory,
        # but the two conversions can differ by ~1e-6 Gyr through rounding → the interpolator wrongly returns NaN (out of range).
        # Clamping back to t_max removes the rounding difference; it is not an extrapolation.
        # The lower boundary (before the start of the trajectory) is not clamped: those stars are truly out of range and stay NaN.
        edge_mask = (tform_valid > t_max_orbit) & (tform_valid < t_max_orbit + T_EDGE_TOL)
        n_edge = int(np.sum(edge_mask))
        tform_for_orbit = tform_valid.copy()
        tform_for_orbit[edge_mask] = t_max_orbit

        n_early = int(np.sum(tform_valid < t_min_orbit))
        print(f"    clamped at upper boundary (born at z=0, rounding): {n_edge}")
        print(f"    born before trajectory start (removed as NaN):   {n_early}")

        # ---------- birth-coordinate transformation (galpos) ----------
        birth = pynbody_decorate.make_star_birth(
            orbits,
            birth_pos=birth_pos_valid,
            mass=mass_valid,
            birth_velocity=birth_vel_valid,
            birth_time=tform_for_orbit,          # [B] clamped times
            birth_pos_units="a kpc",
            birth_velocity_units='kpc Gyr**-1 a**1/2',
            scale_factor=aform_valid
        )
        birth.align_with_galaxy(orientation_align=True)

        # [A] alignment safeguard (verified; if a galpos update changes the behaviour it fails here immediately)
        assert len(birth.s) == n_valid, (
            f"Galaxy {ID}: make_star_birth changed the number of particles "
            f"{n_valid} -> {len(birth.s)}, row matching broken!")
        assert np.allclose(np.asarray(mass_valid, dtype=np.float64),
                           np.asarray(birth.s['mass'], dtype=np.float64),
                           rtol=1e-6), (
            f"Galaxy {ID}: masses returned by make_star_birth differ element-wise from the input; "
            f"particles may have been reordered, row matching broken!")

        # ---------- current coordinates (extracted after face_on) ----------
        testsub.face_on(rmax=8, alignwith="star")   # rotate in place to the z=0 disk frame
        _pos = np.array(testsub.s['pos'])
        current_r_2d = np.sqrt(_pos[valid_mask, 0]**2 + _pos[valid_mask, 1]**2)
        current_r_3d = np.sqrt(np.sum(_pos[valid_mask]**2, axis=1))
        current_z    = _pos[valid_mask, 2]          # [M] vertical coordinate in the z=0 disk frame (physical kpc)

        # ---------- birth radius ----------
        _bpos = np.asarray(birth.s['pos'], dtype=np.float64)
        birth_r_2d = np.sqrt(_bpos[:, 0]**2 + _bpos[:, 1]**2)
        birth_r_3d = np.sqrt(np.sum(_bpos**2, axis=1))

        # ---------- remove NaN (should only be particles born before the trajectory start) ----------
        finite_mask = np.isfinite(birth_r_2d) & np.isfinite(birth_r_3d)
        n_nan = int(np.sum(~finite_mask))
        n_final = int(np.sum(finite_mask))
        print(f"    NaN removed: {n_nan}, remaining: {n_final}")

        # [B] consistency check: NaN may only come from particles born before the trajectory start
        if n_nan != n_early:
            print(f"    *** warning *** number of NaN ({n_nan}) != number born before trajectory start ({n_early}), "
                  f"unexplained source of NaN, please check!")

        # [L] self-check: dr in age bins (real migration → grows with age; coordinate error → independent of age)
        for lo, hi in [(0, 0.02), (0.02, 0.05), (0.05, 0.1)]:
            m = (age_valid >= lo) & (age_valid < hi) & finite_mask
            if np.sum(m) > 10:
                d = birth_r_2d[m] - current_r_2d[m]
                print(f"    [self-check] age {lo:.2f}-{hi:.2f} Gyr ({np.sum(m):>5d} stars): "
                      f"median {np.median(d):+.3f} kpc, "
                      f"|dr|95% = {np.percentile(np.abs(d), 95):.3f} kpc")

        # ---------- [H] particle ledger ----------
        ledger_all[ID] = {
            'total': n_total, 'wind_removed': n_wind, 'valid': n_valid,
            'edge_clipped': n_edge, 'early_nan_removed': n_nan, 'final': n_final,
        }

        # ---------- save (all valid stars; disk selection left to downstream scripts) ----------
        output_file = os.path.join(output_dir, f'{ID}.txt')
        with open(output_file, 'w') as f:
            # [D] ParticleID = real TNG ParticleIDs (iord)
            # [K] Mass = current mass at z=0; InitMass = GFM_InitialMass (mass at formation)
            # [M] Current_z = vertical coordinate in the z=0 face_on disk frame (physical kpc, signed)
            #     example disk selection (downstream): |Current_z| < 0.3
            f.write(f"{'ParticleID':>15} {'Age_Gyr':>10} {'Current_R_2D':>14} {'Birth_R_2D':>12} {'Current_R_3D':>14} {'Birth_R_3D':>12} {'Current_z':>11} {'Mass':>12} {'InitMass':>12} {'FeH':>10} {'MgFe':>10}\n")
            for i in np.where(finite_mask)[0]:
                feh_str  = f"{FeH[i]:>10.6f}"  if np.isfinite(FeH[i])  else f"{'nan':>10}"
                mgfe_str = f"{MgFe[i]:>10.6f}" if np.isfinite(MgFe[i]) else f"{'nan':>10}"
                f.write(f"{iord_valid[i]:>15d} {age_valid[i]:>10.6f} {current_r_2d[i]:>14.6f} {birth_r_2d[i]:>12.6f} {current_r_3d[i]:>14.6f} {birth_r_3d[i]:>12.6f} {current_z[i]:>11.6f} {mass_valid[i]:>12.6f} {init_mass_valid[i]:>12.6f} {feh_str} {mgfe_str}\n")

        print(f"Subhalo {ID}: {n_final} stellar particles saved")
        print(f"  ledger: total {n_total} → wind removed {n_wind} → valid {n_valid} "
              f"→ clamped {n_edge} → NaN removed {n_nan} → final {n_final}")
        print(f"  Age range:          {np.nanmin(age_valid[finite_mask]):.2f} - {np.nanmax(age_valid[finite_mask]):.2f} Gyr")
        print(f"  Current_R_2D range: {np.nanmin(current_r_2d[finite_mask]):.2f} - {np.nanmax(current_r_2d[finite_mask]):.2f} kpc")
        print(f"  Birth_R_2D range:   {np.nanmin(birth_r_2d[finite_mask]):.2f} - {np.nanmax(birth_r_2d[finite_mask]):.2f} kpc")
        print(f"  Current_z range:    {np.nanmin(current_z[finite_mask]):+.3f} ~ {np.nanmax(current_z[finite_mask]):+.3f} kpc")
        n_disk_est = int(np.sum(finite_mask & (np.abs(current_z) < 0.3)))
        print(f"  particles with |z|<0.3 kpc (for reference): {n_disk_est} / {n_final}")
        print(f"  Mass range:         {np.nanmin(mass_valid[finite_mask]):.3f} - {np.nanmax(mass_valid[finite_mask]):.3f} Msol")
        print(f"  InitMass range:     {np.nanmin(init_mass_valid[finite_mask]):.3f} - {np.nanmax(init_mass_valid[finite_mask]):.3f} Msol")
        feh_fin  = FeH[finite_mask]
        mgfe_fin = MgFe[finite_mask]
        print(f"  [Fe/H]:  finite {np.sum(np.isfinite(feh_fin))} / {n_final}, "
              f"median = {np.nanmedian(feh_fin):+.3f}, "
              f"range {np.nanmin(feh_fin):+.2f} ~ {np.nanmax(feh_fin):+.2f}")
        print(f"  [Mg/Fe]: finite {np.sum(np.isfinite(mgfe_fin))} / {n_final}, "
              f"median = {np.nanmedian(mgfe_fin):+.3f}, "
              f"range {np.nanmin(mgfe_fin):+.2f} ~ {np.nanmax(mgfe_fin):+.2f}")
        print()

    except Exception as e:
        error_info = f"Error processing Subhalo ID {ID}: {str(e)}"
        print(f"*** ERROR *** {error_info}")
        processing_errors.append((ID, str(e)))
        continue

# =================== final report ===================
print("\n" + "="*80)
print("PROCESSING COMPLETED - FINAL REPORT")
print("="*80)

if processing_errors:
    print(f"*** ERROR *** {len(processing_errors)} galaxies had PROCESSING ERRORS:")
    for gal_id, error in processing_errors:
        print(f"  - Galaxy ID {gal_id}: {error}")
    print()

print(f"{'ID':>8} {'total':>9} {'wind':>6} {'valid':>9} {'clamp':>5} {'NaN':>6} {'final':>9}")
for gid, L in ledger_all.items():
    print(f"{gid:>8} {L['total']:>9} {L['wind_removed']:>6} {L['valid']:>9} "
          f"{L['edge_clipped']:>5} {L['early_nan_removed']:>6} {L['final']:>9}")

print(f"\nSUMMARY:")
print(f"  Total galaxies: {len(halo_ids)}")
print(f"  Successfully processed: {len(ledger_all)}")
print(f"\nData saved in directory: {output_dir}")
print("="*80)