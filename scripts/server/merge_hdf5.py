# -*- coding: utf-8 -*-
"""
merge_hdf5.py — merge the per-galaxy files data_obtain/<ID>.hdf5 into a single file for download.
Run on the server after data_obtain.py has finished:  python3 merge_hdf5.py

Output: data_obtain_all.hdf5
  /<ID>/header, /<ID>/frame, /<ID>/stars, /<ID>/wind, /<ID>/dm, /<ID>/gas, /<ID>/bh
  The structure inside each galaxy group is identical to the single files (dataset compression settings are kept).
"""
import glob
import os
import h5py

in_dir = './data_obtain/'
out_file = 'data_obtain_all.hdf5'

files = sorted(glob.glob(os.path.join(in_dir, '*.hdf5')))
if not files:
    raise SystemExit(f"no .hdf5 files in {in_dir}; run data_obtain.py first")

tmp = out_file + '.tmp'
with h5py.File(tmp, 'w') as fout:
    for fn in files:
        gid = os.path.splitext(os.path.basename(fn))[0]
        with h5py.File(fn, 'r') as fin:
            grp = fout.create_group(gid)
            for key in fin.keys():
                fin.copy(fin[key], grp, name=key)
        print(f"  merged {gid}")
os.replace(tmp, out_file)
print(f"done: {len(files)} galaxies -> {out_file} ({os.path.getsize(out_file)/1e9:.2f} GB)")
