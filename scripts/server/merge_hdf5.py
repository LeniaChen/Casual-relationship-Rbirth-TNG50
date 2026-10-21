# -*- coding: utf-8 -*-
"""
merge_hdf5.py —— 把 data_obtain/ 里每个星系的 <ID>.hdf5 合并成一个文件, 方便下载。
在服务器上 data_obtain.py 跑完后运行:  python3 merge_hdf5.py

输出: data_obtain_all.hdf5
  /<ID>/header, /<ID>/frame, /<ID>/stars, /<ID>/wind, /<ID>/dm, /<ID>/gas, /<ID>/bh
  每个星系组内部的结构和单个文件完全相同 (数据集的压缩设置也原样保留)。
"""
import glob
import os
import h5py

in_dir = './data_obtain/'
out_file = 'data_obtain_all.hdf5'

files = sorted(glob.glob(os.path.join(in_dir, '*.hdf5')))
if not files:
    raise SystemExit(f"{in_dir} 里没有 .hdf5 文件, 先运行 data_obtain.py")

tmp = out_file + '.tmp'
with h5py.File(tmp, 'w') as fout:
    for fn in files:
        gid = os.path.splitext(os.path.basename(fn))[0]
        with h5py.File(fn, 'r') as fin:
            grp = fout.create_group(gid)
            for key in fin.keys():
                fin.copy(fin[key], grp, name=key)
        print(f"  已合并 {gid}")
os.replace(tmp, out_file)
print(f"完成: {len(files)} 个星系 -> {out_file} ({os.path.getsize(out_file)/1e9:.2f} GB)")
