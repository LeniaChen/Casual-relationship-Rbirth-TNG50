from AnastrisTNG import TNGsimulation, TNGtools, Anatools, TNGgroupcat
import numpy as np
import os
import h5py
from galpos import GalaxyPoseTrajectory
from galpos import pynbody_decorate

# 设置路径 - 适配TNG服务器
run = 'TNG50'
BasePath = '/home/tnguser/sims.TNG/' + run + '-1/output'
snap99 = TNGsimulation.Snapshot(BasePath, 99)
snap99.load_particle_para['star_fields'] = ['GFM_StellarFormationTime', 'BirthPos', 'BirthVel', 'GFM_InitialMass']
snap99.load_particle_para['gas_fields'] = ['ElectronAbundance', 'GFM_CoolingRate', 'GFM_Metallicity', 'GFM_Metals', 'InternalEnergy', 'Machnumber', 'StarFormationRate']

# 轨道数据文件路径 - 需要根据服务器实际位置调整
orbit_file_path = 'TNG50_disk_trajectory.hdf5'  # 假设文件在当前工作目录

# ============ 配置HaloID列表 ============
# 方式1：直接在代码中指定要处理的HaloID（推荐用于服务器环境）
halo_ids = np.array([333424,348901,355724,368843,371126,377655,383975,386269,388544,394621,396628,402555,411449,413372,414917,425719,427211,429471,430864,432106,433289,435752,436932,439099,447914,448830,452031,452978,454171,455291,456326,458470,460722,461785,462710,466549,467415,470345,471248,471996,473329,475016,476266,477328,478216,479290,479938,480802,482155,483594,484448,485056,487742,488530,490814,491426,492244,492876,494011,494709,496788,497557,498522,499704,500576,501208,501725,502371,502995,503437,503987,505586,506720,507292,508538,509091,510273,510585,511303,513845,514272,514829,515695,516760,517271,517899,518682,519311,519720,520310,520885,521429,521803,522983,523889,526029,526478,527309,528322,528836,529365,529855,530330,530852,531320,531910,532760,533060,533590,533998,535050,535410,536146,537236,537488,537941,538370,538905,539333,539667,540082,540452,541218,541497,541847,542252,542662,543376,543729,544408,544697,545437,546114,546474,546870,547293,547844,548151,548418,548790,549090,549516,549746,550149,550475,551190,551541,551973,552414,552581,552879,553230,553489,553837,554189,554523,554798,555013,555287,555601,555815,556304,557076,557396,557721,558067,558410,558731,559036,559711,560082,560421,560751,561001,561325,561676,562029,562742,563471,563732,563943,564268,564825,565089,565947,566365,566661,567004,567382,567607,568923,569251,569452,569653,569877,570319,570598,570841,571072,571633,571908,572121,572328,572840,573062,573413,573738,574037,574286,574615,575190,575356,575596,576062,576254,576516,576705,577125,577372,577576,577719,577924,578294,578541,578834,579232,579509,579759,580035,580250,580406,580641,581308,581656,581912,582137,582344,582612,582809,583027,583256,583745,584286,584511,584875,585282,585839,586066,586169,586420,586670,587019,587256,587819,588004,588178,588647,588831,589496,589655,589873,590218,590385,590792,590926,591210,591469,591641,591796,592282,592437,592662,592840,593259,593480,593694,593904,594098,594731,594887,595100,595246,595442,595661,595887,596244,596660,596866,596986,597141,597530,597792,597934,598112,598320,598696,598840,598988,599465,599865,600096,600363,600555,600946,601123,601522,601693,601861,601997,602330,602509,602653,602820,603185,604241,604371,604516,604667,605054,605708,605844,606405,606820,606907,607310,607530,607712,608061,608386,608745,609001,610017,610361,610532,610819,610988,611152,611385,611780,612241,612699,612825,613006,613808,614180,614326,614489,615182,615292,615462,615607,615792,615903,616200,616526,617027,617149,617324,617559,617719,618120,619081,619223,619381,620348,620463,620832,621043,621167,621255,621372,621707,621855,622298,622437,623042,623122,623575,623877,625041,625185,625281,625423,625715,625800,625925,626198,626697,627459,627572,628949,629266,629409,629891,630760,630870,630975,631253,631558,631677,632310,633322,633437,633614,633963,634631,634707,635905,636201,636640,636729,636848,637631,637971,639432,642220,642303,643642,646733,647190,652267,696559,697948])  # 可以根据需要修改这个列表

# 方式2：如果有数据文件，可以尝试读取
# file_path = '002_data.txt'
# try:
#     data = np.loadtxt(file_path, dtype=str, skiprows=1)
#     if data.ndim == 1:
#         halo_ids = data.astype(int)
#     else:
#         halo_ids = data[:, 0].astype(int)
#     print(f"从文件 {file_path} 成功读取 {len(halo_ids)} 个HaloID")
# except Exception as e:
#     print(f"无法读取文件 {file_path}: {e}")
#     print("使用默认的HaloID列表")
#     halo_ids = np.array([333424, 348901, 355724, 368843, 371126])

print(f"将处理 {len(halo_ids)} 个星系: {halo_ids}")

# 输出目录 - 适配服务器环境
output_dir = './001_orbit/'  # 在当前目录下创建输出文件夹
os.makedirs(output_dir, exist_ok=True)

# 用于记录错误和回退情况的列表
missing_orbit_data = []
processing_errors = []

def load_galaxy_trajectory(orbit_file, galaxy_id):
    """从HDF5文件加载星系轨道数据"""
    try:
        with h5py.File(orbit_file, "r") as file:
            if str(galaxy_id) not in file:
                return None
                
            galaxy_group = file[str(galaxy_id)]
            snap_all = list(galaxy_group.keys())
            snap_all = np.sort(np.int16(snap_all))
            
            # 提取轨道数据
            align = np.array([galaxy_group[f"{i}"]['align'][...] for i in snap_all])
            a = np.array([galaxy_group[f"{i}"]['a'][...] for i in snap_all])
            ang_mom = np.array([galaxy_group[f"{i}"]['ang_mom'][...] for i in snap_all])
            cpos = np.array([galaxy_group[f"{i}"]['cpos'][...] for i in snap_all])
            cvel = np.array([galaxy_group[f"{i}"]['cvel'][...] for i in snap_all])
            t = np.array([galaxy_group[f"{i}"]['t'][...] for i in snap_all])
            re = np.array([galaxy_group[f"{i}"]['re'][...] for i in snap_all])
            
            return {
                'snap_all': snap_all,
                'align': align,
                'a': a,
                'ang_mom': ang_mom,
                'cpos': cpos,
                'cvel': cvel,
                't': t,
                're': re
            }
    except Exception as e:
        print(f"Error loading trajectory for galaxy {galaxy_id}: {str(e)}")
        return None

# 遍历每个ID处理数据
for ID in halo_ids:
    ID = int(ID)  # 确保ID是整数
    print(f"Processing Subhalo ID: {ID}")

    try:
        # 1. 加载轨道数据
        trajectory_data = load_galaxy_trajectory(orbit_file_path, ID)
        
        if trajectory_data is None:
            # =================== 强制报错和记录 ===================
            error_msg = f"*** CRITICAL WARNING *** Galaxy ID {ID}: NO ORBIT DATA FOUND! Using fallback original method!"
            print("="*80)
            print(error_msg)
            print("="*80)
            missing_orbit_data.append(ID)
            
            # 回退到原始方法
            testsub = snap99.load_particle(ID)
            testsub.physical_units()
            trans = testsub.face_on(alignwith='star', rmax=8, retmatrix=True)
            starbir = TNGtools.Star_birth(snap99, ID, usebirthvel=True, usebirthmass=True, useCM=False)
            starbir.transform(trans)
            starbir.physical_units()
            
            print(f"*** FALLBACK COMPLETED for Galaxy ID {ID} ***")
            print("-"*80)
            
        else:
            print(f"✓ Orbit data found for Galaxy ID {ID}, using improved method")
            
            # 2. 构建轨道追踪模型
            # TNG50盒子大小：35000 "a kpc/h" = 35000/0.6774 "a kpc"
            orbits = GalaxyPoseTrajectory(
                trajectory_data['t'], 
                trajectory_data['cpos'], 
                trajectory_data['cvel'], 
                angular_momentum=trajectory_data['ang_mom'], 
                box_size=35000/0.6774
            )
            
            # 3. 加载当前快照的粒子数据
            testsub = snap99.load_particle(ID)
            testsub.physical_units()
            
            # 4. 提取恒星出生属性
            birth_pos = testsub.s['BirthPos'].in_units('a kpc').view(np.ndarray)
            birth_vel = testsub.s['BirthVel'].in_units('kpc Gyr**-1 a**1/2').view(np.ndarray)
            birth_mass = testsub.s['GFM_InitialMass'].in_units('Msol').view(np.ndarray)
            mass = testsub.s['mass'].in_units('Msol').view(np.ndarray)
            age = testsub.s['age'].in_units('Gyr').view(np.ndarray)
            aform = testsub.s['aform'].view(np.ndarray)
            tform = testsub.s['tform'].in_units('Gyr').view(np.ndarray)
            
            # 过滤有效的恒星粒子（避免无效的尺度因子）
            valid_mask = (aform > 0) & (aform <= 1) & (tform > 0) & np.isfinite(tform) & np.isfinite(aform)
            print(f"    过滤前恒星粒子数: {len(mass)}")
            print(f"    过滤后恒星粒子数: {np.sum(valid_mask)}")
            
            if np.sum(valid_mask) == 0:
                print(f"    警告：Galaxy ID {ID} 没有有效的恒星粒子，跳过")
                continue
            
            # 只使用有效的粒子
            birth_pos_valid = birth_pos[valid_mask]
            birth_vel_valid = birth_vel[valid_mask]
            birth_mass_valid = birth_mass[valid_mask]
            mass_valid = mass[valid_mask]
            age_valid = age[valid_mask]
            aform_valid = aform[valid_mask]
            tform_valid = tform[valid_mask]
            
            # 5. 创建动态校正的恒星出生轨迹
            birth = pynbody_decorate.make_star_birth(
                orbits, 
                birth_pos=birth_pos_valid, 
                mass=mass_valid,
                birth_velocity=birth_vel_valid, 
                birth_time=tform_valid,
                birth_pos_units="a kpc", 
                birth_velocity_units='kpc Gyr**-1 a**1/2', 
                scale_factor=aform_valid
            )
            
            # 6. 应用坐标系对齐
            birth.align_with_galaxy(orientation_align=True)
            
            # 7. 设置当前观测视角
            tr = testsub.face_on(rmax=8, alignwith="star")
            
            # 创建starbir兼容对象用于后续处理
            starbir = birth
            
            print(f"✓ Dynamic coordinate correction completed for Galaxy ID {ID}")

        # 检查是否有恒星数据
        if len(testsub.s['mass']) == 0:
            print(f"No stellar particles found for Subhalo ID: {ID}")
            continue
            
        # 提取恒星数据 - 需要处理动态校正后的数据
        if trajectory_data is not None:
            # 使用改进方法时，需要匹配过滤后的粒子
            # 恒星ID - 使用有效粒子的索引
            particle_ids = np.arange(len(mass_valid))
            
            # 恒星年龄 - 使用过滤后的数据
            stellar_age = age_valid
            
            # 恒星质量 - 使用过滤后的质量数据
            stellar_mass = mass_valid
            
            # 现在的半径 - 使用过滤后的当前位置数据
            current_pos = testsub.s['pos'][valid_mask]
            current_r = np.sqrt(np.sum(current_pos**2, axis=1))
        else:
            # 回退方法时使用原始方式
            # 恒星ID
            particle_ids = np.arange(len(testsub.s['mass']))
            
            # 恒星年龄 - 按原始代码方式
            stellar_age = testsub.s['age'].in_units('Gyr')
            
            # 恒星质量 - 按原始代码方式
            stellar_mass = testsub.s['mass'].in_units('Msol')
            
            # 现在的半径 - 按原始代码方式
            current_r = testsub.s['R']
        
        # 出生时的半径 - 使用动态校正后的结果或原始结果
        if hasattr(starbir, 's') and 'R' in starbir.s:
            birth_r = starbir.s['R']
        else:
            # 如果starbir结构不同，计算径向距离
            if hasattr(starbir, 's') and 'pos' in starbir.s:
                birth_pos_corrected = starbir.s['pos']
                birth_r = np.sqrt(np.sum(birth_pos_corrected**2, axis=1))
            elif hasattr(starbir, 's'):
                # 回退情况，使用原始方法的结果
                birth_r = starbir.s['R']
            else:
                raise ValueError(f"Cannot extract birth radius for Galaxy ID {ID}")
        
        # 保存为txt文件，使用空格分隔和右对齐 - 增加质量列
        output_file = os.path.join(output_dir, f'{ID}.txt')
        
        with open(output_file, 'w') as f:
            # 写入表头 - 增加Mass列
            f.write(f"{'ParticleID':>12} {'Age_Gyr':>10} {'Current_R':>12} {'Birth_R':>10} {'Mass':>12}\n")
            
            # 写入数据行 - 增加质量数据
            for i in range(len(particle_ids)):
                f.write(f"{particle_ids[i]:>12} {stellar_age[i]:>10.6f} {current_r[i]:>12.6f} {birth_r[i]:>10.6f} {stellar_mass[i]:>12.6f}\n")
        
        # 打印基本统计信息 - 增加质量统计
        print(f"Subhalo {ID}: {len(particle_ids)} stellar particles saved")
        print(f"  Age range: {stellar_age.min():.2f} - {stellar_age.max():.2f} Gyr")
        print(f"  Current radial range: {current_r.min():.2f} - {current_r.max():.2f}")
        print(f"  Birth radial range: {birth_r.min():.2f} - {birth_r.max():.2f}")
        print(f"  Mass range: {stellar_mass.min():.3f} - {stellar_mass.max():.3f} Msol")
        print()
        
    except Exception as e:
        # 捕获任何异常，打印错误信息，然后继续下一个ID
        error_info = f"Error processing Subhalo ID {ID}: {str(e)}"
        print(f"*** ERROR *** {error_info}")
        processing_errors.append((ID, str(e)))
        continue

# =================== 最终报告 ===================
print("\n" + "="*80)
print("PROCESSING COMPLETED - FINAL REPORT")
print("="*80)

if missing_orbit_data:
    print(f"*** WARNING *** {len(missing_orbit_data)} galaxies had NO ORBIT DATA:")
    print("These galaxies used FALLBACK original method:")
    for gal_id in missing_orbit_data:
        print(f"  - Galaxy ID: {gal_id}")
    print()

if processing_errors:
    print(f"*** ERROR *** {len(processing_errors)} galaxies had PROCESSING ERRORS:")
    for gal_id, error in processing_errors:
        print(f"  - Galaxy ID {gal_id}: {error}")
    print()

total_processed = len(halo_ids) - len(processing_errors)
successful_orbit_correction = total_processed - len(missing_orbit_data)

print(f"SUMMARY:")
print(f"  Total galaxies to process: {len(halo_ids)}")
print(f"  Successfully processed: {total_processed}")
print(f"  With orbit data (improved method): {successful_orbit_correction}")
print(f"  Without orbit data (fallback method): {len(missing_orbit_data)}")
print(f"  Processing errors: {len(processing_errors)}")

print(f"\nData saved in directory: {output_dir}")

# 如果有缺失轨道数据的情况，再次强调
if missing_orbit_data:
    print("\n" + "!"*80)
    print("!!! ATTENTION !!! Some galaxies are missing orbit data:")
    print(f"Missing orbit data for {len(missing_orbit_data)} galaxies: {missing_orbit_data}")
    print("!"*80)

print("="*80)