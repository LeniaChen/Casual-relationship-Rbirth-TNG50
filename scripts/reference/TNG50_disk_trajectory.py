from AnastrisTNG import TNGsimulation
import numpy as np
import h5py
import os
from tqdm import tqdm

#with h5py.File("tng50-1_bar_size.hdf5","r") as disk_file:
#    ID_disk = disk_file["gal_ID"][...]
    
    
def save_dict_to_hdf5(file: h5py.File | h5py.Group | str, data_dict: dict | None = None) ->  h5py.File | h5py.Group | None:
    """Save a dictionary to an HDF5 file or group."""
    if data_dict is None:
        data_dict = {}
    if isinstance(file, str):
        if not os.path.exists(file):
            f = h5py.File(file, "w")
            f.close()
        with h5py.File(file, "r+") as f:
            for i in data_dict:
                f.create_dataset(i, data=data_dict[i])
    if isinstance(file, h5py.File) or isinstance(file, h5py.Group):
        for i in data_dict:
            file.create_dataset(i, data=data_dict[i])
        return file
    return None

def check_path_exists(file_path, path_to_check):
    if not os.path.exists(file_path):
        return False
    
    with h5py.File(file_path, 'r') as f:
        return path_to_check in f

data_file = "TNG50_disk_trajectory.hdf5"  # file name to save galaxy positions

run = 'TNG50' 
BasePath ='/home/tnguser/sims.TNG/'+run + '-1/output'
snap =99
snapshot = TNGsimulation.Snapshot(BasePath,snap)


IDall = [333424,348901,355724,368843,371126,377655,383975,386269,388544,394621,396628,402555,411449,413372,414917,425719,427211,429471,430864,432106,433289,435752,436932,439099,447914,448830,452031,452978,454171,455291,456326,458470,460722,461785,462710,466549,467415,470345,471248,471996,473329,475016,476266,477328,478216,479290,479938,480802,482155,483594,484448,485056,487742,488530,490814,491426,492244,492876,494011,494709,496788,497557,498522,499704,500576,501208,501725,502371,502995,503437,503987,505586,506720,507292,508538,509091,510273,510585,511303,513845,514272,514829,515695,516760,517271,517899,518682,519311,519720,520310,520885,521429,521803,522983,523889,526029,526478,527309,528322,528836,529365,529855,530330,530852,531320,531910,532760,533060,533590,533998,535050,535410,536146,537236,537488,537941,538370,538905,539333,539667,540082,540452,541218,541497,541847,542252,542662,543376,543729,544408,544697,545437,546114,546474,546870,547293,547844,548151,548418,548790,549090,549516,549746,550149,550475,551190,551541,551973,552414,552581,552879,553230,553489,553837,554189,554523,554798,555013,555287,555601,555815,556304,557076,557396,557721,558067,558410,558731,559036,559711,560082,560421,560751,561001,561325,561676,562029,562742,563471,563732,563943,564268,564825,565089,565947,566365,566661,567004,567382,567607,568923,569251,569452,569653,569877,570319,570598,570841,571072,571633,571908,572121,572328,572840,573062,573413,573738,574037,574286,574615,575190,575356,575596,576062,576254,576516,576705,577125,577372,577576,577719,577924,578294,578541,578834,579232,579509,579759,580035,580250,580406,580641,581308,581656,581912,582137,582344,582612,582809,583027,583256,583745,584286,584511,584875,585282,585839,586066,586169,586420,586670,587019,587256,587819,588004,588178,588647,588831,589496,589655,589873,590218,590385,590792,590926,591210,591469,591641,591796,592282,592437,592662,592840,593259,593480,593694,593904,594098,594731,594887,595100,595246,595442,595661,595887,596244,596660,596866,596986,597141,597530,597792,597934,598112,598320,598696,598840,598988,599465,599865,600096,600363,600555,600946,601123,601522,601693,601861,601997,602330,602509,602653,602820,603185,604241,604371,604516,604667,605054,605708,605844,606405,606820,606907,607310,607530,607712,608061,608386,608745,609001,610017,610361,610532,610819,610988,611152,611385,611780,612241,612699,612825,613006,613808,614180,614326,614489,615182,615292,615462,615607,615792,615903,616200,616526,617027,617149,617324,617559,617719,618120,619081,619223,619381,620348,620463,620832,621043,621167,621255,621372,621707,621855,622298,622437,623042,623122,623575,623877,625041,625185,625281,625423,625715,625800,625925,626198,626697,627459,627572,628949,629266,629409,629891,630760,630870,630975,631253,631558,631677,632310,633322,633437,633614,633963,634631,634707,635905,636201,636640,636729,636848,637631,637971,639432,642220,642303,643642,646733,647190,652267,696559,697948]

IDall = np.int32(IDall)# all snap99 galaxy ID to track

for i in tqdm(IDall):
    ID = i
    galaxy_evo = snapshot.galaxy_evolution(ID)
    
    if 'SnapNum' not in galaxy_evo or len(galaxy_evo['SnapNum']) == 0:
        print(f"Warning: No evolution data for galaxy ID {ID}, skipping...")
        continue

    for j in galaxy_evo['SnapNum']:

        snap = j
        
        path_to_check = f"{ID}/{j}"
        if check_path_exists(data_file, path_to_check):
            print(f"Data for galaxy ID {ID}, snap {j} already exists. Skipping...")
            continue
        
        match_indices = galaxy_evo['SnapNum'] == snap
        if not np.any(match_indices):
            print(f"Warning: No match for snapshot {snap} in galaxy {ID}, skipping...")
            continue
        
        subfind_id = int(galaxy_evo['SubfindID'][match_indices][0])
        
        try:
            snapshot_i = TNGsimulation.Snapshot(BasePath,snap)
            
            sub_i=snapshot_i.load_particle(subfind_id)

            if len(sub_i.s)>=100:
                alignwith = "star"
            elif (len(sub_i.g) + len(sub_i.s))>=100:
                alignwith = "baryon"
            else:
                alignwith = "all"

            sub_i.physical_units()
            sub_i.check_boundary()
            mode = 'ssc'
            
            pos_center = sub_i.center(mode=mode)
            sub_i.shift(pos=pos_center)

            re = sub_i.r(0.5,calfor=alignwith)

            vel_center = sub_i.vel_center(pos=[0,0,0.],r_cal=0.5*re)
            sub_i.shift(vel=vel_center)
            
            ang_mom = sub_i.ang_mom_vec(alignwith=alignwith,rmax = 2*re)
            
            
            a = sub_i.properties['a']
            t = sub_i.properties['t'].in_units("Gyr")
            
            galaxy_pos = pos_center.in_units("a kpc",a = sub_i.properties['a'])
            galaxy_vel = vel_center.in_units("a kpc Gyr**-1",a = sub_i.properties['a'])
            
            
            
            
            galaxy_coor = {
                f"{ID}/{j}/a": a,
                f"{ID}/{j}/t": t,
                f"{ID}/{j}/align": alignwith,
                f"{ID}/{j}/re": re,
                f"{ID}/{j}/pos": pos_center,
                f"{ID}/{j}/vel": vel_center,
                f"{ID}/{j}/cpos": galaxy_pos,
                f"{ID}/{j}/cvel": galaxy_vel,
                f"{ID}/{j}/ang_mom": ang_mom,
            }

            save_dict_to_hdf5(data_file,galaxy_coor)
            #print(f"Finished processing galaxy ID {ID}, snap: {snap}")
        except Exception as e:
            print(f"Error processing galaxy ID {ID}, snap {snap}: {str(e)}")
            continue
