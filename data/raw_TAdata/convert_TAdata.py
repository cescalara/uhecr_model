'''
Converts the ASCII formatted TA public data obtained from Abbasi et al. 2014 (hotspot paper)
into h5py format by retaining the format followed by the Auger dataset currently inside data/UHECRdata.h5
and store the resulting grouping into such file.
This can then be accessed easily in the notebooks 4_fits_to_data.
'''
import os
import numpy as np
import h5py
from datetime import datetime
from astropy import units as u
from astropy.coordinates import SkyCoord

# some paths
curr_dir = os.path.abspath(os.path.dirname(__file__))
data_dir = os.path.dirname(curr_dir)
uhecrdata_path = os.path.join(data_dir, "UHECRdata.h5")

def get_keys(dataset="auger2014"):
    '''
    Get the list of keys required for the TA dataset to be added to the h5py file.
    We check this by opening the UHECR datafile and checking the keys for one dataset.

    This should output:
    ['day', 'dec', 'energy', 'glat', 'glon', 'ra', 'theta', 'year']
    where:
    - day: Julian day (day of year) [1, 365]
    - dec: Declination (deg) [-90, 90]
    - energy: EeV
    - glat: galactic latitude b (deg) [-90, 90]
    - glon: galactic longitude l (deg) [-180, 180]
    - ra: right ascension (deg) [0, 360]
    - theta: zenith angle (deg) [0, 90]
    - year: 2004 - 2009

    These information has been verified with Appendix A of Abreu et al. 2010.

    '''
    # read from the currenht h5py file to check how the format is like
    with h5py.File(uhecrdata_path, "r") as f:
        key_list = list(f[dataset].keys())
    
    return key_list

def create_dataset(data_strlist):
    '''
    Create dataset compartementalized into different groups following the 
    Auger datasets.
    '''
    pass



if __name__ == "__main__":
    # dataset types needed to append to h5py file
    key_list = get_keys()

    # path to TA data
    fname = "TAdata_Egt57.txt"
    label = "TA2015"
    TAdata_path = os.path.join(curr_dir, fname)

    # open TA dataset file
    with open(TAdata_path, "r") as f:
        lines = f.readlines()

    # prelim information is contained in lines 1 - 43
    # so only look at lines 44 onwards
    data_strlist = lines[44:]

    # lists / arrays to append data
    day_list = []
    year_list = []
    zen_list = []
    en_list = []
    ra_list = []
    dec_list = []

    
    # split data into relavant sections
    # date (y m d h m s), zen, en, ra, dec
    # ex. '2008 Jun 10 17 05 37 46.91  88.8  93.50  20.82\n'
    for line in data_strlist:
        line_list = line.split(" ")
        
        date_str = "-".join([l for l in line_list[:6]])
        # date containing year, month, day, hour, minute, second
        date = datetime.strptime(date_str, "%Y-%b-%d-%H-%M-%S")
        day = date.timetuple().tm_yday  # day of year
        year = date.timetuple().tm_year

        # remove any spaces contained within list
        while "" in line_list:
            line_list.remove("")

        zenith_angle, energy, ra, dec = line_list[6:]
        dec = dec[:-1]  # to remove the \n symbol

        # append to list / array
        day_list.append(day)
        year_list.append(year)
        zen_list.append(np.float64(zenith_angle))
        en_list.append(np.float64(energy))
        ra_list.append(np.float64(ra))
        dec_list.append(np.float64(dec))
    

    # get galactic longitudes and latitudes 
    data_skycoord = SkyCoord(ra_list, dec_list, frame="icrs", unit="deg")
    glon_list = data_skycoord.galactic.l.deg
    glat_list = data_skycoord.galactic.b.deg

    # add TA data to existing h5py file
    with h5py.File(os.path.join(data_dir, "UHECRdata.h5"), "r+") as f:
        TAdata_group = f.create_group(label)
        TAdata_group.create_dataset("day", data=np.array(day_list))
        TAdata_group.create_dataset("year", data=np.array(year_list))
        TAdata_group.create_dataset("theta", data=np.array(zen_list))
        TAdata_group.create_dataset("energy", data=np.array(en_list))
        TAdata_group.create_dataset("ra", data=np.array(ra_list))
        TAdata_group.create_dataset("dec", data=np.array(dec_list))
        TAdata_group.create_dataset("glat", data=np.array(glat_list))
        TAdata_group.create_dataset("glon", data=np.array(glon_list))












    





