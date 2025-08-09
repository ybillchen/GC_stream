import os
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import ICRS, SkyCoord
from astropy_healpix import HEALPix

def build_healpix_ref(files):
    healpix_8_min  = [int(file[file.find('_')+1:file.rfind('-')])     for file in files]
    healpix_8_max  = [int(file[file.rfind('-')+1:file.rfind('.csv')]) for file in files]
    reference_file = pd.DataFrame({'file':files, 'healpix8_min':healpix_8_min, 'healpix8_max':healpix_8_max}).reset_index(drop=True)
    
    reference_file['healpix7_min'] = [inp >> 2 for inp in reference_file['healpix8_min']]
    reference_file['healpix7_max'] = [inp >> 2 for inp in reference_file['healpix8_max']]
    
    reference_file['healpix6_min'] = [inp >> 2 for inp in reference_file['healpix7_min']]
    reference_file['healpix6_max'] = [inp >> 2 for inp in reference_file['healpix7_max']]
    
    reference_file['healpix9_min'] = [inp << 2 for inp in reference_file['healpix8_min']]
    reference_file['healpix9_max'] = [(inp << 2) + 3 for inp in reference_file['healpix8_max']]
    
    ncols = ['file', 'healpix6_min', 'healpix6_max', 'healpix7_min', 'healpix7_max', 'healpix8_min', 'healpix8_max', 'healpix9_min', 'healpix9_max']
    return reference_file[ncols]

def select_files(ra, dec, radius, files=None, reference_file=None, hpx_level=8, filepath=''):
    center = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame=ICRS())

    if files is None:
        files = sorted([f.name for f in Path(filepath).glob("*.csv")])
    
    if reference_file is None:
        reference_file = build_healpix_ref(files)

    hp = HEALPix(nside=2**hpx_level, order='nested', frame=ICRS())
    hp_cone_search = hp.cone_search_skycoord(center, radius=radius*u.deg)

    subset = []
    for index in reference_file.index:
        row = reference_file.iloc[index]
        hp_min, hp_max = row[f'healpix{hpx_level}_min'], row[f'healpix{hpx_level}_max']
        if np.any(np.logical_and(hp_cone_search >= hp_min, hp_cone_search <= hp_max)):
           subset.append(row['file'])
    return subset

def load_subset(subset, filepath=''):
    ds = []
    for file in subset:
        newfile = file.replace('GaiaSource_', 'GaiaSourceReduced_').replace('csv', 'npy')
        newfilename = os.path.join(filepath, newfile)
        ds.append(np.load(newfilename))
    return np.hstack(ds)
    
def load_phot_subset(subset, filepath=''):
    ds = []
    for file in subset:
        newfile = file.replace('GaiaSource_', 'GaiaSourcePhot_').replace('csv', 'npy')
        newfilename = os.path.join(filepath, newfile)
        ds.append(np.load(newfilename))
    return np.hstack(ds)