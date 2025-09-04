import numpy as np
import pandas as pd

def load_GC():
	# load MW metal and mass

	mass_mw = []
	name_mass_mw = []
	feh_mw = []
	name_feh_mw = []
	ebv_mw = []
	name_ebv_mw = []
	f_mw = '../data/combined_table.txt'
	with open(f_mw, 'rb') as reader:
	    lines = reader.readlines()[3:170]
	    for line in lines:
	        mass_mw.append(np.log10(float(line[80:91])))
	        name_mass_mw.append(str(line[0:13].strip())[2:-1].replace(' ', '_'))
	f_mw = '../data/mwgc.dat'
	with open(f_mw, 'rb') as reader:
	    lines = reader.readlines()[252:409]
	    for line in lines:
	        fehh = str(line[13:19].strip())[2:-1]
	        if len(fehh) == 0:
	            continue
	        feh_mw.append(float(fehh))
	        name_feh_mw.append(str(line[1:13].strip())[2:-1].replace(' ', '_'))
	    for line in lines:
	        ebv = str(line[24:29].strip())[2:-1]
	        if len(ebv) == 0:
	            continue
	        ebv_mw.append(float(ebv))
	        name_ebv_mw.append(str(line[1:13].strip())[2:-1].replace(' ', '_'))
	mass_mw = np.array(mass_mw)
	name_mass_mw = np.array(name_mass_mw)
	df = pd.DataFrame({
		'name': name_mass_mw,
		'logm': mass_mw,
	}).set_index('name')
	feh_mw = np.array(feh_mw)
	name_feh_mw = np.array(name_feh_mw)
	df = df.join(pd.DataFrame({
		'name': name_feh_mw,
		'feh': feh_mw,
	}).set_index('name'))
	ebv_mw = np.array(ebv_mw)
	name_ebv_mw = np.array(name_ebv_mw)
	df = df.join(pd.DataFrame({
		'name': name_ebv_mw,
		'ebv': ebv_mw,
	}).set_index('name'))

	# load mw pos and vel

	dtype5 = {
	    'names': (
	        'name', 'RA', 'DEC', 'l', 'b', 'Rsun', 'ERsun', 'R_GC', '<RV>', 'ERV',
	        'mu_alpha', 'Dmu_alpha', 'mu_delta', 'Dmu_delta', 'rhopmrade', 'X', 'DX', 
	        'Y', 'DY', 'Z', 'DZ', 'U', 'DU', 'V', 'DV', 'W', 'DW', 'RPERI', 'DRPERI', 'RAPO', 'DRAPO'),
	    'formats': (
	        'U16', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 
	        'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
	        'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}
	path = '../data/orbits_table_v4.txt'
	data_mw = np.loadtxt(path, dtype=dtype5)
	df = df.join(pd.DataFrame(data_mw).set_index('name'))

	# load structure

	dtype0 = {
	    'names': (
	        'name', 'RA', 'DEC', 'R_Sun', 'DR_Sun', 'R_GC', 'DR_GC', 'N_RV', 'N_PM', 'Mass', 'DM',
	        'V', 'DV', 'M/L_V', 'DM/L', 'rc', 'rh,l', 'rh,m', 'rt', 'rho_c', 'rho_h,m', 'sig_c', 
	        'sig_h,m', 'lg Trh', 'lg Mini',  'T_Diss', 'M_Low', 'M_High', 'MF', 'symbol', 'DMF', 
	        'sig0', 'vesc', 'etac', 'etah'),
	    'formats': (
	        'U16', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
	        'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',
	        'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U4', 'f8', 
	        'f8', 'f8', 'f8', 'f8')}
	path = '../data/combined_table.txt'
	data_structure_mw = np.loadtxt(path, dtype=dtype0)
	df = df.join(pd.DataFrame(data_structure_mw).set_index('name'), lsuffix='', rsuffix='_structure')

	return df