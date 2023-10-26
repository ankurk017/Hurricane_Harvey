import numpy as np
import xarray as xr

def_folder = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/'

case = 'post'
wrf_input = def_folder+f'{case}/WRF/test/em_real/wrfinput_d02'


postfix = [case+'_'+values for values in ('min1', 'min2', 'min3', 'plu1', 'plu2', 'plu3')]
factor = (0.99, 0.98, 0.97, 1.01, 1.02, 1.03)


for postfix_values, factor_values in zip(postfix, factor):
	outout_file = f'../../perturbation/wrfinput_d02_{postfix_values}'
	print(outout_file)
	A = xr.open_dataset(wrf_input)
	A['QVAPOR'].values[:, A['QVAPOR'].bottom_top == 1, :, :] *= factor_values
	A.to_netcdf(outout_file)



