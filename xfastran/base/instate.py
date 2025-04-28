import os
import numpy as np
import netCDF4
from scipy.interpolate import interp1d, interp2d
from Namelist import Namelist


def slice_0d(nc, key_t, key_f, t0, kind='linear'):
    t = nc.variables[key_t][:]
    f = nc.variables[key_f][:]
    return interp1d(t, f, kind=kind, fill_value=(f[0], f[-1],), bounds_error=False)(t0)


def slice_1d(nc, key_t, key_x, key_f, t0, x0, kind='linear'):
    t = nc.variables[key_t][:]
    x = nc.variables[key_x][:]
    f = nc.variables[key_f][:, :]
    return interp2d(t, x, f.transpose(), kind=kind)(t0, x0).flatten()


def average_0d(nc, key_t, key_f, tmin, tmax):
    t = nc.variables[key_t][:]
    f = nc.variables[key_f][:]
    i = ((t >= tmin) & (t <= tmax))
    return np.average(f[i])


def average_1d(nc, key_t, key_x, key_f, tmin, tmax, x0, kind='linear'):
    t = nc.variables[key_t][:]
    x = nc.variables[key_x][:]
    f = nc.variables[key_f][:, :]
    i = ((t >= tmin) & (t <= tmax))
    f_avg = np.average(f[i], axis=0)
    return interp1d(x, f_avg, kind=kind, fill_value=(f_avg[0], f_avg[-1],), bounds_error=False)(x0)


def nearest_0d(nc, key_t, key_f, t0):
    t = nc.variables[key_t][:]
    f = nc.variables[key_f]
    delta = np.abs(t - t0)
    i = np.argmin(delta)
    print(f'{key_f} at {t[i]}')
    return f[i]


def nearest_1d(nc, key_t, key_x, key_f, t0, kind='linear'):
    t = nc.variables[key_t][:]
    x = nc.variables[key_x][:]
    f = nc.variables[key_f]
    delta = np.abs(t - t0)
    i = np.argmin(delta)
    print(f'{key_f} at {t[i]}')
    return interp1d(x, f[i], kind=kind, fill_value=(f[i][0], f[i][-1],), bounds_error=False)


def wrt_instate(f_timetrace, tmin, tmax, tref=-1, nrho=101, fn_instate='instate', kind='timetrace'):
    nc = netCDF4.Dataset(f_timetrace, 'r', format='NETCDF4')

    if kind == 'snap':
        time_eq = nc.variables['time_eq'][:]
        ip = average_0d(nc, 'time_eq', 'ip', tmin, tmax) 
        b0 = average_0d(nc, 'time_eq', 'bt', tmin, tmax) 
        r0 = average_0d(nc, 'time_eq', 'r0', tmin, tmax) # ok since r0 = const

        rho = np.linspace(0, 1., nrho)
        ne = average_1d(nc, 'time_ne', 'rho_ne', 'ne', tmin, tmax, rho) 
        te = average_1d(nc, 'time_te', 'rho_te', 'te', tmin, tmax, rho) 
        ti = average_1d(nc, 'time_ti', 'rho_ti', 'ti', tmin, tmax, rho) 
        omega = average_1d(nc, 'time_omega', 'rho_omega', 'omega', tmin, tmax, rho) 
        zeff = average_1d(nc, 'time_zeff', 'rho_zeff', 'zeff', tmin, tmax, rho) 

        j_tot = average_1d(nc, 'time_eq', 'rho_eq', 'j_tot', tmin, tmax, rho) 
        p_eq = average_1d(nc, 'time_eq', 'rho_eq', 'p_eq', tmin, tmax, rho) 

        if tref < 0:
            tshape = (tmin + tmax) // 2
        else:
            tshape = tref

        nbdry = nearest_0d(nc, 'time_eq', 'nbdry', tshape)
        rbdry = nearest_0d(nc, 'time_eq', 'rbdry', tshape)[:nbdry]
        zbdry = nearest_0d(nc, 'time_eq', 'zbdry', tshape)[:nbdry]

    elif kind == 'timetrace':
        time_eq = nc.variables['time_eq'][:]
        ip = nearest_0d(nc, 'time_eq', 'ip', tmin) 
        b0 = nearest_0d(nc, 'time_eq', 'bt', tmin) 
        r0 = nearest_0d(nc, 'time_eq', 'r0', tmin)

        rho = np.linspace(0, 1., nrho)
        ne = nearest_1d(nc, 'time_ne', 'rho_ne', 'ne', tmin)(rho)
        te = nearest_1d(nc, 'time_te', 'rho_te', 'te', tmin)(rho) 
        ti = nearest_1d(nc, 'time_ti', 'rho_ti', 'ti', tmin)(rho) 
        omega = nearest_1d(nc, 'time_omega', 'rho_omega', 'omega', tmin)(rho)
        zeff = nearest_1d(nc, 'time_zeff', 'rho_zeff', 'zeff', tmin)(rho) 

        j_tot = nearest_1d(nc, 'time_eq', 'rho_eq', 'j_tot', tmin)(rho) 
        p_eq = nearest_1d(nc, 'time_eq', 'rho_eq', 'p_eq', tmin)(rho) 

        nbdry = nearest_0d(nc, 'time_eq', 'nbdry', tmin)
        rbdry = nearest_0d(nc, 'time_eq', 'rbdry', tmin)[:nbdry]
        zbdry = nearest_0d(nc, 'time_eq', 'zbdry', tmin)[:nbdry]

    rlim = nc.variables['rlim'][:]
    zlim = nc.variables['zlim'][:]

    instate = Namelist()

    instate['instate']['model_shape'] = [0]
    instate['instate']['density_model'] = [0]

    instate['instate']['r0'] = [r0]
    instate['instate']['ip'] = [ip * 1.e-6]
    instate['instate']['b0'] = [abs(b0)]

    instate['instate']['n_ion'] = [1]
    instate['instate']['z_ion'] = [1]
    instate['instate']['a_ion'] = [2]
    instate['instate']['f_ion'] = [1.]

    instate['instate']['n_imp'] = [1]
    instate['instate']['z_imp'] = [6]
    instate['instate']['a_imp'] = [12]
    instate['instate']['f_imp'] = [1.]

    instate['instate']['n_min'] = [0]
    instate['instate']['z_min'] = [1]
    instate['instate']['a_min'] = [1]
    instate['instate']['n_beam'] = [1]
    instate['instate']['z_beam'] = [1]
    instate['instate']['a_beam'] = [2]
    instate['instate']['n_fusion'] = [1]
    instate['instate']['z_fusion'] = [2]
    instate['instate']['a_fusion'] = [4]

    instate['instate']['nrho'] = [nrho]
    instate['instate']['rho'] = rho
    instate['instate']['ne' ] = ne
    instate['instate']['te' ] = te 
    instate['instate']['ti' ] = ti
    instate['instate']['omega' ] = omega
    instate['instate']['zeff' ] = zeff 
    # instate['instate']['p_rad' ] = p_rad 

    instate['instate']['j_tot'] = j_tot
    instate['instate']['p_eq' ] = p_eq
 
    instate['instate']['nbdry' ] = [nbdry]
    instate['instate']['rbdry' ] = rbdry
    instate['instate']['zbdry' ] = zbdry
    instate['instate']['nlim'  ] = [len(rlim)]
    instate['instate']['rlim'  ] = rlim
    instate['instate']['zlim'  ] = zlim

    instate.write(fn_instate)

