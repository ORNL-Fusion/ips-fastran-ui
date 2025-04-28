#!/usr/bin/env python3.8

import os
import sys
import shutil
import numpy as np
import netCDF4
from scipy.interpolate import interp1d, interp2d
from Namelist import Namelist
from xfastran.base import  eqdsk
from xfastran.base.instate import wrt_instate 
from xfastran.d3d import gaprofile
from xfastran.d3d import nbidata
from xfastran.d3d import nubeam
from xfastran.d3d import echdata
import argparse


def input():
    parser = argparse.ArgumentParser()

    parser.add_argument('--snap', 
        action = 'store_true'
        )

    parser.add_argument('--timetrace', 
        action = 'store_true'
        )

    parser.add_argument('--input', 
        action = 'store_true'
        )
    
    parser.add_argument('--template', 
        dest = 'template',
        default = 'timetrace'
        )
    
    parser.add_argument('--shot', 
        dest = 'shot',
        default= 0,
        type=int
        )
    
    parser.add_argument('--time', 
        dest = 'time',
        default= -1,
        type=int
        )

    parser.add_argument('--dtavg', 
        dest = 'dtavg',
        default= -1,
        type=int
        )

    parser.add_argument('--tmin', 
        dest = 'tmin',
        default= -1,
        type=int
        )
    
    parser.add_argument('--tmax', 
        dest = 'tmax',
        default= -1,
        type=int
        )

    parser.add_argument('--tref', 
        dest = 'tref',
        default= -1,
        type=int
        )

    parser.add_argument('--dtbeam', 
        dest = 'dtbeam',
        default= 20,
        type=int
        )

    parser.add_argument('--dtbeam_avg', 
        dest = 'dtbeam_avg',
        default= -1,
        type=int
        )

    parser.add_argument('--dtech', 
        dest = 'dtech',
        default= 20,
        type=int
        )

    parser.add_argument('--dtech_avg', 
        dest = 'dtech_avg',
        default= -1,
        type=int
        )
    
    parser.add_argument('--profdir', 
        dest = 'profdir',
        default = '.'
        )
    
    parser.add_argument('--efitdir', 
        dest = 'efitdir',
        default = '.'
        )

    parser.add_argument('--efit_id', 
        dest = 'efit_id',
        default = ''
        )
    
    parser.add_argument('--rdir', 
        dest = 'rdir',
        default = 'SIMULATION'
        )
    
    parser.add_argument('--nbdry_max', 
        dest = 'nbdry_max',
        default = 101,
        type = int
        )

    parser.add_argument('--ec', 
        action = 'store_true'
        )

    parser.add_argument('--nb', 
        action = 'store_true'
        )

    parser.add_argument('--ifilltime', 
        action = 'store_true'
        )
    
    args = parser.parse_args()

    return args


def main():
    args = input()

    if args.snap and args.timetrace:
       raise Exception('Select one of --snap or --timetrace')
    if not (args.snap or args.timetrace):
       raise Exception('Select one of --snap or --timetrace')

    if args.snap:
       mode = 'snap'
    else:
       mode = 'timetrace'
    print('run mode =', mode)

    shot = args.shot
    efit_id = args.efit_id
    efitdir = args.efitdir
    profdir = args.profdir

    dtbeam = args.dtbeam
    dtbeam_avg = args.dtbeam_avg
    if dtbeam_avg < 0: dtbeam_avg = dtbeam

    dtech = args.dtech
    dtech_avg = args.dtech_avg
    if dtech_avg < 0: dtech_avg = dtech

    if mode == 'snap':
        tmin = args.time - args.dtavg // 2
        tmax = args.time + args.dtavg // 2
        times_nb = np.array([args.time])
        times_ec = np.array([args.time])
    else:
        tmin = args.tmin 
        tmax = args.tmax
        times_nb = np.arange(tmin, tmax + dtbeam, dtbeam)
        times_ec = np.arange(tmin, tmax + dtech, dtech)
    print(f'tmin = {tmin}, tmax = {tmax}')
    print('times_nb =', times_nb)
    print('times_ec =', times_ec)

    tref = args.tref

    nbdry_max = args.nbdry_max

    rdir = os.path.realpath(args.rdir)
    if not os.path.exists(rdir):
        print('create direcotry:', rdir)
        os.makedirs(rdir)

    bin_fluxsurf = os.path.join(os.environ['FLUXSURF_BIN_PATH'], os.environ['FLUXSURF_BIN_NAME'])
    f_oanb_parm_data = os.environ['OANB_PARM_DATA']
    f_nubeam_tmpl = os.environ['INNUBEAM_TEMPLATE']
    xfastran_template_dir = os.environ['XFASTRAN_TEMPLATE_DIR']
    print(xfastran_template_dir)
    
    ncfile = os.path.join(rdir, f't{shot:06}.nc')
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
    nc.close()
    
    print('\n'+72*'=')
    print('= PROCESS EFIT DATA')
    eqdsk.wrt_netcdf(
        shot, 
        efitdir, 
        efit_id, 
        tmin, 
        tmax, 
        nbdry_max, 
        ncfile = ncfile, 
        bin_fluxsurf=bin_fluxsurf
        )

    print('\n'+72*'=')
    print('= PROCESS GAPROFILES')
    prof = gaprofile.readall(
        shot, 
        vars = ['ne', 'te', 'ti', 'omega', 'nz'], 
        tmin = tmin,
        tmax = tmax,
        rdir = profdir,
        ifilltime = args.ifilltime
        )

    gaprofile.wrt_netcdf(
        shot, 
        prof, 
        outdir = '.', 
        ncfile = ncfile
        )
    
    if args.nb:
        print('\n'+72*'=')
        print('= PROCESS NBI DATA')
        nbi = nbidata.get_nbi(
            shot, 
            times_nb,
            dtavg = dtbeam_avg, 
            ncfile = ncfile)

        nubeam.wrt_innubeam(
            shot, 
            outdir = rdir,
            f_oanb_parm_data = f_oanb_parm_data, 
            f_nubeam_tmpl = f_nubeam_tmpl,
            nbi = nbi)
    
    if args.ec:
        print('\n'+72*'=')
        print('= PROCESS EC DATA')
        echdata.get_ech(
            shot, 
            times_ec,
            dtavg = dtech_avg, 
            ncfile = ncfile,
            outdir = rdir
        )


    print('\n'+72*'=')
    print('= WRITE INSTATE')
    fn_instate = os.path.join(rdir, f'instate_{shot:06}')
    wrt_instate(ncfile, tmin, tmax, tref=tref, fn_instate=fn_instate, kind=mode)

    if mode == 'snap' and args.input:
        f_sim = os.path.join(xfastran_template_dir, 'snap', 'fastran_scenario.config')
        print(f_sim)
        shutil.copy(f_sim, rdir)
        input_dir = os.path.join(rdir, 'input')
        if not os.path.exists(input_dir):
            print('create direcotry:', input_dir)
            os.makedirs(input_dir)
        input_files = ['instate']
        if args.nb: input_files += ['innubeam']
        if args.ec: input_files += ['intoray']
        for f in input_files:
            shutil.copyfile( os.path.join(rdir, f'{f}_{shot:06}'), os.path.join(input_dir, f) )
            os.remove( os.path.join(rdir, f'{f}_{shot:06}') )
        input_files = ['infastran', 'infastran0', 'infastran1', 'ingenray_HC', 'inhcd', 'ingenray_LH', 'intglf', 'submitjob.sh', 'submitjob.slurm']
        for f in input_files:
            shutil.copyfile( os.path.join(xfastran_template_dir, 'snap', f), os.path.join(input_dir, f) )
        input_files = ['submitjob.sh', 'submitjob.slurm']
        for f in input_files:
            shutil.copyfile( os.path.join(xfastran_template_dir, 'snap', f), os.path.join(rdir, f) )

if __name__ == '__main__':
   sys.exit(main())
