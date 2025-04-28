import os
import numpy as np
from scipy import interpolate
import netCDF4
import Namelist
from xfastran.d3d import shotinfo
from xfastran.d3d import nbidata
from xfastran.d3d import d3dbeam
from xfastran.d3d import nbigeo


def wrt_innubeam(
        shot,
#       tmin,
#       tmax,
        outdir='.',
        Db=0.,
        f_oanb_parm_data='oanb_parm_data.nc',
        f_nubeam_tmpl='innubeam_tmpl_d3d.dat',
        nbi=None):
    # read nbi data from mdsplus
    if nbi is None:
        nbi = nbidata.get_nbi(shot)

    nb_list = nbi['nb_list']
    print(f'\nnb_list = {nb_list}')
    nnbi = len(nb_list)
    time = nbi['time_nbi'][:]
    einj = nbi['einj'][:]
    pinj = nbi['pnbi'][:]
    btilt_15 = nbi['btilt_15']
    stilt_15L = nbi['stilt_15L']
    stilt_15R = nbi['stilt_15R']

    # time average of injection power, beam species mixes
    pinja = np.zeros(nnbi)
    einja = np.zeros(nnbi)
    ffulla = np.zeros(nnbi)
    fhalfa = np.zeros(nnbi)

    # ind = ((time >= tmin) & (time <= tmax))
    # for k in range(nnbi):
    #     pinja[k] = np.average(pinj[k, :][ind])
    #     einja[k] = einj[k]
    #     if pinja[k] < 1.e3:
    #         einja[k] = 20.e3  # <==to prevent nubeam crash
    #     mix = d3dbeam.beam_species_mix(einja[k] * 1.e-3)
    #     ffulla[k] = mix['after'][0]
    #     fhalfa[k] = mix['after'][1]
    #     print(f'{nb_list[k]:6} , pinja = {pinja[k]*1.e-6:6.2}, einja = {einja[k]*1.e-3:6.2f}, ffulla = {ffulla[k]:4.2}, fhalfa = {fhalfa[k]:4.2}')

    for k in range(nnbi):
        pinja[k] = pinj[k][0]
        einja[k] = einj[k]
        if pinja[k] < 1.e3:
            einja[k] = 20.e3  # <==to prevent nubeam crash
        mix = d3dbeam.beam_species_mix(einja[k] * 1.e-3)
        ffulla[k] = mix['after'][0]
        fhalfa[k] = mix['after'][1]
        print(f'{nb_list[k]:6} , pinja = {pinja[k]*1.e-6:6.2}, einja = {einja[k]*1.e-3:6.2f}, ffulla = {ffulla[k]:4.2}, fhalfa = {fhalfa[k]:4.2}')

    # tbona tboffa

    # innubeam namelist
    out_namelist = Namelist.Namelist()

    out_namelist['nubeam_run']['dt_nubeam'] = [0.02]
    out_namelist['nubeam_run']['nstep'] = [20]
    out_namelist['nubeam_run']['navg'] = [10]

    out_namelist['nbi_config']['nbeam'] = [nnbi]
    out_namelist['nbi_config']['abeama'] = nnbi * [2.0]  # hard-coded
    out_namelist['nbi_config']['xzbeama'] = nnbi * [1.0]  # hard-coded

    nubeam_tmpl = Namelist.Namelist(f_nubeam_tmpl)
    for v in nubeam_tmpl['nbi_init'].keys():
        out_namelist['nbi_init'][v] = nubeam_tmpl['nbi_init'][v]
    for v in nubeam_tmpl['nbi_update'].keys():
        out_namelist['nbi_update'][v] = nubeam_tmpl['nbi_update'][v]

    out_namelist['nbi_config']['pinja'] = pinja
    out_namelist['nbi_config']['einja'] = einja
    out_namelist['nbi_config']['ffulla'] = ffulla
    out_namelist['nbi_config']['fhalfa'] = fhalfa

    LTO = shotinfo.get_LTO(shot)
    print(f'shot = {shot}, LTS = {LTO}')
    nubeam_namelist = nbigeo.nubeam_geo(btilt_15, stilt_15L, stilt_15R, LTO, f_oanb_parm_data=f_oanb_parm_data)

    var_list = nubeam_namelist[nb_list[0]].keys()
    for v in var_list:
        out_namelist['nbi_config'][v] = nubeam_namelist[nb_list[0]][v]
    for v in var_list:
        for key in nb_list[1:]:
            out_namelist['nbi_config'][v] += nubeam_namelist[key][v]

    out_namelist['nbi_model']['difb_0'] = [Db]
    out_namelist['nbi_model']['difb_a'] = [Db]
    out_namelist['nbi_model']['difb_in'] = [2]
    out_namelist['nbi_model']['difb_out'] = [2]
    out_namelist['nbi_model']['nkdifb'] = [3]

    f_nubeam_namelist = os.path.join(outdir, f'innubeam_{shot:06}')

    print(f'writing {f_nubeam_namelist}')
    out_namelist.write(f_nubeam_namelist)
