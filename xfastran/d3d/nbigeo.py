import os
import numpy as np
from scipy import interpolate
import netCDF4
import Namelist


def find_tangent(p0, p1, p2):
    x0 = p0[0]
    y0 = p0[1]
    x1 = p1[0]
    y1 = p1[1]
    x2 = p2[0]
    y2 = p2[1]
    s = (x1 - x0) * (x2 - x1) + (y1 - y0) * (y2 - y1)
    s /= (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)
    s *= -1
    x = x1 + (x2 - x1) * s
    y = y1 + (y2 - y1) * s
    p = [x, y]
    return np.array(p)


def rzt2xyz(p):
    deg2rad = (2.0 * np.pi) / 360.
    r = p[0]
    z = p[1]
    t = deg2rad * p[2]
    pp = [r * np.cos(t), r * np.sin(t), z]
    return np.array(pp)


def length(p1, p2):
    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**0.5


def nbigeo(ps_rzt, pe_rzt, src_angle, src, module):
    # constant
    deg2rad = (2. * np.pi) / 360.
    rad2deg = 360. / (2. * np.pi)

    # 150NB source geometry
    s_w = 5.
    s_h = 6.
    alpha0 = 0.55 * deg2rad
    alpha1 = 1.50 * deg2rad

    # beamlet divergence
    div_w = 0.5
    div_h = {'L': 1.0, 'R': 1.3}[src]

    # focal length
    focl_w = 1.e33
    focl_h = 1.e33

    # 150NB absolute collimator
    l_col = 483.0
    col_rect = [(14.376 - 2.125) * np.cos(deg2rad * 4.33),
                (14.376 + 2.125) * np.cos(deg2rad * 4.33),
                20.599,
                17.381]

    # 150NB port bottom edge
    rp = 244.6
    zp = -27.45
    p_w = 50.0
    p_h = 50.0

    # output structure
    n_module = 4
    varlist = [
        'rtcena', 
        'xlbtna', 
        'xybsca',
        'nbshapa', 
        'bmwidra', 
        'bmwidza',
        'nbapsha',
        'xlbapa', 
        'xybapa',
        'rapedga', 
        'xzpedga',
        'xrapoffa', 
        'xzapoffa',
        'divra', 
        'divza', 
        'foclra', 
        'foclza',
        'nbapsh2',
        'rapedg2', 
        'xzpedg2',
        'xlbapa2', 
        'xrapoff2', 
        'xzapoff2',
        'nlco']
    nubgeo = {}
    for var in varlist:
        nubgeo[var] = []

    # source
    # alpha = [alpha1, alpha0, -alpha0, -alpha1]
    alpha = {
        'up': alpha1,
        'mu': alpha0,
        'ml': -alpha0,
        'lo': -alpha1}[module]

    xoffset = {
        'up': 2.0 * s_h * np.sin(alpha0) + s_h * np.sin(alpha1),
        'mu': s_h * np.sin(alpha0),
        'ml': s_h * np.sin(alpha0),
        'lo': 2.0 * s_h * np.sin(alpha0) + s_h * np.sin(alpha1)
    }[module]

    yoffset = {
        'up': 2.0 * s_h * np.cos(alpha0) + s_h * np.cos(alpha1),
        'mu': s_h * np.cos(alpha0),
        'ml': -s_h * np.cos(alpha0),
        'lo': -2.0 * s_h * np.cos(alpha0) - s_h * np.cos(alpha1)
    }[module]

    # collimator
    da0 = col_rect[0]
    db0 = col_rect[1]
    dc0 = col_rect[2]
    dd0 = col_rect[3]

    dh = l_col * np.tan(deg2rad * src_angle)

    # beamline
    ps = rzt2xyz(ps_rzt)
    pe = rzt2xyz(pe_rzt)
    angle = np.arctan((ps[2] - pe[2]) / length(ps[0:2], pe[0:2]))
    pt = find_tangent([0.0, 0.0], ps[0:2], pe[0:2])

    r_t = length([0.0, 0.0], pt)
    l_st = length(ps[0:2], pt) / np.cos(angle)
    z_s = ps[2]

    # aperture 1: collimator
    x0 = xoffset
    y0 = yoffset
    x1 = l_col / np.cos(deg2rad * src_angle)
    y1 = -np.tan(alpha) * (x1 - x0) + y0
    l_sa = ((x1 - x0)**2 + (y1 - y0)**2)**0.5
    dc = (dc0 - dh + y1) * np.cos(alpha)
    dd = (dd0 + dh - y1) * np.cos(alpha)

    a_h = 0.5 * (dc + dd)
    da_h = a_h - dc
    a_w = 0.5 * (da0 + db0)
    da_w = a_w - da0
    z_a = z_s - l_sa * np.sin(angle)

    # aperture 2: port
    x1 = ps[0]
    y1 = ps[1]
    z1 = ps[2]
    x2 = pe[0]
    y2 = pe[1]
    z2 = pe[2]

    a = (x2 - x1)**2 + (y2 - y1)**2
    b = 2 * ((x2 - x1) * x1 + (y2 - y1) * y1)
    c = x1**2 + y1**2 - rp**2

    l1 = (-b + (b**2 - 4.0 * a * c)**0.5) / (2.0 * a)
    l2 = (-b - (b**2 - 4.0 * a * c)**0.5) / (2.0 * a)
    l = min([l1, l2])

    x = x1 + l * (x2 - x1)
    y = y1 + l * (y2 - y1)
    z = z1 + l * (z2 - z1)
    r = (x**2 + y**2)**0.5
    L = ((x - x1)**2 + (y - y1)**2 + (z - z1)**2)**0.5

    l_sa2 = L
    a_w2 = p_w
    da_w2 = 0.0
    a_h2 = p_h
    da_h2 = a_h2 - (z - zp)

    # output
    nubgeo = {}
    nubgeo['rtcena'] = [r_t]
    nubgeo['xlbtna'] = [l_st]
    nubgeo['xybsca'] = [z_s]
    nubgeo['nbshapa'] = [1]
    nubgeo['bmwidra'] = [s_w]
    nubgeo['bmwidza'] = [s_h]
    nubgeo['divra'] = [deg2rad * div_w]
    nubgeo['divza'] = [deg2rad * div_h]
    nubgeo['foclra'] = [focl_w]
    nubgeo['foclza'] = [focl_h]
    nubgeo['nbapsha'] = [1]
    nubgeo['xlbapa'] = [l_sa]
    nubgeo['xybapa'] = [z_a]
    nubgeo['rapedga'] = [a_w]
    nubgeo['xzpedga'] = [a_h]
    nubgeo['xrapoffa'] = [da_w]
    nubgeo['xzapoffa'] = [da_h]
    nubgeo['nbapsh2'] = [1]
    nubgeo['rapedg2'] = [a_w2]
    nubgeo['xzpedg2'] = [a_h2]
    nubgeo['xlbapa2'] = [l_sa2]
    nubgeo['xrapoff2'] = [da_w2]
    nubgeo['xzapoff2'] = [da_h2]
    nubgeo['nlco'] = [True]

    return nubgeo


def interp2d(x, y, z, x0, y0):
    return interpolate.RectBivariateSpline(x, y, z, kx=1, ky=1)(x0, y0)[0][0]


def nubeam_geo(
        btilt,
        stilt_L,
        stilt_R,
        LTO,
        f_oanb_parm_data='oanb_parm_data.nc'):
    nubeam_namelist = {}

    # -- 15L, 15R
    if LTO >= 2:
        oanb = netCDF4.Dataset(f_oanb_parm_data, 'r', format='NETCDF4')
        x = oanb.variables['stilt'][:]
        y = oanb.variables['btilt'][:]
        for src in ['l', 'r']:
            for segment in ['up', 'mu', 'ml', 'lo']:
                id = 'l' + segment
                rs = oanb.variables['rs_{}_{}'.format(src, segment)][:, :]
                zs = oanb.variables['zs_{}_{}'.format(src, segment)][:, :]
                ts = oanb.variables['ts_{}_{}'.format(src, segment)][:, :]

                ra = oanb.variables['ra_{}_{}'.format(src, segment)][:, :]
                za = oanb.variables['za_{}_{}'.format(src, segment)][:, :]
                ta = oanb.variables['ta_{}_{}'.format(src, segment)][:, :]

                ps = [interp2d(x, y, rs, stilt_L, btilt),
                      interp2d(x, y, zs, stilt_L, btilt),
                      interp2d(x, y, ts, stilt_L, btilt)]

                pa = [interp2d(x, y, ra, stilt_L, btilt),
                      interp2d(x, y, za, stilt_L, btilt),
                      interp2d(x, y, ta, stilt_L, btilt)]

                nubeam_namelist['15{}_{}'.format(src.upper(), segment.upper())] = nbigeo(
                    ps, pa, stilt_L, src.upper(), segment)
        oanb.close()
    else:
        # 15L
        nubeam_namelist['15L'] = {
            'nlco': [True],
            'nbshapa': [1],
            'bmwidra': [6.0],
            'bmwidza': [24.0],
            'foclra': [1.0e33],
            'foclza': [1.0e3],
            'divra': [0.00873],
            'divza': [0.0227],
            'rtcena': [114.6],
            'xlbtna': [802.8],
            'xybsca': [0.0],
            # 'xbzeta'  : [314.289],
            'xlbapa': [186.1],
            'xybapa': [0.0],
            'nbapsha': [1],
            'rapedga': [8.85],
            'xzpedga': [24.0],
            'xrapoffa': [0.0],
            'xzapoffa': [0.0],
            'nbapsh2': [0],
            'rapedg2': [0.0],
            'xzpedg2': [0.0],
            'xlbapa2': [0.0],
            'xrapoff2': [0.0],
            'xzapoff2': [0.0]
        }

        # 15R
        nubeam_namelist['15R'] = {
            'nlco': [True],
            'nbshapa': [1],
            'bmwidra': [6.0],
            'bmwidza': [24.0],
            'foclra': [1.0e33],
            'foclza': [1.0e3],
            'divra': [0.00873],
            'divza': [0.0227],
            'rtcena': [76.2],
            'xlbtna': [817.3],
            'xybsca': [0.0],
            # 'xbzeta'  : [320.16 ],
            'xlbapa': [186.1],
            'xybapa': [0.0],
            'nbapsha': [1],
            'rapedga': [8.85],
            'xzpedga': [24.0],
            'xrapoffa': [0.0],
            'xzapoffa': [0.0],
            'nbapsh2': [0],
            'rapedg2': [0.0],
            'xzpedg2': [0.0],
            'xlbapa2': [0.0],
            'xrapoff2': [0.0],
            'xzapoff2': [0.0]
        }

    # -- 30L
    nubeam_namelist['30L'] = {
        'nlco': [True],
        'nbshapa': [1],
        'bmwidra': [6.0],
        'bmwidza': [24.0],
        'foclra': [1.0e33],
        'foclza': [1.0e3],
        'divra': [0.00873],
        'divza': [0.0227],
        'rtcena': [114.6],
        'xlbtna': [802.8],
        'xybsca': [0.0],
        # 'xbzeta'  : [314.289],
        'xlbapa': [186.1],
        'xybapa': [0.0],
        'nbapsha': [1],
        'rapedga': [8.85],
        'xzpedga': [24.0],
        'xrapoffa': [0.0],
        'xzapoffa': [0.0],
        'nbapsh2': [0],
        'rapedg2': [0.0],
        'xzpedg2': [0.0],
        'xlbapa2': [0.0],
        'xrapoff2': [0.0],
        'xzapoff2': [0.0]
    }

    # -- 30R
    nubeam_namelist['30R'] = {
        'nlco': [True],
        'nbshapa': [1],
        'bmwidra': [6.0],
        'bmwidza': [24.0],
        'foclra': [1.0e33],
        'foclza': [1.0e3],
        'divra': [0.00873],
        'divza': [0.0227],
        'rtcena': [76.2],
        'xlbtna': [817.3],
        'xybsca': [0.0],
        # 'xbzeta'  : [320.16 ],
        'xlbapa': [186.1],
        'xybapa': [0.0],
        'nbapsha': [1],
        'rapedga': [8.85],
        'xzpedga': [24.0],
        'xrapoffa': [0.0],
        'xzapoffa': [0.0],
        'nbapsh2': [0],
        'rapedg2': [0.0],
        'xzpedg2': [0.0],
        'xlbapa2': [0.0],
        'xrapoff2': [0.0],
        'xzapoff2': [0.0]
    }

    # -- 21L, 21R
    if LTO >= 3:
        # 21L
        _nubeam_namelist = {
            'nlco': [True, True, True, True],
            'rtcena': [122.456917576, 122.39843872, 122.291059658, 122.243912327],
            'xlbtna': [850.214028782, 849.524257065, 848.148871449, 847.559317417],
            'xybsca': [208.02, 196.72, 185.34, 173.9],
            'nbshapa': [1, 1, 1, 1],
            'bmwidra': [6.0, 6.0, 6.0, 6.0],
            'bmwidza': [6.0, 6.0, 6.0, 6.0],
            'nbapsha': [1, 1, 1, 1],
            'rapedga': [12.1943563937, 12.1943563937, 12.1943563937, 12.1943563937],
            'xzpedga': [20.8428552257, 20.8490393773, 20.8490393773, 20.8428552257],
            'xlbapa': [453.683214255, 453.76331127, 453.76331127, 453.683214255],
            'xybapa': [52.8362001721, 48.5926639197, 45.4781981306, 41.2352837522],
            'xrapoffa': [4.04777542824, 4.04777542824, 4.04777542824, 4.04777542824],
            'xzapoffa': [-4.46981680614, 0.00602435873612, 3.29382360026, 7.76868597856],
            'divra': [0.00872664625997, 0.00872664625997, 0.00872664625997, 0.00872664625997],
            'divza': [0.0174532925199, 0.0174532925199, 0.0174532925199, 0.0174532925199],
            'foclra': [1e+33, 1e+33, 1e+33, 1e+33],
            'foclza': [1e+33, 1e+33, 1e+33, 1e+33],
            'nbapsh2': [1, 1, 1, 1],
            'rapedg2': [50.0, 50.0, 50.0, 50.0],
            'xzpedg2': [50.0, 50.0, 50.0, 50.0],
            'xlbapa2': [624.883056468, 625.477413129, 625.472461949, 626.013562484],
            'xrapoff2': [0.0, 0.0, 0.0, 0.0],
            'xzapoff2': [28.2732554343, 30.0120056491, 29.9970838718, 31.7070517593],
        }
        for k, segment in enumerate(['up', 'mu', 'ml', 'lo']):
            src = '21L_{}'.format(segment.upper())
            nubeam_namelist[src] = {}
            for key in _nubeam_namelist.keys():
                nubeam_namelist[src][key] = [_nubeam_namelist[key][k]]

        # 21R
        _nubeam_namelist = {
            'nlco': [True, True, True, True],
            'rtcena': [77.6306180199, 77.6900136326, 77.8094778528, 78.8669484228],
            'xlbtna': [867.374512465, 866.495451488, 864.90199983, 864.032624085],
            'xybsca': [208.02, 196.72, 185.34, 173.9],
            'nbshapa': [1, 1, 1, 1],
            'bmwidra': [6.0, 6.0, 6.0, 6.0],
            'bmwidza': [6.0, 6.0, 6.0, 6.0],
            'nbapsha': [1, 1, 1, 1],
            'rapedga': [12.2071139209, 12.2071139209, 12.2071139209, 12.2071139209],
            'xzpedga': [20.8428552257, 20.8490393773, 20.8490393773, 20.8428552257],
            'xlbapa': [453.683214255, 453.76331127, 453.76331127, 453.683214255],
            'xybapa': [52.8453805851, 48.6034330998, 45.4887708784, 41.2492484256],
            'xrapoffa': [-4.08604800995, -4.08604800995, -4.08604800995, -4.08604800995],
            'xzapoffa': [-4.46981680614, 0.00602435873612, 3.29382360026, 7.76868597856],
            'divra': [0.00872664625997, 0.00872664625997, 0.00872664625997, 0.00872664625997],
            'divza': [0.0209439510239, 0.0209439510239, 0.0209439510239, 0.0209439510239],
            'foclra': [1e+33, 1e+33, 1e+33, 1e+33],
            'foclza': [1e+33, 1e+33, 1e+33, 1e+33],
            'nbapsh2': [1, 1, 1, 1],
            'rapedg2': [50.0, 50.0, 50.0, 50.0],
            'xzpedg2': [50.0, 50.0, 50.0, 50.0],
            'xlbapa2': [620.533019759, 621.121100267, 621.141864848, 621.915655011],
            'xrapoff2': [0.0, 0.0, 0.0, 0.0],
            'xzapoff2': [26.7727547458, 28.5751817186, 28.6478071129, 30.489610682],
        }
        for k, segment in enumerate(['up', 'mu', 'ml', 'lo']):
            src = '21R_{}'.format(segment.upper())
            nubeam_namelist[src] = {}
            for key in _nubeam_namelist.keys():
                nubeam_namelist[src][key] = [_nubeam_namelist[key][k]]
    elif LTO >= 1:
        # 21L
        nubeam_namelist['21L'] = {
            'nlco': [False],
            'nbshapa': [1],
            'bmwidra': [6.0],
            'bmwidza': [24.0],
            'foclra': [1.0e33],
            'foclza': [1.0e3],
            'divra': [0.00873],
            'divza': [0.0227],
            'rtcena': [76.2],
            'xlbtna': [817.3],
            'xybsca': [0.0],
            # 'xbzeta'  : [159.84 ],
            'xlbapa': [186.1],
            'xybapa': [0.0],
            'nbapsha': [1],
            'rapedga': [8.85],
            'xzpedga': [24.0],
            'xrapoffa': [0.0],
            'xzapoffa': [0.0],
            'nbapsh2': [0],
            'rapedg2': [0.0],
            'xzpedg2': [0.0],
            'xlbapa2': [0.0],
            'xrapoff2': [0.0],
            'xzapoff2': [0.0]
        }

        # 21R
        nubeam_namelist['21R'] = {
            'nlco': [False],
            'nbshapa': [1],
            'bmwidra': [6.0],
            'bmwidza': [24.0],
            'foclra': [1.0e33],
            'foclza': [1.0e3],
            'divra': [0.00873],
            'divza': [0.0227],
            'rtcena': [114.6],
            'xlbtna': [802.8],
            'xybsca': [0.0],
            # 'xbzeta'  : [165.711],
            'xlbapa': [186.1],
            'xybapa': [0.0],
            'nbapsha': [1],
            'rapedga': [8.85],
            'xzpedga': [24.0],
            'xrapoffa': [0.0],
            'xzapoffa': [0.0],
            'nbapsh2': [0],
            'rapedg2': [0.0],
            'xzpedg2': [0.0],
            'xlbapa2': [0.0],
            'xrapoff2': [0.0],
            'xzapoff2': [0.0]
        }
    else:
        # 21L
        nubeam_namelist['21L'] = {
            'nlco': [True],
            'nbshapa': [1], 
            'bmwidra': [6.0], 
            'bmwidza': [24.0],
            'foclra': [1.0e33], 
            'foclza': [1.0e3],
            'divra': [0.00873], 
            'divza': [0.0227],
            'rtcena': [114.6], 
            'xlbtna': [802.8], 
            'xybsca': [0.0],
            # 'xbzeta'  : [165.711],
            'xlbapa': [186.1], 
            'xybapa': [0.0],
            'nbapsha': [1], 
            'rapedga': [8.85], 
            'xzpedga': [24.0],
            'xrapoffa': [0.0], 
            'xzapoffa': [0.0],
            'nbapsh2': [0], 
            'rapedg2': [0.0], 
            'xzpedg2': [0.0],
            'xlbapa2': [0.0], 
            'xrapoff2': [0.0], 
            'xzapoff2': [0.0]
        }

        # 21R
        nubeam_namelist['21R'] = {
            'nlco': [True],
            'nbshapa': [1], 
            'bmwidra': [6.0], 
            'bmwidza': [24.0],
            'foclra': [1.0e33], 
            'foclza': [1.0e3],
            'divra': [0.00873], 
            'divza': [0.0227],
            'rtcena': [76.2], 
            'xlbtna': [817.3], 
            'xybsca': [0.0],
          # 'xbzeta'  : [159.84 ],
            'xlbapa': [186.1], 
            'xybapa': [0.0],
            'nbapsha': [1], 
            'rapedga': [8.85], 
            'xzpedga': [24.0],
            'xrapoffa': [0.0], 
            'xzapoffa': [0.0],
            'nbapsh2': [0], 
            'rapedg2': [0.0], 
            'xzpedg2': [0.0],
            'xlbapa2': [0.0], 
            'xrapoff2': [0.0], 
            'xzapoff2': [0.0]
        }

    # -- 33L
    nubeam_namelist['33L'] = {
        'nlco': [True],
        'nbshapa': [1], 
        'bmwidra': [6.0], 
        'bmwidza': [24.0],
        'foclra': [1.0e33], 
        'foclza': [1.0e3],
        'divra': [0.00873], 
        'divza': [0.0227],
        'rtcena': [114.6], 
        'xlbtna': [802.8], 
        'xybsca': [0.0],
      # 'xbzeta'  : [14.2889],
        'xlbapa': [186.1], 
        'xybapa': [0.0],
        'nbapsha': [1], 
        'rapedga': [8.85], 
        'xzpedga': [24.0],
        'xrapoffa': [0.0], 
        'xzapoffa': [0.0],
        'nbapsh2': [0], 
        'rapedg2': [0.0], 
        'xzpedg2': [0.0],
        'xlbapa2': [0.0], 
        'xrapoff2': [0.0], 
        'xzapoff2': [0.0]
    }

    # -- 33R
    nubeam_namelist['33R'] = {
        'nlco': [True],
        'nbshapa': [1], 
        'bmwidra': [6.0], 
        'bmwidza': [24.0],
        'foclra': [1.0e33], 
        'foclza': [1.0e3],
        'divra': [0.00873], 
        'divza': [0.0227],
        'rtcena': [76.2], 
        'xlbtna': [817.3], 
        'xybsca': [0.0],
      # 'xbzeta'  : [20.1602],
        'xlbapa': [186.1], 
        'xybapa': [0.0],
        'nbapsha': [1], 
        'rapedga': [8.85], 
        'xzpedga': [24.0],
        'xrapoffa': [0.0], 
        'xzapoffa': [0.0],
        'nbapsh2': [0], 
        'rapedg2': [0.0], 
        'xzpedg2': [0.0],
        'xlbapa2': [0.0], 
        'xrapoff2': [0.0], 
        'xzapoff2': [0.0]
    }

    # return
    return nubeam_namelist
