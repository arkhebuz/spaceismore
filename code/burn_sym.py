# -*- coding: utf-8 -*-
from math import sqrt, cos, sin
from numpy import matrix, longdouble as ldbl
from scipy.integrate import odeint
from scipy.optimize import fmin as simplex
from sgp4.ext import rv2coe
from consts import *


def diffeq2d(w, t, p):
    """
    Defines the movement equations.
    Arguments:
        w :  vector of the state variables:
                  w = [x, y, vx, vy]
        t :  time
        p :  vector of the parameters:
                  p = [T, m0, A, GMz]
    """
    x,y, vx,vy = w
    T, m0, A, GMz = p

    # Create f = (x', y', vx', vy'):
    FG = GMz/(x*x + y*y)**1.5
    FT = T/( (m0 -A*t)*sqrt(vx*vx + vy*vy) )

    f = [vx,
         vy,
         FT*vx - FG*x,
         FT*vy - FG*y ]
    return f


def do_burn(lus, vecinit, t):
    # Pack up the parameters and initial conditions:
    p = [lus.Thrust*gstd, lus.Stack_mass, lus.Prop_flow, GMz]

    solution = odeint(diffeq2d, vecinit, t, args=(p,), atol=1.0e-8, rtol=1.0e-8)
    vecend = solution[-1]

    xf,yf, vxf,vyf, = vecend
    velocity_sq = vxf**2 + vyf**2
    height = sqrt(xf**2 + yf**2)
    Ep = -GMz/height
    C3 = (velocity_sq + 2*Ep)/1000000
    return C3, vecend


def parkingorbit(pe, long_of_asc_nod, arg_of_peri, fi0, smax, ecc):
    h = pe/(1.0+ecc*cos(fi0))

    r = matrix([ [h*cos(fi0)],
                 [h*sin(fi0)],
                 [0] ])

    v = matrix([ [-sqrt(GMz/pe)*sin(fi0)],
                 [sqrt(GMz/pe)*(ecc+cos(fi0))],
                 [0] ])

    row1 = [ cos(long_of_asc_nod)*cos(arg_of_peri) -sin(long_of_asc_nod)*sin(arg_of_peri),
            -cos(long_of_asc_nod)*sin(arg_of_peri) -sin(long_of_asc_nod)*cos(arg_of_peri),
             0 ]

    row2 = [ sin(long_of_asc_nod)*cos(arg_of_peri) +cos(long_of_asc_nod)*sin(arg_of_peri),
            -sin(long_of_asc_nod)*sin(arg_of_peri) +cos(long_of_asc_nod)*cos(arg_of_peri),
             0 ]

    row3 = [ 0,
             0,
             1 ]

    mactran = matrix([row1, row2, row3])
    rijk = mactran.dot(r).reshape((1,3)).tolist()[0]
    vijk = mactran.dot(v).reshape((1,3)).tolist()[0]

    wekinit = rijk[:2] + vijk[:2]
    return wekinit


def c3(ext_LUS_payload, lus, target_C3 = 0, burns_division = 0.535, rich_return = False):
    lus.set_payload(ext_LUS_payload)
    lus.compute_stack_mass()

    numpoints = 1001
    stoptime = (lus.Usable_prop -lus.Boiled_prop)/lus.Prop_flow
    timetuple = [stoptime*float(i)/(numpoints-1) for i in range(numpoints)]
    brn1_times = timetuple[:int(burns_division*numpoints)]
    brn2_times = timetuple[int(burns_division*numpoints)-1:]

    vecinit = parkingorbit(Rz+lus.Orbit_height, 0, 0, 0, Rz+lus.Orbit_height, 0)
    __, vecend = do_burn(lus, vecinit, brn1_times)
    Rbr = [vecend[0]/1000, vecend[1]/1000, 0]   # m → km
    Vbr = [vecend[2]/1000, vecend[3]/1000, 0]   # m/s → km/s

    p, smax, ecc, __, omega, argp, nu, m, _, _, _ = rv2coe(Rbr, Vbr, GMz_km)
    p = p*1000          # km → m
    smax = smax*1000    # km → m

    #@profile
    def fi_opt(fi0):
        #fi_norm = eph.degrees(fi0[0]).norm
        fi_norm = fi0[0] % (2*pi)
        vecinit_tr = parkingorbit(p, omega, argp, fi_norm, smax, ecc)
        C3, __ = do_burn(lus, vecinit_tr, brn2_times)
        return -1*C3

    a = simplex(fi_opt, -nu, disp=0)
    fi_aft_opt = a[0] % (2*pi)

    vecinit_tr = parkingorbit(p, omega, argp, fi_aft_opt, smax, ecc)
    C3, __ = do_burn(lus, vecinit_tr, brn2_times)

    if rich_return:
        Rpery_tr = smax*(1-ecc)-Rz
        Rapo_tr = smax*(1+ecc)-Rz
        elements = [p, smax, ecc, omega, argp, fi_aft_opt, m]
        period_h = 2*pi*sqrt( (smax**3)/(GMz))/3600
        brn1_lenght = timetuple[int(burns_division*numpoints)-1]
        brn2_lenght = timetuple[-1] - brn1_lenght
        return elements, C3, Rpery_tr, Rapo_tr, period_h, brn1_lenght, brn2_lenght
    else:
        return abs(C3 - target_C3)


if __name__ == '__main__':
    from PyKEP import epoch, epoch_from_string
    from stage import ext_LUS

    lus = ext_LUS('rl')
    lus.set_boiloff(0.002)
    lus.set_launch_date(epoch(3).mjd2000)
    lus.set_time(epoch(8).mjd2000)
    lus.print_info()
    pload = simplex(c3, 27000, args=(lus, 40))[0]
    elem, c3, pery, apo, t, br1, br2 = c3(pload, lus, target_C3=40, rich_return = 1)
    print pload, c3, pery, apo, t, br1, br2
