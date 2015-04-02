# -*- coding: utf-8 -*-
from PyKEP import planet_ss,epoch,epoch_from_string, planet
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import random
import numpy as np
import pandas as pd
pd.set_option('display.width', 250)

from trajectories import hdf_plot


def get_series(dataframe, parameters, isdate, bounds1, bounds2, mXY, doublefront = None):
    HDF_file = 'protolod.hdf'
    prt = pd.read_hdf(HDF_file, dataframe)

    if isdate[0]:
        maskdf = prt[parameters[0]].map(lambda x : epoch_from_string(x).mjd)
        croped1 = prt[ (maskdf > epoch_from_string(bounds1[0]).mjd)
                     & (maskdf < epoch_from_string(bounds1[1]).mjd) ]
    else:
        croped1 = prt[ (prt[parameters[0]] > bounds1[0])
                     & (prt[parameters[0]] < bounds1[1]) ]

    if isdate[1]:
        maskdf = croped1[parameters[1]].map(lambda x : epoch_from_string(x).mjd)
        croped2 = croped1[ (maskdf > epoch_from_string(bounds2[0]).mjd)
                         & (maskdf < epoch_from_string(bounds2[1]).mjd) ]
    else:
        croped2 = croped1[ (prt[parameters[1]] > bounds2[0])
                         & (prt[parameters[1]] < bounds2[1]) ]

    x = croped2[parameters[0]]
    y = croped2[parameters[1]]

    for z in x,y:
        z.index = range(len(z))

    if isdate[0]:
        x = x.convert_objects(convert_dates='coerce').astype(datetime)
    if isdate[1]:
        y = y.convert_objects(convert_dates='coerce').astype(datetime)

    a = hdf_plot(dataframe)
    xp1, yp1 = a.pareto_frontier(x,y, mXY[0], mXY[1])

    if doublefront is not None:
        if doublefront == 'x':
            xp2, yp2 = a.pareto_frontier(x,y, not mXY[0], mXY[1])
            xp = xp1 + xp2
            yp = yp1 + yp2
            # Sort, in case of broken order in doublefront
            myList = sorted([[xp[i], yp[i]] for i in range(len(xp))], reverse=0)
            xp = [pair[0] for pair in myList]
            yp = [pair[1] for pair in myList]

        elif doublefront == 'y':
            xp2, yp2 = a.pareto_frontier(x,y, mXY[0], not mXY[1])
            xp = xp1 + xp2
            yp = yp1 + yp2
            # Sort, in case of broken order in doublefront
            myList = sorted([[yp[i], xp[i]] for i in range(len(xp))], reverse=0)
            xp = [pair[1] for pair in myList]
            yp = [pair[0] for pair in myList]

    else:
        xp, yp = xp1, yp1

    return x, y, xp, yp


def scalar_vs_date(t, l, ticks, save = False, filename = 'plot', lloc = 'upper left', yminor = False):
    x, y, xp, yp = t
    # plot formating ahead
    plt.figure(figsize=(10, 6), dpi=80)

    plt.plot(x, y,
             'o',
             color='blue',
             label=l[0])

    plt.plot(xp, yp,
             '-',
             color='green',
             label=l[1])

    plt.xlabel(l[2])
    plt.ylabel(l[3])
    plt.legend(loc=lloc, prop={'size':13})

    # Axis and grid modyfications
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator(minticks=ticks[0], maxticks=ticks[1]))
    plt.gca().xaxis.set_minor_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()

    #~ plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    #~ plt.gca().xaxis.set_major_formatter(plt.FormatStrFormatter('%y%m%d'))

    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(ticks[2]))
    plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(ticks[3]))

    plt.gca().xaxis.grid(True,'minor')
    if yminor:
        plt.gca().yaxis.grid(True,'minor')
    plt.gca().xaxis.grid(True,'major', linewidth=1)
    plt.gca().yaxis.grid(True,'major', linewidth=1)

    plt.tight_layout()
    if save:
        plt.savefig(filename+'.png', dpi=200)
    else:
        plt.show()


def scalar_vs_scalar(t, l, ticks, save = False, filename = 'plot', lloc = 'upper left', yminor = False):
    x, y, xp, yp = t
    # plot formating ahead
    plt.figure(figsize=(10, 6), dpi=80)

    plt.plot(x, y,
             'o',
             color='blue',
             label=l[0])

    plt.plot(xp, yp,
             '-',
             color='green',
             label=l[1])

    plt.xlabel(l[2])
    plt.ylabel(l[3])
    plt.legend(loc=lloc, prop={'size':13})

    # Axis and grid modyfications
    #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(ticks[0]))
    plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(ticks[1]))
    plt.gcf().autofmt_xdate()

    #~ plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    #~ plt.gca().xaxis.set_major_formatter(plt.FormatStrFormatter('%y%m%d'))

    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(ticks[2]))
    plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(ticks[3]))

    plt.gca().xaxis.grid(True,'minor')
    if yminor:
        plt.gca().yaxis.grid(True,'minor')
    plt.gca().xaxis.grid(True,'major', linewidth=1)
    plt.gca().yaxis.grid(True,'major', linewidth=1)

    plt.tight_layout()
    if save:
        plt.savefig(filename+'.png', dpi=200)
    else:
        plt.show()


if __name__ == '__main__':
    label_dict = { 'C3' : 'Earth departure C3 [$km^2s^{-2}$]',
                   'DEC' : 'Earth departure asymptote declination [$deg$]',
                   'Ldate' : 'Launch date',
                   'VRE' : 'Earth reentry velocity [$m s^{-1}$]',
                   'MDist' : 'Minimal Mars flyby distance [$km$]',
                   'VDist' : 'Minimal Venus flyby distance [$km$]' }

    # 2018 700-day class Mars Flyby
    ss, sh = 0,0
    if sh:
        parameters = ('C3', 'DEC')
        isdate = (0, 0)
        bounds1 = (20,40)
        bounds2 = (-105, -5)
        t = get_series('eme_700day_2018_single', parameters, isdate, bounds1, bounds2, (0,0), doublefront = 'y')
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [1, 0.25, 1, 0.5]
        scalar_vs_scalar(t, labels, ticks, save = ss, filename = 'eme700.c3.dec', lloc = 'upper right', yminor = 1)

    if sh:
        parameters = ('Ldate', 'C3')
        isdate = (1, 0)
        bounds1 = ('2018-04-01 00:00:00', '2018-06-06 00:00:00')
        bounds2 = (18, 46)
        t = get_series('eme_700day_2018_single', parameters, isdate, bounds1, bounds2, (0,0), doublefront = 'x')
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [3, 12, 0.5, 1]
        scalar_vs_date(t, labels, ticks, save = ss, filename = 'eme700.ldate.c3')

    if sh:
        parameters = ('Ldate', 'DEC')
        isdate = (1, 0)
        bounds1 = ('2018-04-01 00:00:00', '2018-06-06 00:00:00')
        bounds2 = (-60, -5)
        t = get_series('eme_700day_2018_single', parameters, isdate, bounds1, bounds2, (1,1))
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [3, 12, 1,0.5]
        scalar_vs_date(t, labels, ticks, save = ss, filename = 'eme700.ldate.dec', lloc='lower left', yminor = 1)

    if sh:
        parameters = ('Ldate', 'VRE')
        isdate = (1, 0)
        bounds1 = ('2018-04-01 00:00:00', '2018-06-06 00:00:00')
        bounds2 = (12000, 12900)
        t = get_series('eme_700day_2018_single', parameters, isdate, bounds1, bounds2, (1,0))
        labels = [ "Trajectory solutions", "(None)"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [3, 12, 20, 5]
        scalar_vs_date(t, labels, ticks, save = ss, filename = 'eme700.ldate.vre', lloc='upper right')

    if sh:
        parameters = ('C3', 'VRE')
        isdate = (0, 0)
        bounds1 = (0, 46)
        bounds2 = (11600, 14200)
        t = get_series('eme_700day_2018_single', parameters, isdate, bounds1, bounds2, (0,0), doublefront = 'y')
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [1, 0.25, 20, 5]
        scalar_vs_scalar(t, labels, ticks, save = ss, filename = 'eme700.c3.vre', lloc='upper right')

    if sh:
        parameters = ('C3', 'MDist')
        isdate = (0, 0)
        bounds1 = (0, 46)
        bounds2 = (0, 250000)
        t = get_series('eme_700day_2018_single', parameters, isdate, bounds1, bounds2, (0,0))
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [1, 0.25, 10000, 5000]
        scalar_vs_scalar(t, labels, ticks, save = ss, filename = 'eme700.c3.mdist', lloc='upper right')


    # 2017/2018 Mars Flyby
    ss, sh = 0,0
    if sh:
        parameters = ('DEC','C3')
        isdate = (0, 0)
        bounds2 = (38,46)
        bounds1 = (-15, -5)
        t = get_series('eme_500day_2018_single', parameters, isdate, bounds1, bounds2, (0,0))
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [1, 0.25, 1, 0.5]
        scalar_vs_scalar(t, labels, ticks, save = ss, filename = 'eme500.dec.c3', lloc = 'upper right', yminor = 1)

    if sh:
        parameters = ('Ldate', 'C3')
        isdate = (1, 0)
        bounds1 = ('2017-12-01 00:00:00', '2018-01-06 00:00:00')
        bounds2 = (38, 46)
        t = get_series('eme_500day_2018_single', parameters, isdate, bounds1, bounds2, (0,0))
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [3, 12, 1, 0.5]
        scalar_vs_date(t, labels, ticks, save = ss, filename = 'eme500.ldate.c3', lloc = 'upper right', yminor = 1)

    if sh:
        parameters = ('Ldate', 'DEC')
        isdate = (1, 0)
        bounds1 = ('2017-12-01 00:00:00', '2018-01-06 00:00:00')
        bounds2 = (-15, -5)
        t = get_series('eme_500day_2018_single', parameters, isdate, bounds1, bounds2, (0,1))
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [3, 12, 1, 0.5]
        scalar_vs_date(t, labels, ticks, save = ss, filename = 'eme500.ldate.dec', yminor = 1)

    if sh:
        parameters = ('Ldate', 'VRE')
        isdate = (1, 0)
        bounds1 = ('2017-12-01 00:00:00', '2018-01-06 00:00:00')
        bounds2 = (13600, 14200)
        t = get_series('eme_500day_2018_single', parameters, isdate, bounds1, bounds2, (1,0))
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [3, 12, 100, 50]
        scalar_vs_date(t, labels, ticks, save = ss, filename = 'eme500.ldate.vre', yminor = 1)

    if sh:
        parameters = ('VRE', 'C3')
        isdate = (0, 0)
        bounds2 = (38, 46)
        bounds1 = (13600, 14200)
        t = get_series('eme_500day_2018_single', parameters, isdate, bounds1, bounds2, (0,0))
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [100, 25, 1, 0.5]
        scalar_vs_scalar(t, labels, ticks, save = ss, filename = 'eme500.vre.c3', lloc = 'upper right')

    if sh:
        parameters = ('C3', 'MDist')
        isdate = (0, 0)
        bounds1 = (38, 46)
        bounds2 = (0, 250000)
        t = get_series('eme_500day_2018_single', parameters, isdate, bounds1, bounds2, (0,0))
        labels = [ "Trajectory solutions", "(None)"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [1, 0.5, 100, 25]
        scalar_vs_scalar(t, labels, ticks, save = ss, filename = 'eme500.c3.mdist', yminor = 1)


    # 2021 Venus-Mars Flyby
    ss, sh = 0,0
    if sh:
        parameters = ('C3', 'DEC')
        isdate = (0, 0)
        bounds1 = (19,25)
        bounds2 = (-60, 25)
        t = get_series('evme_2021_single', parameters, isdate, bounds1, bounds2, (0,1), doublefront = 'y')
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [1, 0.25, 1,0.5]
        scalar_vs_scalar(t, labels, ticks, save = ss, filename = 'evme.c3.dec', yminor = 1)

    if sh:
        parameters = ('Ldate', 'C3')
        isdate = (1, 0)
        bounds1 = ('2021-11-15 00:00:00', '2021-12-01 00:00:00')
        bounds2 = (19, 25)
        t = get_series('evme_2021_single', parameters, isdate, bounds1, bounds2, (1,0), doublefront = 'x')
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [3, 8, 1, 0.5]
        scalar_vs_date(t, labels, ticks, save = ss, filename = 'evme.ldate.c3', yminor = 1)

    if sh:
        parameters = ('Ldate', 'DEC')
        isdate = (1, 0)
        bounds1 = ('2021-11-10 00:00:00', '2021-12-30 00:00:00')
        bounds2 = (-65, 25)
        t = get_series('evme_2021_single', parameters, isdate, bounds1, bounds2, (0,1), doublefront = 'x')
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [3, 8, 1, 0.5]
        scalar_vs_date(t, labels, ticks, save = ss, filename = 'evme.ldate.dec', yminor = 1)

    if sh:
        parameters = ('Ldate', 'VRE')
        isdate = (1, 0)
        bounds1 = ('2021-11-15 00:00:00', '2022-12-01 00:00:00')
        bounds2 = (19, 250000)
        t = get_series('evme_2021_single', parameters, isdate, bounds1, bounds2, (0,0), doublefront = 'x')
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [3, 8, 100, 50]
        scalar_vs_date(t, labels, ticks, save = ss, filename = 'evme.ldate.vre', yminor = 1)

    if sh:
        parameters = ('C3', 'VRE')
        isdate = (0, 0)
        bounds1 = (19, 25)
        bounds2 = (0, 250000)
        t = get_series('evme_2021_single', parameters, isdate, bounds1, bounds2, (0,0), doublefront = 'x')
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [1, 0.25, 100, 50]
        scalar_vs_scalar(t, labels, ticks, save = ss, filename = 'evme.c3.vre', yminor = 1)

    if sh:
        parameters = ('VDist', 'C3')
        isdate = (0, 0)
        bounds2 = (19, 25)
        bounds1 = (0, 250000)
        t = get_series('evme_2021_single', parameters, isdate, bounds1, bounds2, (0,0), doublefront = 'x')
        labels = [ "Trajectory solutions", "Pareto front approximation"] + [label_dict[parameters[i]] for i in (0,1)]
        ticks = [500, 100, 1, 0.5]
        scalar_vs_scalar(t, labels, ticks, save = ss, filename = 'evme.vdist.c3', lloc = 'lower left', yminor = 1)
