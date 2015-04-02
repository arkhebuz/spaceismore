# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import rc
from scipy.optimize import fmin as simplex
from datetime import datetime, timedelta

from burn_sym import *
from stage import ext_LUS


# IM C3, list start @ 15.12.2017, list end @ 09.01.2018, 04.01.2018 (day 20) -> renetry speed cutoff
allc3req = [43.11, 42.77, 42.45, 42.12, 41.81,
            41.51, 41.22, 40.94, 40.68, 40.42,
            40.19, 39.96, 39.76, 39.58, 39.40,
            39.26, 39.14, 39.02, 38.93, 38.86,
            38.81, 38.79, 38.80, 38.85, 38.98, 39.24]
alltc3req = range(26)


def swp_time(x):
    return datetime(2017, 12, 15) + timedelta(x)


def all_the_rest(lus, offset, loiter, save = 0):
    lus.set_launch_date(offset)
    lus.set_time(offset+loiter)
    lus.print_info()

    print "Launching numerical simualtion..."
    C3_target = allc3req[offset+loiter]
    ext_LUS_payload = simplex(c3, 27000, args=(lus, C3_target))[0]
    print "         Payload: {0} kg @ [{1}off, {2}loi] ".format(round(ext_LUS_payload,2), offset, loiter)

    # Get parking orbit parameters
    __, __, pery, apo, transorb_period, br_time_1, br_time_2 = c3(ext_LUS_payload, lus, target_C3 = C3_target, rich_return = 1)

    # Plot generation
    print "\nDeploying C3@{0}% boiloff rate chart...".format(100*lus.Boiloff)

    # C3 plot data generation.
    allc3 = []
    newt = []
    for day in range(loiter+1):
        lus.set_time(offset+day)
        allc3.append(c3(ext_LUS_payload, lus))
        newt.append(day+offset)

    # Chart name.
    name = "C3 @ [{0}%/d boiloff rate, {1}-day loiter time]".format(lus.Boiloff*100, loiter)

    # ext_LUS C3 curve annotation text at start.
    anno2_text = r"$Diff:\ {0} \ km^{{2}}s^{{-2}}$".format(round(allc3[0]-allc3req[offset], 2))

    # ext_LUS C3 curve annotation text at end.
    anno3_text = r"$Diff:\ {0} \ km^{{2}}s^{{-2}}$".format(round(allc3[loiter]-allc3req[offset+loiter], 2))

    # Additional information.
    anno4_text = ("Payload: {0} kg \n".format(round(ext_LUS_payload,-2))
                 +"Engines: {0}\n".format(lus.Engines)
                 +"Initial T/W: {0}\n".format(round(lus.Thrust/lus.Stack_mass, 3))
                 +"Burns: {0} s + {1} s\n".format(round(br_time_1,2), round(br_time_2,2))
                 +"Parking orbit: {0} x {1} km\n".format(round(pery/1000,1), round(apo/1000,1))
                 +"Orbital period: {0} hour".format(round(transorb_period,2)) )

    # Name of chart file.
    filename = ("C3-{0}-[{1}, {2}]-".format(lus.Boiloff, offset, loiter)
               +"{0}kg-{1}-{1}kg_add".format(round(ext_LUS_payload,2), lus.Engines, lus.Stage_additions) )

    # Creating plot.
    plt.figure(figsize=(7.5, 4.5), dpi=180) # (11, 8) 80dpi dla full, (7.5, 4.5) 180dpi dla small


    # Plot ext_LUS C3
    plt.plot([swp_time(i) for i in newt],
             allc3,
             label="Instantaneous C3")

    # Plot required C3
    plt.plot([swp_time(i) for i in alltc3req],
             allc3req,
             label='Required C3')

    # Plot reentry speed barrier line
    plt.plot([swp_time(20)]*2,
             [38, 44],
             color='red',
             linewidth=1.5,
             linestyle="--")

    # Annotate reentry speed barrier line
    plt.annotate(r'Capsule Ventry barrier',
                 xy=(swp_time(20), 43.2),
                 xycoords='data',
                 xytext=(-150, +20),
                 textcoords='offset points',
                 fontsize=12,
                 bbox=dict(facecolor='white', edgecolor='None', alpha=0.65),
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2") )

    # Annotate ext_LUS C3 curve at start
    plt.annotate(anno2_text,
                 xy=(swp_time(newt[0]), allc3[0]),
                 xycoords='data',
                 xytext=(swp_time(8.5), 41.5),
                 textcoords='data',
                 fontsize=12,
                 bbox=dict(facecolor='white', edgecolor='None', alpha=0.65),
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.1") )

    # Annotate ext_LUS C3 curve at end
    plt.annotate(anno3_text,
                 xy=(swp_time(newt[loiter]), allc3[loiter] ),
                 xycoords='data',
                 xytext=(swp_time(15), 41.0),
                 textcoords='data',
                 fontsize=12,
                 bbox=dict(facecolor='white', edgecolor='None', alpha=0.65),
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.1") )

    # Put additional information in place
    plt.text(swp_time(0.6), 38.2,
             anno4_text,
             style='italic',
             fontsize=10,
             bbox={'facecolor':'white', 'alpha':0.5, 'pad':10} )

    # Add labels, tittle, legend
    #plt.xlabel('Date')
    plt.ylabel('C3 [$km^2s^{-2}$]')
    #plt.title(name)
    plt.legend(loc='upper left', prop={'size':12})

    # Axis and grid modyfications
    rc('xtick', labelsize=8)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator(minticks=6, maxticks=10))
    plt.gca().xaxis.set_minor_locator(mdates.DayLocator())
    plt.gcf().autofmt_xdate()

    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(1.0))
    #plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(0.1)

    plt.gca().xaxis.grid(True,'minor')
    #plt.gca().yaxis.grid(True,'minor')
    plt.gca().xaxis.grid(True,'major', linewidth=1)
    plt.gca().yaxis.grid(True,'major', linewidth=1)

    plt.tight_layout()

    if save:
        plt.savefig(filename+'.png', dpi=200)

    plt.show()


if __name__ == '__main__':
    offset = 16
    loiter = 4

    lus = ext_LUS('rl')
    lus.set_boiloff(0.002)

    all_the_rest(lus, offset, loiter, save = 0)
