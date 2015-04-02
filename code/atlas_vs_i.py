# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin as simplex

from burn_sym import *
from stage import ext_LUS


centaur_dry = 2138 # kg
centaur_prop = 20830 # kg

# Data, from ATLAS LAUNCH SYSTEM MISSION PLANNER’S GUIDE, Revision 10a, s ~130
CCAFS_286deg_200km = 20500 # kg
CCAFS_516deg_200km = 19000 # kg
VAFB_634deg_200km = 17255 # kg
inclinations = [28.6, 51.6, 63.4]
payloads = [CCAFS_286deg_200km, CCAFS_516deg_200km, VAFB_634deg_200km]
# 12120 kg for CCAFS 63.4deg -> heavy dogleg

payloads = [i+centaur_dry for i in payloads]
payloads = [i*1./payloads[0] for i in payloads]
x_to_plot = range(28, 65)

#fitnij polynomial
z = np.polyfit(inclinations, payloads, 2)
fof = np.poly1d(z) #zbuduj polynimoala

plt.figure(figsize=(10, 6), dpi=80)
plt.plot(x_to_plot, fof(x_to_plot), label='$2^{nd}$ order polynomial')

plt.text(29.5, 0.983,
         "CCAFS \n20500 kg at 28.6 deg",
         style='italic',
         fontsize=10,
         horizontalalignment='center',
         bbox={'facecolor':'white', 'alpha':0.5, 'pad':5} )

plt.text(55, 0.94,
         "CCAFS \n19000 kg at 51.6 deg",
         style='italic',
         fontsize=10,
         horizontalalignment='center',
         bbox={'facecolor':'white', 'alpha':0.5, 'pad':5} )

plt.text(59, 0.845,
         "VAFB \n17255 kg at 63.4 deg",
         style='italic',
         fontsize=10,
         horizontalalignment='center',
         bbox={'facecolor':'white', 'alpha':0.5, 'pad':5} )

plt.plot(inclinations, payloads, 'o', label='Mission Planner\'s Guide values')
plt.xlabel('Inclination of a 200 km high LEO orbit [$deg$]')
plt.ylabel('Normalized Atlas 552 IMLEO mass (to $1$ at $28.6 deg$)')
plt.legend(loc='lower left', prop={'size':12})

plt.gca().xaxis.set_major_locator(plt.MultipleLocator(5))
plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))
plt.gca().yaxis.set_major_locator(plt.MultipleLocator(.05))
plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(.01))
plt.gca().xaxis.grid(True,'minor')
#plt.gca().yaxis.grid(True,'minor')
plt.gca().xaxis.grid(True,'major', linewidth=1)
plt.gca().yaxis.grid(True,'major', linewidth=1)

print fof(56.9)

plt.tight_layout()
#plt.savefig('atlas5.png', dpi=200)
plt.show()


if 1:
    C3 = 19.79 #km2s2
    DEC = 56.917 # deg

    lus = ext_LUS('mb')
    lus.set_boiloff(0.002)
    lus.set_time(4)
    lus.compute_stack_mass()
    lus.hit_from_incl(fof(DEC))
    #lus.print_info()
    pl = simplex(c3, 37000, args=(lus, C3))[0]
    print "2x MB-60 ext-LUS:", round(pl, -2)

    lus2 = ext_LUS('rl')
    lus2.set_boiloff(0.002)
    lus2.set_time(4)
    lus2.compute_stack_mass()
    lus2.hit_from_incl(fof(DEC))
    #lus.print_info()
    pl = simplex(c3, 37000, args=(lus2, C3))[0]
    print "2x RL-10 ext-LUS:", round(pl, -2)


