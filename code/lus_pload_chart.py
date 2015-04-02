# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from scipy.optimize import fmin as simplex

from burn_sym import *
from stage import ext_LUS


C3_target = 20
C3s = range(-2,61,2)
rl10_ext_LUS_payload = [71817.069724202156, 68288.214755058289, 64989.727959036827, 61900.086379051208, 59000.368613004684, 56273.861977458,
                        53705.765506625175, 51282.913953065872, 48993.559885025024, 46827.188131213188, 44774.358066916466, 42826.569959521294,
                        40976.150304079056, 39216.153654456139, 37540.279501676559, 35942.797440290451, 34418.484646081924, 32962.571239471436,
                        31570.691925287247, 30238.843989372253, 28963.350605964661, 27740.82792699337, 26568.156999349594, 25442.457854747772,
                        24361.067381501198, 23321.519047021866, 22321.524953842163, 21358.960309624672, 20431.849265098572, 19538.351958990097,
                        18676.752689480782, 17845.450016856194]
mb60_ext_LUS_payload = [76371.628886461258, 72677.338966727257, 69224.009478092194, 65989.132878184319, 62952.919372916222, 60097.897803783417,
                        57408.586460351944, 54871.217080950737, 52473.504316806793, 50204.451325535774, 48054.18489575386, 46013.815757632256,
                        44075.31860768795, 42231.429916620255, 40475.559255480766, 38801.713496446609, 37204.429864883423, 35678.719773888588,
                        34220.017164945602, 32824.135619401932, 31487.22892999649, 30205.758190155029, 28976.460814476013, 27796.324306726456,
                        26662.562844157219, 25572.596919536591, 24524.032339453697, 23514.647108316422, 22542.373886704445, 21605.288082361221,
                        20701.594975590706, 19829.619097709656]

#~ lus = ext_LUS('rl')
#~ lus.set_boiloff(0.002)
#~ lus.set_launch_date(0)
#~ lus.set_time(4)
#~ lus.print_info()
#~ ext_LUS_payload = [simplex(c3, 27000, args=(lus, C3_target))[0] for C3_target in C3s]
#~ print ext_LUS_payload


plt.figure(figsize=(11, 8), dpi=80)
plt.plot(C3s, [i/1000 for i in rl10_ext_LUS_payload], label="RL-10 ext-LUS")
plt.plot(C3s, [i/1000 for i in mb60_ext_LUS_payload], label="MB-60 ext-LUS")

plt.xlabel('C3 [$km^2s^{-2}$]')
plt.ylabel('Payload mass [$metric\ tons$]')
plt.legend(loc='upper right', prop={'size':13})

plt.text(0, 12.2,
         "Figure assumes circular, 240 km high parking orbit with an inclination of 28.5 deg",
         style='italic',
         fontsize=10,
         bbox={'facecolor':'white', 'alpha':0.5, 'pad':10} )

plt.gca().xaxis.set_major_locator(plt.MultipleLocator(4))
plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))

plt.gca().yaxis.set_major_locator(plt.MultipleLocator(5))
plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(2.5))

plt.gca().xaxis.grid(True,'minor')
plt.gca().yaxis.grid(True,'minor')
plt.gca().xaxis.grid(True,'major', linewidth=1)
plt.gca().yaxis.grid(True,'major', linewidth=1)

plt.ylim(10, 80)
plt.xlim(-4, 64)

plt.tight_layout()
#plt.savefig('payload_chart.png', dpi=200)
plt.show()
