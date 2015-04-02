# -*- coding: utf-8 -*-
from scipy.optimize import fmin as simplex
from PyKEP import epoch,epoch_from_string
import pandas as pd

from planets import de421_planets
from mga import mga_1dsm
from burn_sym import *
from stage import ext_LUS
from trajectories import hdf_plot


# Read partials dataframe
prt = pd.read_hdf('protolod.hdf', 'eme_500day_2018_single')
fin = pd.read_hdf('protolod.hdf', 'eme_500day_2018_multi')
prtv = pd.read_hdf('protolod.hdf', 'evme_2021_single')
finv = pd.read_hdf('protolod.hdf', 'evme_2021_multi')

#
a_prt = hdf_plot('eme_500day_2018_single')
a_fin = hdf_plot('eme_500day_2018_multi')
a_prtv = hdf_plot('evme_2021_single')
a_finv = hdf_plot('evme_2021_multi')


# Introduce lus to global namespace
print "lusrl"
lusrl = ext_LUS('rl')
lusrl.set_boiloff(0.002)
lusrl.set_launch_date(0)
lusrl.set_time(4)
#lusrl.print_info()

print "lusmb"
lusmb = ext_LUS('mb')
lusmb.set_boiloff(0.002)
lusmb.set_launch_date(0)
lusmb.set_time(4)
#lusmb.print_info()

print "lusj"
lusj = ext_LUS('j')
lusj.set_boiloff(0.002)
lusj.set_launch_date(0)
lusj.set_time(4)
#lusj.print_info()


# Initialize problem (to use its description/plot functions)
a = de421_planets()
FBseq = a.eme_500day_2018(100*1000)
FBseq_v = a.evme_2021(400*1000, 100*1000)
prob = mga_1dsm(seq = FBseq)
probv = mga_1dsm(seq = FBseq_v)


# Warp optimalization routine in a more function
def pload(C3_target, lus = lusrl):
    return simplex(c3, 27000, args=(lus, C3_target))[0]


if __name__ == '__main__':
    C3 = 39
    ldate_mjd2000_df = prt['Ldate'].map(lambda x : epoch_from_string(x).mjd2000)
    redate_mjd2000_df = prt['REDate'].map(lambda x : epoch_from_string(x).mjd2000)
    tof_days = redate_mjd2000_df - ldate_mjd2000_df
    print ldate_mjd2000_df[0]

    print type(list(prt[ (prt['C3'] < 38.2)
             & (6600 > ldate_mjd2000_df)
             & (ldate_mjd2000_df > 6500.)
             & (prt['VRE'] < 14200)
             & (prt['MDist'] >= 100)
             & (prt['dVmag'] < 0.01)
             & (tof_days < 1.6*365.25)].ix[:, :10].to_records(index=False)))

    newdf = pd.DataFrame()
    for i in range(prtv.shape[0]):
        # ...extract data as (hopefully) python's list...
        vector = prtv.ix[i:i, :].to_records(index=False)[0]
        opis = probv.AIO(list(vector)[:14], doplot=0, doprint=0)
        #opis['comm'] = '-'

        #~ if opis['dVmag'] != vector[-1]:
            #~ print i, opis['dVmag'], vector[-1], opis['Ldate'], 'drh'

        if opis['dVmag'] < 0.01:
            #print i, opis['dVmag'], vector[-1], opis['Ldate']
            StoreVector = pd.DataFrame(opis, index = (0,))
            newdf = newdf.append(StoreVector, ignore_index=True)
    print newdf
    #newdf.to_hdf("protolod.hdf", "evme_2021_single")

    tdtbhbttbbb = ("""
    # if launched as a script, do...
    start = epoch_from_string('2018-01-01 00:00:00').mjd2000
    end = epoch_from_string('2018-01-06 00:00:00').mjd2000
    a = hdf_plot('partials')
    a.fit_period_mjd2000(start, end)
    __, c3fun = a.fit_to_pareto('launch', 'c3', mX=0, mY=0, chart = False)

    t0 = epoch_from_string('2018-01-05 00:00:00').mjd2000
    C3_target = c3fun(t0)

    print "Launching numerical simualtion..."
    ext_LUS_payload = pload(C3_target)
    print "         Payload: {0} kg @ {1}loi ".format(round(ext_LUS_payload,2), 4)

    # Get parking orbit parameters
    __, __, pery, apo, transorb_period, br_time_1, br_time_2 = c3(ext_LUS_payload, lusrl, target_C3 = C3_target, rich_return = 1)
    """)

# b = prt[prt[11]> 180].sort([10], ascending=1).iloc[0]
#  prob.pretty(b.values.tolist()[:10])

