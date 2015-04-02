# -*- coding: utf-8 -*-
from PyGMO import archipelago, problem, migration, population, island
from PyGMO.algorithm import jde,nsga_II, mde_pbx, de_1220
from PyGMO.topology import ring,fully_connected, unconnected, pan, erdos_renyi, ageing_clustered_ba, one_way_ring, clustered_ba
from PyKEP import planet_ss,epoch,epoch_from_string, planet
from datetime import datetime
from sys import stdout
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
pd.set_option('display.width', 250)

from planets import de421_planets
from mga import mga_1dsm


class Evolution_Handler(object):
    """
    Helper class for running multi and single-objective evolution and
    saving the results.
    """
    HDF_file = 'protolod.hdf'

    def run(self, seed='old'):
        """
        Runs evolution with earlier set parameters.

        'old' seed == use suitable chromosomes from `partials` dataframe,
                      generete new only if there's not enough old ones.
        'fresh' seed == evolve new chromosomes to use as initial guesses
        """
        if seed == 'old':
            founds, number_found = self.find_in_base()
            param = number_found - self.M_N

            if param < 0:
                print "We have only {0} usable chromosomes in the database, per {1} required.".format(number_found, self.M_N)
                l, __ = self.evolve_partials(abs(param))
                combined = founds+[l[i].x for i in range(len(l))]

            elif param > 0:
                combined = random.sample(founds, self.M_N)

            else:
                combined = founds

        if seed == 'fresh':
            print "Evolving fresh chromosomes..."
            l, __ = self.evolve_partials(self.M_N)
            combined = [l[i].x for i in range(len(l))]

        if len(combined) != self.M_N:      raise ValueError
        print "\nLaunching Multi-Objective evolution..."
        isl, prob = self.mlt_obj_evo(combined)
        self.writing_finals(isl, prob)

    def evolve_partials(self, quantity, autowrite = False):
        print "Launching Single-Objective evolution..."
        l, prob = self.sngl_obj_evo(quantity)
        self.writing_partials(l, prob, autowrite)
        return l, prob

    def probinit(self, aaa, n_obj):
        """Problem initialization

        self.probinit(algorithm name, n_obj -> 0 or 2)
        """
        # Set algorithm...
        if aaa == 'nsga':
            algo = nsga_II(m=0.05)
        else:
            algo = jde(memory=True)
            #algo = mde_pbx()
            #algo = de_1220()

        # ...and initialize problem with instance atributes
        prob = mga_1dsm(seq = self.FBseq,
                        multi_objective = n_obj,
                        dsm_dv_barrier = self.MAX_DV)

        prob.set_vinf((self.C3)**0.5)
        prob.set_tof(self.TOF[0], self.TOF[1])
        prob.set_entry_barrier(self.entry_barrier)
        prob.set_launch_window(self.EPOCHSTART, self.EPOCHEND)
        return prob, algo

    def find_in_base(self):
        """Find in a database chromosomes suitable for use as initial guesses
        """
        prt = pd.read_hdf(self.HDF_file, self.partials)

        ldate_mjd_df = prt['Ldate'].map(lambda x : epoch_from_string(x).mjd)
        redate_mjd_df = prt['REDate'].map(lambda x : epoch_from_string(x).mjd)
        tof_days = redate_mjd_df - ldate_mjd_df

        c = 4*(self.N-3)
        founds =  prt[ (prt['C3']          <=  self.C3)
                     & (ldate_mjd_df       >=  self.EPOCHSTART.mjd)
                     & (self.EPOCHEND.mjd  >=  ldate_mjd_df)
                     & (prt['VRE']         <=  self.entry_barrier)
                     & (prt['MDist']       >=  self.MARS_DIST/1000.)
                     & (prt['dVmag']       <=  self.MAX_DV)
                     & (tof_days           >=  self.TOF[0]*365.25)
                     & (tof_days           <=  self.TOF[1]*365.25) ].ix[:, :10+c].to_records(index=False)

        number_found = len(founds)
        return list(founds), number_found

    def sngl_obj_evo(self, lacking):
        """Single-objective (dV only) evolution
        `lacking` - number of chromosomes to return
        """
        prob, algo = self.probinit('jde', 0)
        l = list()
        u = 6+(self.N-3)*4
        for i in range(lacking):
            archi = archipelago(algo,prob,8,16, topology=fully_connected())
            for j in range(u):
                archi.evolve(5)
                stdout.write("\r{0} / {1}".format(i*u+j+1, lacking*u))
                stdout.flush()
            tmp = [isl for isl in archi]
            tmp.sort(key = lambda x: x.population.champion.f[0]);
            l.append(tmp[0].population.champion)
        stdout.write(" Done. ")
        return l, prob

    def mlt_obj_evo(self, l):
        """ewolucja multi-objective
        """
        prob, algo = self.probinit('nsga', 2) # paaaaanie, jak ewolujemy...
        sel = migration.best_s_policy(0.1, migration.rate_type.fractional)
        rep = migration.fair_r_policy(0.1, migration.rate_type.fractional)
        archi = archipelago(topology=clustered_ba())

        for g in range(self.M_N/2):
            pop = population(prob,14)
            pop.push_back(l[g])
            pop.push_back(l[g+self.M_N/2])
            isl = island(algo,pop,s_policy=sel, r_policy=rep)
            archi.push_back(isl)

        nn = 40
        for i in xrange(nn):
            archi.evolve(5)
            stdout.write("\r{0} / {1}".format(5*(i+1), 5*nn))
            stdout.flush()

        archi.join()
        stdout.write(" Done. Printing champion...\n\n")
        isl = min(archi, key=self.fitness_key)
        #print isl.population.champion.f
        return isl, prob

    def fitness_key(self, x):
        """Fitness key for multi objective optimalization solution selection
        """
        dv, c3, ventr = x.population.champion.f
        #print dv, c3, ventr
        if dv == 0 and ventr == 0:
            return -10000.0/c3
        else:
            return dv+ c3+ ventr

    def writing_partials(self, l, prob, autowrite):
        this_evo = pd.DataFrame()
        appended = 0

        for i in range(len(l)):
            chromosome = l[i].x
            desc = prob.AIO(chromosome, 0,0,1)
            if desc['dVmag'] <= self.MAX_DV:
                StoreVector = pd.DataFrame(desc, index = (0,))
                this_evo = this_evo.append(StoreVector, ignore_index=True)
                appended = appended + 1

        if appended == 0:
            writeornot = 'n'
        else:
            print "\n", this_evo.ix[:, 6+(self.N-2)*4:]
            print "\n", this_evo.ix[:, :6+(self.N-2)*4]
            if autowrite:
                writeornot = 'y'
            else:
                writeornot = raw_input('Save above results? (y/n) ')

        if writeornot == 'y':
            partials = pd.read_hdf(self.HDF_file, self.partials)
            partials = partials.append(this_evo, ignore_index=True)
            partials.to_hdf(self.HDF_file, self.partials)

    def writing_finals(self, isl, prob):
        chromosome = isl.population.champion.x
        desc = prob.AIO(chromosome, 0,0,1)
        desc['comm'] = '-'
        StoreVector = pd.DataFrame(desc, index = (0,))
        print StoreVector.ix[0:, 6+(self.N-2)*4:]

        writeornot = raw_input('Save above result to main dataframe/print chromosome/discard? (y/p/d) ')

        if writeornot == 'y':
            getcomm = raw_input('Comment (leave empy if none): ')
            if getcomm != '':
                StoreVector['comm'] = getcomm
            a = pd.read_hdf(self.HDF_file, self.finals)
            b = a.append(StoreVector, ignore_index=True)
            #print StoreVector.ix[0:, 10:].sort(columns=10)
            b.to_hdf(self.HDF_file, self.finals)

        elif writeornot == 'p':
            print StoreVector.ix[0:, 0:(6+(self.N-2)*4)]

    def df_precursor(self, dfname, addcomm = 0):
        l, prob = self.sngl_obj_evo(1)
        chromosome = l[0].x
        desc = prob.AIO(chromosome, 0,0,1)

        if addcomm:
            desc['comm'] = 'df_precursor'

        if desc['dVmag'] <= self.MAX_DV:
            StoreVector = pd.DataFrame(desc, index = (0,))
            StoreVector.to_hdf(self.HDF_file, dfname)
        else:
            print "dVmag exceeds the limit"
        print desc

    def main_num(self, main_num):
        self.M_N = main_num

    def c3(self, C3):
        self.C3 = C3

    def mars_dist(self, planet_dist):
        self.MARS_DIST = planet_dist

    def venus_dist(self, planet_dist):
        self.VENUS_DIST = planet_dist

    def max_dv(self, max_dv):
        self.MAX_DV = max_dv

    def v_entry(self, entry_barrier):
        self.entry_barrier = entry_barrier

    def launch_period(self, start, end):
        self.EPOCHSTART = epoch_from_string(start)
        self.EPOCHEND = epoch_from_string(end)

    def launch_period_mjd2000(self, start, end):
        self.EPOCHSTART = epoch(start, epoch.epoch_type.MJD2000)
        self.EPOCHEND = epoch(end, epoch.epoch_type.MJD2000)

    def tof(self, tof):
        """Set time of flight bonduaries"""
        self.TOF = tof

    def eme_500day_2018(self):
        """Define planets sequence with apropiate Mars flyby minimal distance"""
        a = de421_planets()
        self.FBseq = a.eme_500day_2018(self.MARS_DIST)
        self.N = len(self.FBseq)
        self.partials = 'eme_500day_2018_single'
        self.finals = 'eme_500day_2018_multi'

    def eme_700day_2018(self):
        """Define planets sequence with apropiate Mars flyby minimal distance"""
        a = de421_planets()
        self.FBseq = a.eme_700day_2018(self.MARS_DIST)
        self.N = len(self.FBseq)
        self.partials = 'eme_700day_2018_single'
        self.finals = 'eme_700day_2018_multi'

    def evme_2021(self):
        """Define planets sequence with apropiate Mars flyby minimal distance"""
        a = de421_planets()
        self.FBseq = a.evme_2021(self.VENUS_DIST, self.MARS_DIST)
        self.N = len(self.FBseq)
        self.partials = 'evme_2021_single'
        self.finals = 'evme_2021_multi'


class hdf_plot(object):
    """
    Helper class, for anal-yzing my shitty results.

    Keys used in ploting methods:
    'c3'  - C3 value
    'launch' - launch date
    'marsflyby' - Mars flyby date
    'venusflyby' - Venus flyby date
    'entry' - Earth reentry date
    'marsdist' - minimal Mars distance
    'venusdist' - minimal Venus distance
    'dec' - declination of a launch asymptohe ****
    'ra' - of a launch asymptote
    'ventr' - Earth reentry velocity
    """
    HDF_file = 'protolod.hdf'

    def __init__(self, dataframe):
        prt = pd.read_hdf(self.HDF_file, dataframe)
        self.prt = prt

        #dates needs a conversion from oryginal (something numpy's???)
        self.launch = prt.convert_objects(convert_dates='coerce')['Ldate'].astype(datetime)
        self.marsflyby = prt.convert_objects(convert_dates='coerce')['MFBDate'].astype(datetime)
        self.entry = prt.convert_objects(convert_dates='coerce')['REDate'].astype(datetime)

        try:
            self.venusflyby = prt.convert_objects(convert_dates='coerce')['VFBDate'].astype(datetime)
            self.vdist = prt['VDist']
        except:
            self.venusflyby = None
            self.vdist = None

        self.mdist = prt['MDist']
        self.c3 = prt['C3']
        self.dec = prt['DEC']
        self.ra = prt['RA']
        self.ventr = prt['VRE']

        self.xy = { 'c3': self.c3,
                    'launchd': self.launch,
                    'marsflybyd': self.marsflyby,
                    'venusflybyd': self.venusflyby,
                    'entryd': self.entry,
                    'marsdist': self.mdist,
                    'venusdist': self.vdist,
                    'dec': self.dec,
                    'ra': self.ra,
                    'ventr': self.ventr }

        self.datestart = 0
        self.dateend = 100000

    def plot(self, x, y, mX = 1, mY = 1):
        """Plots choosen properities of solutions in 2D.
        """
        x = self.xy[x]
        y = self.xy[y]

        # Find pareto front:
        xp, yp = self.pareto_frontier(x,y, mX, mY)

        plt.plot(x, y, 'o')
        plt.plot(xp, yp, '-')
        plt.tight_layout()
        plt.show()

    def fit_to_pareto(self, x, y, mX = 1, mY = 1, chart = True, order = 3, doublefront = None):
        """
        Fits polynomial to pareto front.
        Stupid, doesnt know where pareto front is... use mX/mY and check.
        At the moment it can take dates only as x variable, and convert
        it to mjd2000 epoch.
        """
        x = self.xy[x]
        y = self.xy[y]
        # Find pareto front:
        if doublefront == 'x':
            xp1, yp1 = self.pareto_frontier(x,y, mX, mY)
            xp2, yp2 = self.pareto_frontier(x,y, not mX, mY)
            xp = xp1 + xp2
            yp = yp1 + yp2
        else:
            xp, yp = self.pareto_frontier(x,y, mX, mY)

        # Conversion to float is required for fitting
        try:
            xp = [float(i) for i in xp]
            yp = [float(i) for i in yp]
        except:
            # Launch date needs additional handling
            xp = [epoch_from_string(i.strftime("%Y-%m-%d %H-%M-%S.%f")).mjd2000 for i in xp]
            yp = [float(i) for i in yp]

        # Sort, in case of broken order in doublefront
        myList = sorted([[xp[i], yp[i]] for i in range(len(xp))], reverse=0)
        xp = [pair[0] for pair in myList]
        yp = [pair[1] for pair in myList]

        x_to_fit = []
        y_to_fit = []
        # select fitting peroid
        for i in range(len(xp)):
            if self.dateend >= xp[i] >= self.datestart:
                x_to_fit.append(xp[i])
                y_to_fit.append(yp[i])

        z = np.polyfit(x_to_fit, y_to_fit, order) #fitnij polynomial
        f = np.poly1d(z) #zbuduj polynimoala
        values = f(x_to_fit)

        if order == 3:
            # https://stackoverflow.com/questions/24065904/numpy-calculate-polynom-efficiently
            fast_f = lambda x: z[3] + x*(z[2] + x*(z[1] + x*z[0]))
            fast_fitted = [fast_f(i) for i in x_to_fit]
            #print x_to_fit[0]
        else:
            fast_f = False
            fast_fitted = values

        if chart:
            plt.plot(x_to_fit, y_to_fit, 'o')
            plt.plot(x_to_fit, values, '--')
            plt.plot(x_to_fit, fast_fitted, '-')
            plt.tight_layout()
            plt.show()

        return f, fast_f   # returnij polynomiala

    def fit_period(self, fitstart, fitend):
        self.datestart = epoch_from_string(fitstart).mjd2000
        self.dateend = epoch_from_string(fitend).mjd2000

    def fit_period_mjd2000(self, fitstart, fitend):
        self.datestart = epoch(fitstart, epoch.epoch_type.MJD2000).mjd2000
        self.dateend = epoch(fitend, epoch.epoch_type.MJD2000).mjd2000

    def pareto_frontier(self, Xs, Ys, maxX = True, maxY = False):
        '''Method to take two equally-sized lists and return just the elements which lie
        on the Pareto frontier, sorted into order.
        Default behaviour is to find the maximum for both X and Y, but the option is
        available to specify maxX = False or maxY = False to find the minimum for either
        or both of the parameters.
        http://oco-carbon.com/metrics/find-pareto-frontiers-in-python/
        '''
        # Sort the list in either ascending or descending order of X
        myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
        # Start the Pareto frontier with the first value in the sorted list
        p_front = [myList[0]]
        # Loop through the sorted list
        for pair in myList[1:]:
            if maxY:
                if pair[1] >= p_front[-1][1]: # Look for higher values of Y…
                    p_front.append(pair) # … and add them to the Pareto frontier
            else:
                if pair[1] <= p_front[-1][1]: # Look for lower values of Y…
                    p_front.append(pair) # … and add them to the Pareto frontier
        # Turn resulting pairs back into a list of Xs and Ys
        p_frontX = [pair[0] for pair in p_front]
        p_frontY = [pair[1] for pair in p_front]
        return p_frontX, p_frontY


def populate_front():
    eh = Evolution_Handler()
    eh.mars_dist(100000) # m
    eh.venus_dist(400000) # m
    eh.main_num(5)
    eh.max_dv(0.01)
    eh.v_entry(12900)
    eh.tof([1.0, 2.2])
    eh.evme_2021()

    # launch strt
    start = epoch_from_string('2021-11-16 7:00:00').mjd2000

    # timespansd for fitting
    fitstart = epoch_from_string('2021-11-16 00:00:00').mjd2000
    fitend = epoch_from_string('2021-11-18 12:00:00').mjd2000
    h168 = 1.0/7
    h72 = 1.0/3
    h24 = 1
    h8 = 3
    h3 = 8
    h2 = 12
    iterations_vector = [ [h8, h3, 1, 5],] #,
    #                      [h24, h3, 1, 5] ]
    for s in iterations_vector:
        dt = 1.0/s[0]
        span = 1.0/s[1]
        for i in range(int(2*1/dt)):
            a = hdf_plot('evme_2021_single')
            a.fit_period_mjd2000(fitstart, fitend)
            c4fun, __ = a.fit_to_pareto('launchd', 'c3', mX=1, mY=0, chart = 0, order = s[3], doublefront = 'x')
            print c4fun(start + i*dt + 0.0*span), epoch(start + i*dt, epoch.epoch_type.MJD2000)

            eh.launch_period_mjd2000(start + i*dt,
                                     start + i*dt + span)
            eh.c3( c4fun(start + i*dt + 0.5*span) +0.00)
            eh.evolve_partials(s[2], autowrite = 1)
        LL = hdf_plot('evme_2021_single')
        LL.plot('launchd', 'c3', mX=1, mY=0)
        LL.fit_to_pareto('launchd', 'c3', mX=1, mY=0, chart = 1, order = 3, doublefront = 'x')


if __name__ == '__main__':
    # ustawienia
    eh = Evolution_Handler()
    eh.tof([1.0, 2.2])
    eh.mars_dist(100000) # m
    eh.venus_dist(400000) # m
    eh.c3(25.6) # km^2/s^2
    eh.main_num(32)
    eh.max_dv(0.01)
    eh.v_entry(13900)

    #eh.launch_period('2017-12-25 02:00:01.000',
    #                 '2018-01-05 23:30:01.000')
    #eh.launch_period('2018-05-02 00:00:01.000',
    #                 '2018-05-02 03:00:01.000')
    #eh.launch_period('2021-11-16 00:00:01.000',
    #                 '2021-11-17 00:00:01.000')

    #eh.evme_2021()
    #eh.eme_500day_2018()
    #eh.eme_700day_2018()

    #eh.evolve_partials(1)
    #eh.run(seed='old')
    #eh.df_preGGGGGcursor(eh.finals, addcomm=1)

    #populate_front()

    if 1:
        a = hdf_plot('evme_2021_single')
        start = epoch_from_string('2021-11-17 00:00:00').mjd2000
        end = epoch_from_string('2021-11-29 12:00:00').mjd2000
        a.fit_period_mjd2000(start, end)
        #print type(a.launch[0])
        #a.plot('launchd', 'c3', mX=1, mY=0)
        #c4fun, __ = a.fit_to_pareto('launchd', 'c3', mX=1, mY=0, chart = 1, order = 4, doublefront = 'x')
        #print c4fun(end-1-1.0/24)
        #a.plot('c3', 'entry', mX=0, mY=0)
        #a.plot('ventr', 'c3', mX=0, mY=0)
        #a.plot('ventr', 'launch', mX=1, mY=1)
        #a.plot('launch', 'ra', mX=0, mY=1)
        #a.plot('launch', 'dec', mX=0, mY=1)

    # .astype(datetime)
    # onast = pandastests.a.convert_objects(convert_dates='coerce')[18]
    # 57 stopni
