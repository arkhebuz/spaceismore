# -*- coding: utf-8 -*-
from numpy import matrix, cos, sin
from PyKEP import DEG2RAD


class pos(object):
    def __init__(self):
        self.trajectory_positions = []

        axtl = 23.43929*DEG2RAD # Earth axial tlit
        # Macierz przejscia z  ECI do helio
        self.mactran = matrix([ [1, 0, 0],
                                [0, cos(axtl), sin(axtl)],
                                [0, -sin(axtl), cos(axtl)] ])
        self.day = 86400.0

    def positions_lambert(self, l, sol=0, units = 1.0, index = 0):
        """
        Plots a particular solution to a Lambert's problem

        USAGE: plot_lambert(ax,l, N=60, sol=0, units = 'PyKEP.AU', legend = 'False')
          * l:      PyKEP.lambert_problem object
          * sol:    solution to the Lambert's problem we want to plot (must be in 0..Nmax*2)
                where Nmax is the maximum number of revolutions for which there exist a solution.
          * units:  the length unit to be used in the plot
        """
        from PyKEP import propagate_lagrangian, AU

        if sol > l.get_Nmax()*2:
            raise ValueError("sol must be in 0 .. NMax*2 \n * Nmax is the maximum number of revolutions for which there exist a solution to the Lambert's problem \n * You can compute Nmax calling the get_Nmax() method of the lambert_problem object")

        #We extract the relevant information from the Lambert's problem
        r = l.get_r1()
        v = l.get_v1()[sol]
        T = l.get_tof()
        mu = l.get_mu()

        #We define the integration time ...
        if T/86400. < 1:
            N = int(T) #...compute number of points...
            dt = T / (N-1.)
        else:
            N = int(2*T/86400.) #...compute number of points...
            dt = T / (N-1.)
        timetuple = [i*dt for i in range(N)]

        #... and alocate the cartesian components for r
        x = [0.0]*N
        y = [0.0]*N
        z = [0.0]*N

        #We calculate the spacecraft position at each dt
        for i in range(N):
            x[i] = r[0]/units
            y[i] = r[1]/units
            z[i] = r[2]/units
            r,v = propagate_lagrangian(r,v,dt,mu)

        self.trajectory_positions.append([index, x, y, z, timetuple])

    def positions_kepler(self, r,v,t,mu, units = 1, index = 0):
        """
        Plots the result of a keplerian propagation

        USAGE: plot_kepler(ax,r,v,t,mu, N=60, units = 1, color = 'b', legend = False):
          * r:      initial position (cartesian coordinates)
          * v:      initial velocity (cartesian coordinates)
          * t:      propagation time
          * mu:     gravitational parameter
          * units:  the length unit to be used in the plot
        """

        from PyKEP import propagate_lagrangian

        #We define the integration time ...
        if t/86400. < 1:
            N = int(t) #...compute number of points...
            dt = t / (N-1.)
        else:
            N = int(2*t/86400.) #...compute number of points...
            dt = t / (N-1.)
        timetuple = [i*dt for i in range(N)]

        #... and calcuate the cartesian components for r
        x = [0.0]*N
        y = [0.0]*N
        z = [0.0]*N

        #We calculate the spacecraft position at each dt
        for i in range(N):
            x[i] = r[0]/units
            y[i] = r[1]/units
            z[i] = r[2]/units
            r,v = propagate_lagrangian(r,v,dt,mu)

        self.trajectory_positions.append([index, x, y, z, timetuple])

    def return_sc_positions(self):
        return self.trajectory_positions

    def set_launch_epoch(self, mjd2000_epoch):
        self.Lepoch = mjd2000_epoch

    def rework_time(self):
        T = self.trajectory_positions
        n = len(T)

        ep = self.Lepoch
        timespans = [T[i][4][-1]/86400. for i in range(n)]
        beginings = [ep + sum(timespans[:i]) for i in range(n)]

        # initialize lists...
        t = []
        x = []
        y = []
        z = []

        # ...and join values, correcting for time
        for g in range(n):
            t = t + [beginings[g]+i/86400. for i in T[g][4]]
            x = x + T[g][1]
            y = y + T[g][2]
            z = z + T[g][3]

        # save spacecraft state
        self.sc_state = [x, y, z, t]

    def eq2eclipt(self, xyz):
        macxyz = matrix(xyz)
        return self.mactran.dot(macxyz)

    def planets_pos(self):
        from jplephem import Ephemeris
        import de421
        from PyKEP import epoch
        self.eph = Ephemeris(de421)

        earthpos = []
        marspos = []
        venuspos = []

        for ep in self.sc_state[3]:
            posSun, __ = self.eph.position_and_velocity('sun', epoch(ep, epoch.epoch_type.MJD2000).jd)

            positione, __ = self.eph.position_and_velocity('earthmoon', epoch(ep, epoch.epoch_type.MJD2000).jd)
            positione = self.eq2eclipt(positione - posSun)
            earthpos.append(positione)

            positionm, __ = self.eph.position_and_velocity('mars', epoch(ep, epoch.epoch_type.MJD2000).jd)
            positionm = self.eq2eclipt(positionm - posSun)
            marspos.append(positionm)

            positionv, __ = self.eph.position_and_velocity('venus', epoch(ep, epoch.epoch_type.MJD2000).jd)
            positionv = self.eq2eclipt(positionv - posSun)
            venuspos.append(positionv)

        self.earthpos_km = [earthpos[i].reshape((1,3)).tolist()[0] for i in range(len(earthpos))]
        self.marspos_km = [marspos[i].reshape((1,3)).tolist()[0] for i in range(len(earthpos))]
        self.venuspos_km = [venuspos[i].reshape((1,3)).tolist()[0] for i in range(len(earthpos))]

    def distance_to_planets(self):
        dre = []
        drm = []
        drv = []

        for i in range(len(self.sc_state[3])):
            ep = self.earthpos_km[i]
            dist = [ep[0] - self.sc_state[0][i] / 1000.,
                    ep[1] - self.sc_state[1][i] / 1000.,
                    ep[2] - self.sc_state[2][i] / 1000.]
            dre.append((sum([g**2 for g in dist]))**.5)

            mp = self.marspos_km[i]
            dist = [mp[0] - self.sc_state[0][i] / 1000.,
                    mp[1] - self.sc_state[1][i] / 1000.,
                    mp[2] - self.sc_state[2][i] / 1000.]
            drm.append((sum([g**2 for g in dist]))**.5)

            vp = self.venuspos_km[i]
            dist = [vp[0] - self.sc_state[0][i] / 1000.,
                    vp[1] - self.sc_state[1][i] / 1000.,
                    vp[2] - self.sc_state[2][i] / 1000.]
            drv.append((sum([g**2 for g in dist]))**.5)

        return dre,drm, drv, self.sc_state[3]


def pl3_plot(earth, mars, venus, time, filename = False):
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime, timedelta
    from PyKEP import AU
    AU = AU/1000

    time = [datetime(2000, 01, 01) +timedelta(x) for x in time]

    plt.plot(time, [i/AU for i in earth], label='Earth')
    plt.plot(time, [i/AU for i in venus], label='Venus')
    plt.plot(time, [i/AU for i in mars], label='Mars')

    plt.legend(loc='upper left', prop={'size':13})
    plt.ylabel('Distance to spacecraft [$AU$]')

    # Axis and grid modyfications
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator(minticks=6, maxticks=10))
    plt.gca().xaxis.set_minor_locator(mdates.MonthLocator())
    plt.gcf().autofmt_xdate()

    #plt.gca().yaxis.set_major_locator(plt.MultipleLocator(1.0))
    #plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(0.1)

    plt.gca().xaxis.grid(True,'minor')
    #plt.gca().yaxis.grid(True,'minor')
    plt.gca().xaxis.grid(True,'major', linewidth=1)
    plt.gca().yaxis.grid(True,'major', linewidth=1)

    plt.tight_layout()
    if filename:
        plt.savefig(filename+'.png', dpi=200)
    plt.show()


if __name__ == '__main__':
    from planets import de421_planets
    from mga import mga_1dsm
    import pandas as pd

    prt = pd.read_hdf('protolod.hdf', 'eme_500day_2018_single')
    prt = pd.read_hdf('protolod.hdf', 'eme_700day_2018_single')
    prt = pd.read_hdf('protolod.hdf', 'evme_2021_single')

    a = de421_planets()
    FBseq = a.evme_2021(400*1000,100*1000)
    prob = mga_1dsm(seq = FBseq)

    b = prt.ix[139:139, :]
    gen = b.values.tolist()[0][:14]

    posinst = pos()
    posinst = prob.AIO(gen, 0,0,0, posinst)

    posinst.rework_time()
    posinst.planets_pos()
    drz, drm, drv, t = posinst.distance_to_planets()
    pl3_plot(drz, drm, drv, t, filename = 0)


