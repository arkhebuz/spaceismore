from PyGMO.problem import base as base_problem
from PyKEP import epoch,DAY2SEC,planet_ss,MU_SUN,lambert_problem,propagate_lagrangian,fb_prop, AU, RAD2DEG, DEG2RAD
from math import pi, cos, sin, acos, sqrt, asin
from scipy.linalg import norm


"""
skopiowane z /usr/local/lib/python2.7/dist-packages/PyKEP/interplanetary/_mga_1dsm.py (oryginal)

This class represents an mga-1DSM global optimization problem (single and multi-objective)

SEE : Izzo: "Global Optimization and Space Pruning for Spacecraft Trajectory Design, Spacecraft Trajectory Optimization, Conway, B. (Eds.), Cambridge University Press, pp.178-199, 2010)

The decision vector is [t0,u,v,Vinf,eta1,T] + [beta, rp/rV, eta2,a2] ..... in the units: [mjd2000,nd,nd,km/s,nd,years] + [rad,nd,nd,nd] + ....
where Vinf = Vinf_mag*(cos(theta)*cos(phi)i+cos(theta)*sin(phi)j+sin(phi)k) and theta = 2*pi*u and phi = acos(2*v-1)-pi/2

Each leg time-of-flight can be obtained as Tn = T*an, T(n-1) = (T - Tn)*a(n-1), .... , Ti = (T-T(i+1)-T(i+2)- .... - Tn)*ai

NOTE: The resulting problem is box-bounded (unconstrained). The resulting trajectory is time-bounded.

"""


def f(x):
    """Return number of non zero elements in three elements vector"""
    if x[0] == 0: a = 0
    else: a = 1
    #if x[1] == 0: pass
    #else: a = a + 1
    if x[2] == 0: pass
    else: a = a + 1
    return a


class mga_1dsm(base_problem):
    """
    Constructs a global optimization problem (box-bounded, continuous) representing an interplanetary trajectory modelled
    as a Multiple Gravity Assist trajectory that allows one only Deep Space Manouvre between each leg.

    USAGE: opt_prob = mga_1dsm(seq = [planet_ss('earth'),planet_ss('venus'),planet_ss('earth')], t0 = [epoch(0),epoch(1000)], tof = [1.0,5.0], vinf = 2.5, multi_objective = False)

    * seq:  list of PyKEP.planet defining the encounter sequence for the trajectoty (including the initial planet)
    * t0:   list of PyKEP epochs defining the launch window
    * tof:  minimum and maximum time of flight allowed (in years)
    * vinf: maximum launch hyperbolic velocity allowed
    * multi-objective: when True defines the problem as a multi-objective problem, returning total DV and time of flight
    """
    def __init__(self, seq = [planet_ss('earth'),planet_ss('mars'),planet_ss('earth'),planet_ss('mars')], t0 = [epoch(0),epoch(1000)], tof = [1.0,5.0], vinf = 2.5, multi_objective = False, dsm_dv_barrier = 20):
        self.__n = len(seq) - 1
        dim = 6 + (self.__n-1) * 4
        obj_dim = multi_objective + 1 #tutaj zmienia sie wymiar !!!!
        #First we call the constructor for the base PyGMO problem
        #As our problem is n dimensional, box-bounded (may be multi-objective), we write
        #(dim, integer dim, number of obj, number of con, number of inequality con, tolerance on con violation)
        super(mga_1dsm,self).__init__(dim,0,obj_dim,0,0,0)

        #We then define all planets in the sequence as data members
        self.seq = seq

        #moje
        self.dsm_dv_barrier = dsm_dv_barrier

        #And we compute the bounds
        lb = [t0[0].mjd2000,0.0,0.0,1e-5     ,0.0     ,tof[0]*365.25] + [0   ,1.05 ,1e-5    ,1e-5]     * (self.__n-1)
        ub = [t0[1].mjd2000,1.0,1.0,vinf*1000,1.0-1e-5,tof[1]*365.25] + [2*pi,30.0,1.0-1e-5,1.0-1e-5] * (self.__n-1)

        #Accounting that each planet has a different safe radius......
        for i,pl in enumerate(seq[1:-1]):
            lb[7+4*i] = pl.safe_radius / pl.radius

        #And we set them
        self.set_bounds(lb,ub)

    #Objective function
    def _objfun_impl(self,x):
        #1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
        T = list([0]*(self.__n))
        #a[-i] = x[-1-(i-1)*4]
        for i in xrange(self.__n-1):
            j = i+1;
            T[-j] = (x[5] - sum(T[-(j-1):])) * x[-1-(j-1)*4]
        T[0] = x[5] - sum(T)

        #2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.__n+1))
        r_P = list([None] * (self.__n+1))
        v_P = list([None] * (self.__n+1))
        DV = list([None] * (self.__n+1))

        for i,planet in enumerate(self.seq):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i],v_P[i] = self.seq[i].eph(t_P[i])

        #3 - We start with the first leg
        theta = 2*pi*x[1]
        phi = acos(2*x[2]-1)-pi/2

        Vinfx = x[3]*cos(phi)*cos(theta)
        Vinfy = x[3]*cos(phi)*sin(theta)
        Vinfz = x[3]*sin(phi)

        v0 = [a+b for a,b in zip(v_P[0], [Vinfx,Vinfy,Vinfz])]
        r,v = propagate_lagrangian(r_P[0], v0, x[4]*T[0]*DAY2SEC, MU_SUN)

        #Lambert arc to reach seq[1]
        dt = (1-x[4])*T[0]*DAY2SEC
        l = lambert_problem(r,r_P[1],dt,MU_SUN)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        #First DSM occuring at time nu1*T1
        DV[0] = norm([a-b for a,b in zip(v_beg_l, v)])

        #4 - And we proceed with each successive leg
        for i in range(1,self.__n):
            #Fly-by
            v_out = fb_prop(v_end_l, v_P[i] ,x[7+(i-1)*4]*self.seq[i].radius, x[6+(i-1)*4], self.seq[i].mu_self)
            #s/c propagation before the DSM
            r,v = propagate_lagrangian(r_P[i], v_out, x[8+(i-1)*4]*T[i]*DAY2SEC, MU_SUN)
            #Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1-x[8+(i-1)*4])*T[i]*DAY2SEC
            l = lambert_problem(r, r_P[i+1], dt, MU_SUN)
            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]
            #DSM occuring at time nu2*T2
            DV[i] = norm([a-b for a,b in zip(v_beg_l, v)])

        #Last Delta-v
        DV[-1] = norm([a-b for a,b in zip(v_end_l,v_P[-1])])

        #moje
        Eh56500 = 61933310.95
        #Eh100000 = 61517435.56
        entry_Vel = (DV[-1]**2 + 2*Eh56500)**0.5
        entry_Vel2 = entry_Vel
        #axtl = 23.43929*DEG2RAD
        #DEC = abs(asin( sin(axtl)*cos(phi)*sin(theta) + cos(axtl)*sin(phi) ))*RAD2DEG # deklinacja asymptoty ucieczkowej
        sum_dv = sum(DV[:-1])
        eff_C3 = (x[3])**2

        if entry_Vel < self.entry_vel_barrier:
            entry_Vel2 = 0.0
            del DV[-1]

        #~ if eff_C3 < self.C3_barrier:
            #~ #eff_C3 = 0
            #~ pass

        if sum_dv < self.dsm_dv_barrier:
            sum_dv = 0+entry_Vel2
        else:
            sum_dv = sum(DV)

        if self.f_dimension == 1:
            return (sum_dv,)
        else:
            return (sum_dv, eff_C3, entry_Vel2) #,

    def _compare_fitness_impl(self, f1, f2):
        if self.f_dimension == 1:   return f1[0] < f2[0]
        elif self.f_dimension == 3:
            a = f(f1)
            b = f(f2)
            if a < b:   return True
            elif a == b:
                #if a == 1:  return f1[1] < f2[1]
                #else:
                c = (f1[0]/self.dsm_dv_barrier) + (f1[2]/self.entry_vel_barrier)
                d = (f2[0]/self.dsm_dv_barrier) + (f2[2]/self.entry_vel_barrier)
                return c < d
            else:   return False

    def AIO(self,x, doplot = True, doprint = True, rtrn_desc = True, dists_class = None):
        P = doprint
        plots_datas = list([None] * (self.__n))
        PLDists = list([None] * (self.__n-1))
        FBDates = list([None] * (self.__n-1))

        #1 -  we 'decode' the chromosome recording the various times of flight (days) in the list T
        T = list([0]*(self.__n))
        #a[-i] = x[-1-(i-1)*4]
        for i in xrange(self.__n-1):
            j = i+1;
            T[-j] = (x[5] - sum(T[-(j-1):])) * x[-1-(j-1)*4]
        T[0] = x[5] - sum(T)

        #2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.__n+1))
        r_P = list([None] * (self.__n+1))
        v_P = list([None] * (self.__n+1))
        DV = list([None] * (self.__n+1))

        for i,planet in enumerate(self.seq):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i],v_P[i] = self.seq[i].eph(t_P[i])

        #3 - We start with the first leg
        if P: print "First Leg:        " + self.seq[0].name + " to " + self.seq[1].name

        theta = 2*pi*x[1]
        phi = acos(2*x[2]-1)-pi/2

        Vinfx = x[3]*cos(phi)*cos(theta)
        Vinfy = x[3]*cos(phi)*sin(theta)
        Vinfz = x[3]*sin(phi)

        if P:
            print("Departure:        " + str(t_P[0]) + " (" + str(t_P[0].mjd2000) + " mjd2000) \n"
                  "Duration:         " + str(T[0]) + " days\n"
                  "VINF:             " + str(x[3] / 1000) + " km/sec\n"
                  "C3:               " + str((x[3] / 1000)**2) + " km^2/s^2")

        v0 = [a+b for a,b in zip(v_P[0], [Vinfx,Vinfy,Vinfz])]
        r,v = propagate_lagrangian(r_P[0], v0, x[4]*T[0]*DAY2SEC, MU_SUN)
        if P: print "DSM after         " + str(x[4]*T[0]) + " days"

        #Lambert arc to reach seq[1]
        dt = (1-x[4])*T[0]*DAY2SEC
        l = lambert_problem(r, r_P[1], dt, MU_SUN)
        v_end_l = l.get_v2()[0]
        v_beg_l = l.get_v1()[0]

        # Append data needed for potential plot generation
        plots_datas[0] = [v0, x[4]*T[0]*DAY2SEC, l]

        #First DSM occuring at time nu1*T1
        DV[0] = norm([a-b for a,b in zip(v_beg_l,v)])
        if P: print "DSM magnitude:    " + str(DV[0]) + "m/s"

        #4 - And we proceed with each successive leg
        for i in range(1,self.__n):
            if P:
                print("\nleg no.           " + str(i+1) + ": " + self.seq[i].name + " to " + self.seq[i+1].name + "\n"
                      "Duration:         " + str(T[i]) + "days")

            #Fly-by
            v_out = fb_prop(v_end_l, v_P[i], x[7+(i-1)*4]*self.seq[i].radius, x[6+(i-1)*4], self.seq[i].mu_self)
            PLDists[i-1] = (x[7+(i-1)*4] -1)*self.seq[i].radius/1000.
            FBDates[i-1] = t_P[i]

            if P:
                print("Fly-by epoch:     " + str(t_P[i]) + " (" + str(t_P[i].mjd2000) + " mjd2000) \n"
                      "Fly-by radius:    " + str(x[7+(i-1)*4]) + " planetary radii\n"
                      "Fly-by distance:  " + str( (x[7+(i-1)*4] -1)*self.seq[i].radius/1000.) + " km")

            #s/c propagation before the DSM
            r,v = propagate_lagrangian(r_P[i], v_out, x[8+(i-1)*4]*T[i]*DAY2SEC, MU_SUN)
            if P: print "DSM after         " + str(x[8+(i-1)*4]*T[i]) + " days"

            #Lambert arc to reach Earth during (1-nu2)*T2 (second segment)
            dt = (1-x[8+(i-1)*4])*T[i]*DAY2SEC
            l = lambert_problem(r, r_P[i+1], dt, MU_SUN)
            v_end_l = l.get_v2()[0]
            v_beg_l = l.get_v1()[0]

            # Append data needed for potential plot generation
            plots_datas[i] = [v_out, x[8+(i-1)*4]*T[i]*DAY2SEC, l]

            #DSM occuring at time nu2*T2
            DV[i] = norm([a-b for a,b in zip(v_beg_l,v)])
            if P: print "DSM magnitude:    " + str(DV[i]) + "m/s"

        #Last Delta-v
        if P:  print "\nArrival at " + self.seq[-1].name
        DV[-1] = norm([a-b for a,b in zip(v_end_l,v_P[-1])])

        if P:
            print("Arrival Vinf:     " + str(DV[-1]) + "m/s  \n"
                  "Total mission time: " + str(sum(T)/365.25) + " years (" + str(sum(T)) + " days) \n"
                  "DSMs mag:         " + str(sum(DV[:-1])) + "m/s  \n"
                  "Entry Vel:        " + str(sqrt(2*(0.5*(DV[-1])**2 + 61933310.95))) + "m/s  \n"
                  "Entry epoch:      " + str(epoch(t_P[0].mjd2000 + sum(T))) )

        if doplot:
            import matplotlib as mpl
            from mpl_toolkits.mplot3d import Axes3D
            import matplotlib.pyplot as plt
            from PyKEP.orbit_plots import plot_planet, plot_lambert, plot_kepler

            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure(figsize=(8, 8), dpi=80)
            ax = fig.gca(projection='3d')
            ax.scatter(0,0,0, color='y')

            for i,planet in enumerate(self.seq):
                plot_planet(ax, planet, t0=t_P[i], color=(0.8,0.6,0.8), legend=True, units = AU)

            for i, tpl in enumerate(plots_datas):
                plot_kepler(ax, r_P[i], tpl[0], tpl[1], MU_SUN ,N = 100, color='b', legend=False, units = AU)
                plot_lambert(ax, tpl[2], sol = 0, color='r', legend=False, units = AU)

            fig.tight_layout()
            plt.show()

        if dists_class is not None:
            for i, tpl in enumerate(plots_datas):
                dists_class.positions_kepler(r_P[i], tpl[0], tpl[1], MU_SUN ,index = 2*i)
                dists_class.positions_lambert(tpl[2], sol = 0, index = 2*(i+.5))
            dists_class.set_launch_epoch(t_P[0].mjd2000)
            return dists_class

        if rtrn_desc:
            import numpy as np
            from math import atan2

            desc = dict()
            for i in range(len(x)):
                desc[i] = x[i]

            theta = 2*pi*x[1]
            phi = acos(2*x[2]-1)-pi/2

            axtl = 23.43929*DEG2RAD # Earth axial tlit
            mactran = np.matrix([ [1,         0,          0],
                                  [0, cos(axtl), -sin(axtl)],
                                  [0, sin(axtl),  cos(axtl)] ]) # Macierz przejscia z helio do ECI

            macdegs = np.matrix([ [cos(phi)*cos(theta)],
                                  [cos(phi)*sin(theta)],
                                  [sin(phi)] ])
            aaa = mactran.dot(macdegs)
            theta_earth = atan2(aaa[1,0], aaa[0,0])    # Rektascensja asym
            phi_earth = asin(aaa[2,0])                 # Deklinacja asym
            desc['RA'] = 360+RAD2DEG*theta_earth
            desc['DEC'] = phi_earth*RAD2DEG

            desc['Ldate'] = str(t_P[0])
            desc['C3'] = (x[3] / 1000)**2

            desc['dVmag'] = sum(DV[:-1])
            desc['VinfRE'] = DV[-1]
            desc['VRE'] = (DV[-1]**2 + 2*61933310.95)**.5
            desc['REDate'] = str(epoch(t_P[0].mjd2000 + sum(T)))

            for i in range(1,self.__n):
                desc[self.seq[i].name[0].capitalize()+'Dist'] = PLDists[i-1]
                desc[self.seq[i].name[0].capitalize()+'FBDate'] = str(FBDates[i-1])

            return desc

    def set_tof(self, minimum, maximum):
        """
        Sets the minimum and maximum time of flight allowed (in days)
        """
        lb = list(self.lb)
        ub = list(self.ub)
        lb[5] = minimum*365.25
        ub[5] = maximum*365.25
        self.set_bounds(lb,ub)

    def set_launch_window(self, start, end):
        """
        Sets the launch window allowed in terms of starting and ending epoch
        """
        lb = list(self.lb)
        ub = list(self.ub)
        lb[0] = start.mjd2000
        ub[0] = end.mjd2000
        self.set_bounds(lb,ub)

    def set_vinf(self, vinf):
        """
        Sets the allowed launch vinf (in km/s)
        """
        lb = list(self.lb)
        ub = list(self.ub)
        lb[3] = 0
        ub[3] = vinf * 1000
        self.set_bounds(lb,ub)
        self.C3_barrier = (vinf*1000)**2

    def set_entry_barrier(self, v):
        """
        Sets the earth reentry velocity barrier (in m/s)
        """
        self.entry_vel_barrier = v

