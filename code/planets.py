# -*- coding: utf-8 -*-
from PyKEP import planet_ss, epoch_from_string, epoch, planet, DEG2RAD, RAD2DEG
from datetime import datetime, timedelta
from numpy import matrix, cos, sin, pi
from scipy.optimize import fmin as simplex
from jplephem import Ephemeris
from sgp4.ext import rv2coe
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import de421


class de421_planets(object):
    """Helper class, for building PyKEP's planet objects representing
    solar system planets with positions being taken from JPL DE421 at
    specifed epoch.
    """
    def __init__(self):
        """initialize constants
        """
        axtl = 23.43929*DEG2RAD # Earth axial tlit
        # Macierz przejscia z  ECI do helio
        self.mactran = matrix([ [1, 0, 0],
                                [0, cos(axtl), sin(axtl)],
                                [0, -sin(axtl), cos(axtl)] ])

        self.eph = Ephemeris(de421)
        self.day = 86400.0

    def eq2eclipt(self, xyz):
        macxyz = matrix(xyz)
        return self.mactran.dot(macxyz)

    def errors_plot(self, planet, epoch_center_str, quality = False):
        ep_center = epoch_from_string(epoch_center_str)
        position, __ = self.return_pos_vel(planet, ep_center)

        F = lambda A, B: ( (A[0]-B[0])**2 + (A[1]-B[1])**2 + (A[2]-B[2])**2 )**0.5
        range_tuple = range(-30, 1)
        pos_errors = []
        weighted_pos_errors = []
        for i in range_tuple:
            time = ep_center.jd+i
            posKEP, __ = planet.eph(epoch(time, epoch.epoch_type.JD))
            pos_errors.append(F(position[30+i], posKEP)*10**-3)
            weighted_pos_errors.append(self.weights[30+i]*F(position[30+i], posKEP)*10**-3)

        if not quality:
            #plt.plot(range_tuple, weighted_pos_errors, 'o')
            plt.plot(range_tuple, pos_errors, 'o')
            plt.tight_layout()
        else:
            range_tuple = [datetime(2000, 01, 01) +timedelta(ep_center.mjd2000+d) for d in range(-30,1)]
            plt.figure(figsize=(7.5, 4.5), dpi=180)
            plt.plot(range_tuple, pos_errors, 'o', label=planet.name.capitalize())
            plt.ylabel('Position error magnitude [$km$]')
            plt.legend(loc='upper right', prop={'size':12})

            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
            plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator(minticks=6, maxticks=10))
            plt.gca().xaxis.set_minor_locator(mdates.DayLocator())
            plt.gcf().autofmt_xdate()

            plt.gca().yaxis.set_major_locator(plt.MultipleLocator(50.0))
            plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(25))

            plt.gca().xaxis.grid(True,'minor')
            #plt.gca().yaxis.grid(True,'minor')
            plt.gca().xaxis.grid(True,'major', linewidth=1)
            plt.gca().yaxis.grid(True,'major', linewidth=1)

            plt.tight_layout()
            plt.savefig('planet_position_errors_plot.png', dpi=200)
        plt.show()

    def compute_weights(self):
        """Compute the weights by date in the default 30-day fitting period
        """
        weights = [sin(pi/( 2+(0.04*day)**4 ))**8 for day in range(-30,1)]
        #  Last days are more important
        weights[30] = weights[30]+3
        weights[29] = weights[29]+2
        weights[28] = weights[28]+1
        self.weights = weights

    def return_pos_vel(self, fitted_planet, fit_end):
        if fitted_planet.name == 'earth':
            # there's no moon, just barycenter
            plstring = 'earthmoon'
        else:
            plstring = fitted_planet.name

        positions = []
        velocities = []
        for i in range(-30, 1):
            time = fit_end.jd+i
            posSun, velSun = self.eph.position_and_velocity('sun', time)
            pre_position, pre_velocity = self.eph.position_and_velocity(plstring, time)
            pre_position = self.eq2eclipt(pre_position - posSun)
            pre_velocity = self.eq2eclipt(pre_velocity - velSun)

            position = (1000*pre_position).reshape((1,3)).tolist()[0]
            velocity = (pre_velocity/self.day).reshape((1,3)).tolist()[0]
            velocities.append(velocity)
            positions.append(position)
        return positions, velocities

    def fit(self, fitted_planet, fit_end_str, pl_safe_distance):
        # Planet properities
        name = fitted_planet.name
        muself = fitted_planet.mu_self
        musun = fitted_planet.mu_central_body
        radius = fitted_planet.radius
        print "\nFitting keplerian approximation to:", [name]
        # Set cut-off window date
        fit_end = epoch_from_string(fit_end_str)
        # Initialize self.weights variable
        self.compute_weights()

        # Calculate planet's positions and velocities over defaut 30-day timespan
        positions, velocities = self.return_pos_vel(fitted_planet, fit_end)
        # Make initial guess (x0) for simplex algorithm
        pos_km = [ positions[-1][0]/1000,
                   positions[-1][1]/1000,
                   positions[-1][2]/1000 ]
        __, a, ecc, incl, omega, argp, __, m, __, __, __ = rv2coe(pos_km, velocities[-1], musun/(1000**3))
        x0 = (a*1000, ecc, incl, omega, argp, m)

        def f(x, fit_end, weights, poss):
            # Function to minimalize, returns sum of squares of position errors
            # over 30 day period, modified with weights.
            body = planet(fit_end, x.tolist(), musun, muself, radius, radius+pl_safe_distance, name)
            squares = []
            for i in range(-30, 1):
                r, __ = body.eph(epoch(fit_end.mjd2000+i, epoch.epoch_type.MJD2000))
                rorg = poss[30+i]
                w = weights[30+i]
                squares.append( w*sum([(r[n] - rorg[n])**2 for n in [0,1,2]]) )
            return sum(squares)

        # Minimalization sometimes doesn't end with "Optimization terminated successfully."
        # communicate, but the results are generally good enough even then.
        fitted_elements = simplex(f, x0, args=(fit_end, self.weights, positions), full_output=0)
        body = planet(fit_end, fitted_elements.tolist(), musun, muself, radius, radius+pl_safe_distance, name)
        return body, fitted_elements.tolist()

    def eme_500day_2018(self, flyby_dist):
        """
        Method returning planet sequence for 2017/2018 window for mars
        free-return flyby mission. Planets are constructed using hard-coded orbital
        elements, being the results of optimization coded above.
        """
        elements_earth_launch = [149596388450.0783, 0.016678335161958293, 4.225442260171425e-05, 3.0483146025310095, 5.03066988767222, 0.0416636903599645]
        elements_mars = [227935414087.26007, 0.0933324133158771, 0.03225613287661293, 0.8640335870435378, 5.003255854143598, 6.049811097247842]
        elements_earth_entry = [149599006432.44122, 0.016714632681253273, 4.563992533184915e-05, 3.082582152662022, 4.998724154778759, 2.3744754415524087]

        launchepoch = '2018-01-06 00:00:00'
        marsflybyepoch = '2018-Aug-22 00:00:00'
        entryepoch = '2019-May-22 00:00:00'

        sequence = [ (launchepoch, planet_ss('earth'), elements_earth_launch, 200*1000),
                     (marsflybyepoch, planet_ss('mars'), elements_mars, flyby_dist),
                     (entryepoch, planet_ss('earth'), elements_earth_entry, 0) ]

        planets = []
        for pl in sequence:
            name = pl[1].name
            muself = pl[1].mu_self
            musun = pl[1].mu_central_body
            radius = pl[1].radius
            body = planet(epoch_from_string(pl[0]), pl[2], musun, muself, radius, radius+pl[3], name)
            planets.append(body)

        return planets

    def eme_700day_2018(self, flyby_dist):
        """
        Method returning planet sequence for 2018 window for mars
        free-return 700-day flyby mission. Planets are constructed using hard-coded orbital
        elements, being the results of optimization coded above.
        """
        elements_earth_launch = [149600147887.68948, 0.016678913612078988, 4.2914948114489875e-05, 3.0511667130642532, 5.029905567754205, 2.3790693213784753]
        elements_mars = [227935414087.26007, 0.0933324133158771, 0.03225613287661293, 0.8640335870435378, 5.003255854143598, 6.049811097247842]
        elements_earth_entry = [149597795100.11346, 0.016733410527094122, 4.642665726484085e-05, 3.095208080733212, 4.985102945894168, 2.3882414620556034]

        launchepoch = '2018-05-22 00:00:00'
        marsflybyepoch = '2018-Aug-22 00:00:00'
        entryepoch = '2020-May-22 00:00:00'

        sequence = [ (launchepoch, planet_ss('earth'), elements_earth_launch, 200*1000),
                     (marsflybyepoch, planet_ss('mars'), elements_mars, flyby_dist),
                     (entryepoch, planet_ss('earth'), elements_earth_entry, 0) ]

        planets = []
        for pl in sequence:
            name = pl[1].name
            muself = pl[1].mu_self
            musun = pl[1].mu_central_body
            radius = pl[1].radius
            body = planet(epoch_from_string(pl[0]), pl[2], musun, muself, radius, radius+pl[3], name)
            planets.append(body)

        return planets

    def evme_2021(self, venus_fb_dist, mars_fb_dist):
        elements_earth_launch = [149597363832.46185, 0.016681697815423065, 4.851093253043914e-05, 3.069957419844765, 5.011961120616479, 5.7709834002754645]
        elements_venus = [108208096620.31639, 0.006763298656078312, 0.05924349187589259, 1.3372372757968485, 0.9548235480680097, 2.2058859198013865]
        elements_mars = [227938437328.67865, 0.09348154468802035, 0.032252632078193955, 0.8637491543029134, 5.003176772196014, 1.1391254143939635]
        elements_earth_entry = [149598150947.04254, 0.01670128416625776, 5.2674057801829184e-05, 3.043625307657902, 5.040664796636818, 3.2138772414785812]

        launchepoch = '2021-Dec-05 00:00:00'
        venusflybyepoch = '2022-Apr-12 00:00:00'
        marsflybyepoch = '2022-Oct-24 00:00:00'
        entryepoch = '2023-Jul-10 00:00:00'

        sequence = [ (launchepoch, planet_ss('earth'), elements_earth_launch, 200*1000),
                     (venusflybyepoch, planet_ss('venus'), elements_venus, venus_fb_dist),
                     (marsflybyepoch, planet_ss('mars'), elements_mars, mars_fb_dist),
                     (entryepoch, planet_ss('earth'), elements_earth_entry, 0) ]

        planets = []
        for pl in sequence:
            name = pl[1].name
            muself = pl[1].mu_self
            musun = pl[1].mu_central_body
            radius = pl[1].radius
            body = planet(epoch_from_string(pl[0]), pl[2], musun, muself, radius, radius+pl[3], name)
            planets.append(body)

        return planets


if __name__ == '__main__':
    a = de421_planets()
    fb = a.evme_2021(400*1000, 100*1000)

