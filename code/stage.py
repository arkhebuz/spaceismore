# -*- coding: utf-8 -*-
from PyKEP import epoch, epoch_from_string
from consts import *


class ext_LUS(object):
    """Constructs and reprsents an ext-LUS stage (its parameters).
    """
    def __init__(self, engines, stage_additions = 702):
        """
        Ext-LUS stage constructor: ext_LUS(engines, stage_additions = 702)

        engines: 'rl' for RL10
                 'mb' for MB60
                 'j' for J2X

        stage_additions: mass (in kg)
        """
        # --------------------------------------------- LUS stage data [START]
        LH2_prop_base = 34334*lbm
        LO2_prop_base = 200853*lbm
        prop_total_base = LH2_prop_base + LO2_prop_base
        RCS_prop = 1432*lbm
        RL_MB_dry_mass_base = 26133*lbm
        RL_MB_total_mass_base = 262752*lbm
        usable_TLI_prop_base = 105000
        RL_MB_ascent_prop = 63000
        residual_prop_base = prop_total_base -usable_TLI_prop_base
        LEO_orbit = 130*nm
        prop_load_max = 125000

        LEO_payload_RL10 = 93100
        Isp_RL10 = 462.5
        thrust_RL10 = 99000*lbm

        LEO_payload_MB60 = 97000
        Isp_MB60 = 465
        thrust_MB60 = 120000*lbm

        J2X_ascent_prop = 95000
        Isp_J2X = 448
        LEO_payload_J2X = 105200
        thrust_J2X = 294000*lbm
        J2X_dry_mass_base = 27593*lbm
        # --------------------------------------------- LUS stage data [END]

        # --------------------------------------------- EUS-rl stage data [START]
        EUS_usable_prop = 276000*lbm
        EUS_prop_wth_unusables = 283000*lbm
        EUS_pseudodry_mass = 32000*lbm
        EUS_total_mass = 308000*lbm
        EUS_Isp_RL10 = 462
        EUS_thrust_RL10 = 4*24750*lbm
        EUS_LEO_payload_RL10 = 97500
        EUS_LEO_orbit = 130*nm
        EUS_RCS_prop = 1400*lbm
        # --------------------------------------------- EUS-rl stage data [END]


        if engines == 'rl':
            self.Isp = Isp_RL10
            self.Thrust = thrust_RL10
            LEO_payload = LEO_payload_RL10
            ascent_prop = RL_MB_ascent_prop
            dry_mass_base = RL_MB_dry_mass_base
            self.Engines = "4xRL-10"
            mass_const = 3000+0

        elif engines == 'mb':
            self.Isp = Isp_MB60
            self.Thrust = thrust_MB60
            LEO_payload = LEO_payload_MB60
            ascent_prop = RL_MB_ascent_prop
            dry_mass_base = RL_MB_dry_mass_base -2*200
            self.Engines = "2xMB-60"
            mass_const = 3000 -2*200

        elif engines == 'j':
            self.Isp = Isp_J2X
            self.Thrust = thrust_J2X
            LEO_payload = LEO_payload_J2X
            ascent_prop = J2X_ascent_prop
            dry_mass_base = J2X_dry_mass_base
            self.Engines = "1xJ2X"
            mass_const = 3000 +662

        elif engines == 'eus':
            self.Isp = EUS_Isp_RL10
            self.Thrust = EUS_thrust_RL10
            LEO_payload = EUS_LEO_payload_RL10
            ascent_prop = RL_MB_ascent_prop*EUS_total_mass/RL_MB_total_mass_base
            self.Engines = "4xRL-10"
            mass_const = 3000 +0
            RCS_prop = 0
            residual_prop_base = 0
            dry_mass_base = EUS_pseudodry_mass
            prop_load_max = EUS_prop_wth_unusables

        prop_load = prop_load_max
        brnout_mass_base = dry_mass_base + RCS_prop + residual_prop_base
        brt_mass_int_param = brnout_mass_base
        resid_mass_int_param = residual_prop_base

        for i in range(50):
            usable_LEO_prop = LEO_payload + brnout_mass_base -brt_mass_int_param    # kg
            prop_stretch = (usable_LEO_prop -(prop_load -ascent_prop -resid_mass_int_param))/prop_load    # %
            resid_mass_int_param = residual_prop_base*(1 + prop_stretch)    # kg
            brt_mass_int_param = dry_mass_base + (dry_mass_base -mass_const)*prop_stretch +stage_additions +resid_mass_int_param +RCS_prop    # kg
            #print resid_mass_int_param, brt_mass_int_param, usable_LEO_prop, prop_stretch

        self.Resid_mass = round(resid_mass_int_param, 2)
        self.Brt_mass = round(brt_mass_int_param, 2)
        self.Usable_prop = round(usable_LEO_prop, 2)
        self.Prop_stretch = round(prop_stretch, 5)
        self.Prop_flow = self.Thrust/self.Isp # kg/s
        self.Orbit_height = LEO_orbit
        self.Stage_additions = stage_additions
        #initialize some other variables
        self.Consumed_prop = 0
        self.Launch_date = 0
        self.Actual_payload = 0

    def vol(self):
        self.Resid_mass = (34334+200853)*lbm-105000
        self.Brt_mass = 26133*lbm + self.Resid_mass +1432*lbm
        self.Usable_prop = 90000 - self.Brt_mass
        self.Prop_stretch = 0.0
        self.Prop_flow = 99000*lbm/462.5
        self.Stage_additions = 702

    def set_boiloff(self, percentage):
        self.Boiloff = percentage

    def set_time(self, mjd2000_epoch):
        """Sets number of days since stage orbital insertion"""
        self.L_time = mjd2000_epoch - self.Launch_date
        self.Boiled_prop = self.L_time * self.Boiloff * self.Usable_prop

    def set_launch_date(self, mjd2000_epoch):
        self.Launch_date = mjd2000_epoch

    def set_payload(self, payload):
        """
        Stopień powinien wiedzieć o swoim ładunku, całkowitej masie, zuzytym już paliwie, paliwie ktore wyparowalo?
        """
        self.Actual_payload = payload

    def set_consumed_prop(self, prop_kg):
        self.Consumed_prop = abs(prop_kg)

    def compute_stack_mass(self):
        if (self.Usable_prop -self.Consumed_prop -self.Boiled_prop) < 0:    raise ValueError
        self.Stack_mass = self.Actual_payload + self.Brt_mass + self.Usable_prop -self.Consumed_prop -self.Boiled_prop

    def hit_from_incl(self, hit_ratio):
        self.Usable_prop = self.Usable_prop - self.Stack_mass*(1-hit_ratio)

    def print_info(self):
        print("Stage:               {0} based ext-LUS \n".format(self.Engines)
             +"LEO loiter time:     {0} days \n".format(self.L_time)
             +"Boiloff rate:        {0} %  \n".format(self.Boiloff*100)
             +"Stage mods:          {0} kg \n".format(self.Stage_additions)
             +"Stage burnout mass:  {0} kg \n".format(self.Brt_mass)
             +"Residual fuel:       {0} kg \n".format(self.Resid_mass)
             +"Usable prop:         {0} kg \n".format(self.Usable_prop)
             +"Tank stretch:        {0} %  \n".format(100*self.Prop_stretch) )


if __name__ == '__main__':
    offset = 16
    window = 4
    phasing = 4

    lus = ext_LUS('rl')
    lus.set_boiloff(0.002)
    lus.set_launch_date(epoch(3).mjd2000)
    lus.set_time(epoch(4).mjd2000)
    lus.print_info()
    print lus.Prop_flow, lus.Boiled_prop #*5/(24*lus.prop_flow)
    lus.set_payload(0)
    lus.compute_stack_mass()
    print lus.Usable_prop - lus.Stack_mass*0.1

    #print ext_LUS_brt_mass+ext_LUS_usable_prop
