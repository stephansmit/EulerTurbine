import numpy as np
from turbine_state import TurbineState
from scipy.optimize import minimize
#import matplotlib.pyplot as plt
import collections
import math

class RadialTurbine():
    def __init__(self, hz, h1, h2, h3, r1, r2, r3, P01, T01, alpha1, alpha2, zeta3, P3, massflow, optimizedh):

        omega = hz*2*math.pi
        self.optimized = optimizedh
        self.desired_massflow = massflow
        self.states = []
        self.omega = omega
        self.state1 = TurbineState(0, r1, h1)
        self.state1.thermodynamic.set_totalPT(P01,T01)
        minimize(self.calc_static_given_massflow_alpha1, 0.001, args=(alpha1))

        self.state2 = TurbineState(omega, r2, h2)
        self.state2.set_previous(self.state1)
        self.state2.set_first(self.state1)

        self.state2.thermodynamic.set_total_with_etatotal_statetotal_ptotal(1.0,
                                                                            self.state1.thermodynamic.total,
                                                                            self.state1.thermodynamic.total.P)
        self.state2.massflow = self.state1.massflow

        minimize(self.calc_P2_given_massflow_alpha2, 0.01e5, args=(alpha2, self.state2.thermodynamic.total.S), method='Nelder-Mead' )

        self.state2.set_rothalpy()

        self.state3 = TurbineState(omega, r3, h3)
        self.state3.massflow = self.state1.massflow

        self.state2.set_next(self.state3)
        self.state3.set_previous(self.state2)
        self.state3.set_first(self.state1)

        self.state3.kinematic.zeta = zeta3
        self.state3.thermodynamic.set_staticPS(P3, self.state2.thermodynamic.static.S)
        w3 = np.sqrt(2 * (self.state2.rothalpy - self.state3.thermodynamic.static.H + (self.state3.kinematic.u.theta ** 2) / 2))
        #w3thetar = w3 * np.cos(np.radians(zeta3))

        if self.optimized:
            cmag = np.sqrt(w3**2-self.state3.kinematic.u.mag**2)
            self.state3.kinematic.set_state_alpha_cmag(0.0,cmag )
            h3= self.state3.massflow/(self.state3.kinematic.c.mag*self.state3.thermodynamic.static.D*self.state3.circumferential)
            self.state3.set_area_h(h3)
        else:
            c3r = self.state1.massflow / (self.state3.thermodynamic.static.D * self.state3.area)
            self.state3.kinematic.w.set_vector_with_rcomponent_mag(c3r, w3)
            self.state3.kinematic.set_c_with_w_u()

        self.state3.thermodynamic.total.H = self.state2.rothalpy + self.state3.kinematic.u.theta * self.state3.kinematic.c.theta
        self.state3.thermodynamic.set_totalHS(self.state3.thermodynamic.total.H, self.state3.thermodynamic.static.S)

        self.state3.set_massflow()
        self.state1.set_rothalpy()
        self.state3.set_rothalpy()

        self.states.append(self.state1)
        self.states.append(self.state2)
        self.states.append(self.state3)
        for i in self.states:
            i.kinematic.set_angles()
            i.thermodynamic.static.set_speedofsound()
            i.kinematic.set_mach_numbers(i.thermodynamic.static.A)

        self.calc_DOR_velocity2()
        self.calc_DOR_enthalpy()
        self.set_work_enthalpy()
        self.set_work_velocity()
        self.calc_hdiff()

    # def calc_error3(self,h3):
    #     h2 = h2[0]
    #     self.state2.set_area_h(h2)
    #     self.state2.set_massflow()
    #     return (self.state2.massflow - self.state1.massflow)**2
    #
    def calc_static_given_massflow_alpha1(self,c1, alpha1):
        c1 = c1[0]
        self.state1.kinematic.set_state_alpha_cmag(alpha1, c1)
        self.state1.thermodynamic.set_static_total_c(c1)
        self.state1.set_massflow()
        return (self.desired_massflow - self.state1.massflow)**2

    def calc_P2_given_massflow_alpha2(self,P2, alpha2, s2):
        p2 = P2[0]
        self.state2.thermodynamic.set_staticPS(p2,s2)
        cmag = self.state2.thermodynamic.get_c_static_total()
        self.state2.kinematic.set_state_alpha_cmag(alpha2, cmag)
        massflow = self.state2.area*self.state2.thermodynamic.static.D*self.state2.kinematic.c.r
        return (self.state2.massflow-massflow)**2

    def calc_H3_setting_alpha3tozero(self,H3):
        h3=H3[0]
        self.state3.set_area_h(h3)
        minimize(self.calc_P3_rothalpy, [.1e5], args=(self.state1.thermodynamic.static.S), method='Nelder-Mead')
        # c3r = self.state1.massflow/(self.state3.thermodynamic.static.D*self.state3.area)
        # self.state3.kinematic.w.set_vector_with_rcomponent_mag(c3r, w3thetar)
        # self.state3.kinematic.set_c_with_w_u()
        self.state3.kinematic.set_angles()
        return (self.state3.kinematic.alpha)**2

    def calc_P3_rothalpy(self,P3, s3):
        p3=P3[0]
        #calc static conditions
        self.state3.thermodynamic.set_staticPS(p3,s3)
        #calc mag. relative velocity
        w3 = np.sqrt(2*(self.state2.rothalpy - self.state3.thermodynamic.static.H + (self.state3.kinematic.u.theta**2)/2))
        #calc radial component
        c3r = self.state1.massflow / (self.state3.thermodynamic.static.D * self.state3.area)
        self.state3.kinematic.w.set_vector_with_rcomponent_mag(c3r, w3)
        self.state3.kinematic.set_c_with_w_u()
        h03 = self.state3.thermodynamic.static.H +0.5*self.state3.kinematic.c.mag**2
        #self.state3.thermodynamic.total.H = self.state2.thermodynamic.total.H - \
        #                                    (self.state2.kinematic.c.theta * self.state2.kinematic.u.theta-
        #                                    self.state3.kinematic.c.theta *self.state3.kinematic.u.theta)
        self.state3.thermodynamic.total.H = self.state2.rothalpy + self.state3.kinematic.c.theta*self.state3.kinematic.u.theta

        print(h03)
        print(self.state3.thermodynamic.total.H)
        return (self.state3.thermodynamic.total.H-h03)**2


    def calc_error2(self,h2):
        h2 = h2[0]
        self.state2.set_area_h(h2)
        self.state2.set_massflow()
        return (self.state2.massflow - self.state1.massflow)**2

    def calc_error(self,c, alpha2, P3, DOR):
        c2 = c[0]
        self.state2.kinematic.set_state_alpha_cmag(alpha2, c2)
        self.state2.thermodynamic.set_static_total_c(c2)
        self.state2.set_rothalpy()
        self.state2.set_massflow()
        return (self.state2.massflow - self.state1.massflow)**2

    def set_work_enthalpy(self):
        self.work_enthalpy = self.state2.thermodynamic.total.H - self.state3.thermodynamic.total.H
        return

    def set_work_velocity(self):
        # if self.state3.kinematic.c.theta > 0:
            self.work_velocity = self.state2.kinematic.u.theta*self.state2.kinematic.c.theta - \
                self.state3.kinematic.u.theta*self.state3.kinematic.c.theta
        # else:
        #     self.work_velocity = self.state2.kinematic.u.theta * self.state2.kinematic.c.theta + \
        #                          self.state3.kinematic.u.theta * self.state3.kinematic.c.theta
        #
        # return




    def calc_hdiff(self):
        self.hdiff = self.state2.kinematic.u.theta*self.state2.kinematic.omega - self.state3.kinematic.u.theta*self.state2.kinematic.omega
    def calc_DOR_velocity2(self):
        self.DOR_velocity2 = ( self.state3.kinematic.w.mag**2
                              -self.state2.kinematic.w.mag**2
                              -self.state3.kinematic.u.mag**2 +
                              +self.state2.kinematic.u.mag**2)/\
                             ( self.state2.kinematic.c.mag**2
                              -self.state1.kinematic.c.mag**2
                              +self.state3.kinematic.w.mag ** 2
                              -self.state2.kinematic.w.mag ** 2
                              -self.state3.kinematic.u.mag ** 2
                              +self.state2.kinematic.u.mag ** 2
                              )
        return
    def calc_DOR_enthalpy(self):
        self.DOR_enthalpy = (self.state2.thermodynamic.static.H
                             - self.state3.thermodynamic.static.H)/ \
                            (self.state1.thermodynamic.static.H
                             -self.state3.thermodynamic.static.H)
        return
        #state2.set_state_with_etastatic_pstatic_alpha(etaTS12, P2, alpha2)
        #state2.kinematic.set_mach_numbers(state2.thermodynamic.static.A)


    def print_stages(self):
        for index, state in enumerate(self.states):
            print("--------------------------------------")
            print("State "+ str(index))
            print("--------------------------------------")
            state.print_state()

    def print_turbine(self):
        print("--------------------------------------")
        print("DOR velocity2: " + str(    self.DOR_velocity2))
        print("DOR enthalpy: " + str(    self.DOR_enthalpy))
        print("HDIFF: " + str(self.hdiff))
        print("Rotational speed: " + str(self.omega))
        print("Work enthalpy: " + str(self.work_enthalpy))
        print("Work velocity: " + str(self.work_velocity))
        print("test: " +str(self.state3.thermodynamic.static.H+0.5*self.state3.kinematic.c.mag**2))

        print("--------------------------------------")

#    def draw_velocity_triangle(self):
#        originstatorinlet = [0, 0]
#        originrotorinlet = [150, 0]
#        originrotoroutlet = [300, 0]
#        fig, ax = plt.subplots()
#        self.state1.kinematic.draw_with_c(ax, originstatorinlet)
#        self.state2.kinematic.draw_with_c(ax, originrotorinlet)
#        self.state3.kinematic.draw_with_w(ax, originrotoroutlet)
#        ax.set_xlim([0, 1000])
#        ax.set_ylim([-500, 500])
#        plt.show()
#



    def get_turbine_info(self):
        return self.flatten(dict({
               "optimized": self.optimized,
               "R": self.DOR_enthalpy,
               "omega": self.omega,
               "work": self.work_enthalpy,
               "state1": self.state1.get_turbinestate_info(),
               "state2": self.state2.get_turbinestate_info(),
               "state3": self.state3.get_turbinestate_info()}))


    def flatten(self,d, parent_key='', sep='_'):
        items = []
        for k, v in d.items():
            new_key = parent_key + sep + k if parent_key else k
            if isinstance(v, collections.MutableMapping):
                items.extend(self.flatten(v, new_key, sep=sep).items())
            else:
                items.append((new_key, v))
        return dict(items)

column_names = [u'R', u'omega', u'state1_area', u'state1_height', u'state1_kinematic_alpha',
       u'state1_kinematic_beta', u'state1_kinematic_c_mag',
       u'state1_kinematic_c_r', u'state1_kinematic_c_theta',
       u'state1_kinematic_u_mag', u'state1_kinematic_u_r',
       u'state1_kinematic_u_theta', u'state1_kinematic_w_mag',
       u'state1_kinematic_w_r', u'state1_kinematic_w_theta',
       u'state1_massflow', u'state1_radius', u'state1_rothalpy',
       u'state1_thermodynamic_static_D', u'state1_thermodynamic_static_H',
       u'state1_thermodynamic_static_P', u'state1_thermodynamic_static_S',
       u'state1_thermodynamic_static_T', u'state1_thermodynamic_total_D',
       u'state1_thermodynamic_static_A',
       u'state1_thermodynamic_total_H', u'state1_thermodynamic_total_P',
       u'state1_thermodynamic_total_S', u'state1_thermodynamic_total_T',
       u'state1_thermodynamic_total_A',   
       u'state2_area', u'state2_height', u'state2_kinematic_alpha',
       u'state2_kinematic_beta', u'state2_kinematic_c_mag',
       u'state2_kinematic_c_r', u'state2_kinematic_c_theta',
       u'state2_kinematic_u_mag', u'state2_kinematic_u_r',
       u'state2_kinematic_u_theta', u'state2_kinematic_w_mag',
       u'state2_kinematic_w_r', u'state2_kinematic_w_theta',
       u'state2_massflow', u'state2_radius', u'state2_rothalpy',
       u'state2_thermodynamic_static_D', u'state2_thermodynamic_static_H',
       u'state2_thermodynamic_static_P', u'state2_thermodynamic_static_S',
       u'state2_thermodynamic_static_T', u'state2_thermodynamic_total_D',
       u'state2_thermodynamic_static_A',
       u'state2_thermodynamic_total_H', u'state2_thermodynamic_total_P',
       u'state2_thermodynamic_total_S', u'state2_thermodynamic_total_T',
       u'state2_thermodynamic_total_A',
       u'state3_area', u'state3_height', u'state3_kinematic_alpha',
       u'state3_kinematic_beta', u'state3_kinematic_c_mag',
       u'state3_kinematic_c_r', u'state3_kinematic_c_theta',
       u'state3_kinematic_u_mag', u'state3_kinematic_u_r',
       u'state3_kinematic_u_theta', u'state3_kinematic_w_mag',
       u'state3_kinematic_w_r', u'state3_kinematic_w_theta',
       u'state3_massflow', u'state3_radius', u'state3_rothalpy',
       u'state3_thermodynamic_static_D', u'state3_thermodynamic_static_H',
       u'state3_thermodynamic_static_P', u'state3_thermodynamic_static_S',
       u'state3_thermodynamic_static_T', u'state3_thermodynamic_total_D',
       u'state3_thermodynamic_static_A',
       u'state3_thermodynamic_total_H', u'state3_thermodynamic_total_P',
       u'state3_thermodynamic_total_S', u'state3_thermodynamic_total_T',
       u'state3_thermodynamic_total_A',
       u'work']

empty_dict = dict()
for i in column_names:
    empty_dict[i]= np.nan

