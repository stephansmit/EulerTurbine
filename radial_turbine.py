import numpy as np
from turbine_state import TurbineState
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import collections


class RadialTurbine():
    def __init__(self, omega, h1, h2, h3, r1, r2, r3, P1, T1, C1, alpha1, alpha2, P3, DOR):

        self.states = []
        self.omega = omega
        self.state1 = TurbineState(0, r1, h1)
        self.state1.thermodynamic.set_staticPT(P1, T1)
        self.state1.thermodynamic.set_total_static_c(C1)
        self.state1.kinematic.set_state_alpha_cmag(alpha1, C1)
        self.state1.set_massflow()
        self.state1.kinematic.set_mach_numbers(self.state1.thermodynamic.static.A)


        self.state2 = TurbineState(omega, r2, h2)
        self.state2.set_previous(self.state1)
        self.state2.set_first(self.state1)
        self.state3 = TurbineState(omega, r3, h3)
        self.state3.set_previous(self.state2)
        self.state3.set_first(self.state1)

        self.state2.thermodynamic.set_total_with_etatotal_statetotal_ptotal(1.0,
                                                                            self.state1.thermodynamic.total,
                                                                            self.state1.thermodynamic.total.P)



        #determine total enthalpy
        # self.state3.thermodynamic.total.H = self.state2.thermodynamic.total.H -10000
        # self.state3.thermodynamic.set_htotal_with_etastatic_statetotal_pstatic(etats,
        #                                                                        self.state1.thermodynamic.total,
        #                                                                        P3)


        c0 = 470
        c2 = c0
        #set total and static conditions using assumption


        minimize(self.calc_error, c0, args=(alpha2, P3, DOR))





        self.states.append(self.state1)
        self.states.append(self.state2)
        self.states.append(self.state3)
        for i in self.states:
            i.kinematic.set_angles()
            i.thermodynamic.static.set_speedofsound()
            i.kinematic.set_mach_numbers(i.thermodynamic.static.A)
        #state2 = TurbineState(omega, r2, h2)
        #state2.set_previous(state1)
        self.calc_DOR_velocity()
        self.calc_DOR_velocity2()
        self.calc_DOR_enthalpy()
        self.set_work_enthalpy()
        self.set_work_velocity()

    def calc_error(self,c, alpha2, P3, DOR):
        #assume c2

        c2 = c[0]
        #set total and static conditions using assumption

        self.state2.kinematic.set_state_alpha_cmag(alpha2, c2)
        self.state2.thermodynamic.set_static_total_c(c2)
        self.state2.set_rothalpy()

        #set static properties with DOR and P3
        self.state3.set_static_dor_pstatic(DOR, P3)
        self.print_stages()

        #calculate radial absolute velocity using massflow
        c3r = self.state1.massflow/(self.state3.thermodynamic.static.D*self.state3.area)
        #using conservation of rothalpy calculate w3
        # print(self.state2.rothalpy - self.state3.thermodynamic.static.H + (self.state3.kinematic.u.mag**2/2))
        w3 = np.sqrt(2*(self.state2.rothalpy - self.state3.thermodynamic.static.H + (self.state3.kinematic.u.mag**2)/2))
        #make velocity triangle
        self.state3.kinematic.w.set_vector_with_rcomponent_mag(c3r, w3)
        self.state3.kinematic.set_c_with_w_u()

        #calculate total conditions using velocity triangle and static
        self.state3.set_htotal_work()
        # self.state3.thermodynamic.set_total_static_c(self.state3.kinematic.c.mag)

        self.state2.set_massflow()
        self.state3.set_massflow()
        self.state1.set_rothalpy()
        self.state3.set_rothalpy()
        return (self.state2.massflow - self.state1.massflow)**2

    def set_work_enthalpy(self):
        self.work_enthalpy = self.state2.thermodynamic.total.H -self.state3.thermodynamic.total.H
        return

    def set_work_velocity(self):
        self.work_velocity = self.state2.kinematic.u.theta*self.state2.kinematic.c.theta - \
            abs(self.state3.kinematic.u.theta*self.state3.kinematic.c.theta)
        return

    def calc_DOR_velocity(self):
        self.DOR_velocity = (self.state3.kinematic.w.theta ** 2 - self.state2.kinematic.w.theta ** 2
                    - self.state3.kinematic.u.theta ** 2 + self.state2.kinematic.u.theta ** 2
                    ) / \
                   (self.state2.kinematic.c.theta ** 2 - self.state1.kinematic.c.theta ** 2
                    + self.state3.kinematic.w.theta ** 2 - self.state2.kinematic.w.theta ** 2
                    - self.state3.kinematic.u.theta ** 2 + self.state2.kinematic.u.theta ** 2
                    )
        return

    def calc_DOR_velocity2(self):
        self.DOR_velocity2 = (self.state3.kinematic.w.mag**2 -
                             self.state2.kinematic.w.mag**2 -
                             self.state3.kinematic.u.mag**2 +
                             self.state2.kinematic.u.mag**2)/ \
                             (self.state2.kinematic.c.mag**2 -
                              self.state1.kinematic.c.mag**2 +
                              self.state3.kinematic.w.mag ** 2 -
                              self.state2.kinematic.w.mag ** 2 -
                              self.state3.kinematic.u.mag ** 2 +
                              self.state2.kinematic.u.mag ** 2

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
        print("DOR velocity: " + str(    self.DOR_velocity))
        print("DOR velocity2: " + str(    self.DOR_velocity2))
        print("DOR enthalpy: " + str(    self.DOR_enthalpy))
        print("Work velocity: " + str(self.work_velocity))
        print("Work enthalpy: " + str(self.work_enthalpy))

        print("--------------------------------------")

    def draw_velocity_triangle(self):
        originstatorinlet = [0, 0]
        originrotorinlet = [150, 0]
        originrotoroutlet = [300, 0]
        fig, ax = plt.subplots()
        self.state1.kinematic.draw_with_c(ax, originstatorinlet)
        self.state2.kinematic.draw_with_c(ax, originrotorinlet)
        self.state3.kinematic.draw_with_w(ax, originrotoroutlet)
        ax.set_xlim([0, 1000])
        ax.set_ylim([-500, 500])
        plt.show()




    def get_turbine_info(self):
        return self.flatten(dict({
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

column_names = [u'R', u'state1_area', u'state1_height', u'state1_kinematic_alpha',
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

