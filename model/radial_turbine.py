import numpy as np
from turbine_state import TurbineState
from scipy.optimize import minimize
#import matplotlib.pyplot as plt
import collections
import math

class RadialTurbine():
    def __init__(self, hz, h1, h2, h3, r1, r2, r3, P01, T01, alpha1, P2, zeta3, P3, P4,eta1, eta2,massflow, optimizedh):

        omega = hz*2*math.pi

        self.optimized = optimizedh
        self.desired_massflow = massflow
        self.states = []
        self.omega = omega
        self.eta1 = eta1
        self.eta2 = eta2

        #-----------------------------------------
        #      turbine state 1 (inlet stator)
        # -----------------------------------------

        self.state1 = TurbineState(0, r1, h1)
        self.state1.thermodynamic.set_totalPT(P01,T01)
        self.state1.thermodynamic.isentropic_total = self.state1.thermodynamic.total

        self.state1.set_rothalpy()
        minimize(self.calc_cmag_given_massflow_alpha1, 0.001, args=(alpha1))
        self.state1.thermodynamic.isentropic_static = self.state1.thermodynamic.static

        #-----------------------------------------
        #      turbine state 2 (stator/rotor interface)
        # -----------------------------------------

        self.state2 = TurbineState(omega, r2, h2)
        self.state2.massflow = self.state1.massflow
        self.state2.set_previous(self.state1)
        self.state2.set_first(self.state1)
        #self.state2.thermodynamic.total = self.state1.thermodynamic.isentropic_total

        self.state2.thermodynamic.isentropic_total = self.state1.thermodynamic.isentropic_total
        self.state2.thermodynamic.set_isentropic_staticPS(P2,self.state1.thermodynamic.total.S)
        self.state2.thermodynamic.total.H = self.state1.thermodynamic.total.H
        c2s = np.sqrt(2*(self.state2.thermodynamic.total.H-self.state2.thermodynamic.isentropic_static.H))
        c2mag = eta1*c2s
        minimize(self.calc_S2_given_c2P2, self.state1.thermodynamic.total.S, args=(P2, c2mag), method='Nelder-Mead' )
        self.state2.thermodynamic.set_totalHS(self.state2.thermodynamic.total.H, self.state2.thermodynamic.static.S)
        c2r = self.state2.massflow/self.state2.area*self.state2.thermodynamic.static.D
        self.state2.kinematic.set_state_with_cr_cmag(c2r,c2mag)
        self.state2.set_rothalpy()

        #-----------------------------------------
        #      turbine state 3 (rotor outlet)
        # -----------------------------------------

        self.state3 = TurbineState(omega, r3, h3)
        self.state3.massflow = self.state1.massflow
        self.state3.set_previous(self.state2)
        self.state3.set_first(self.state1)
        # self.state3.kinematic.zeta = zeta3
        self.state3.thermodynamic.set_isentropic_staticPS(P3, self.state2.thermodynamic.static.S)

        w3s = np.sqrt(2 * (self.state2.rothalpy - self.state3.thermodynamic.isentropic_static.H + (self.state3.kinematic.u.theta ** 2) / 2))
        w3 =eta2 *w3s
        print(w3s)

        minimize(self.calc_S3_given_w3P3, self.state3.thermodynamic.isentropic_static.S, args=(P3,w3),method='Nelder-Mead')
        # #w3thetar = w3 * np.cos(np.radians(zeta3))
        if self.optimized:
            #pass
            cmag = np.sqrt(w3**2-self.state3.kinematic.u.mag**2)
            self.state3.kinematic.set_state_alpha_cmag(0.0,cmag )
            h3= self.state3.massflow/(self.state3.kinematic.c.mag*self.state3.thermodynamic.static.D*self.state3.circumferential)
            self.state3.set_area_h(h3)
        else:
            c3r = self.state1.massflow / (self.state3.thermodynamic.static.D * self.state3.area)
            self.state3.kinematic.w.set_vector_with_rcomponent_mag2(c3r, w3)
            self.state3.kinematic.set_c_with_w_u()
        self.state3.thermodynamic.total.H = self.state2.rothalpy + self.state3.kinematic.u.theta * self.state3.kinematic.c.theta
        self.state3.thermodynamic.set_totalHS(self.state3.thermodynamic.total.H, self.state3.thermodynamic.static.S)
        self.state3.thermodynamic.set_isentropic_totalHS(self.state3.thermodynamic.total.H, self.state3.thermodynamic.isentropic_static.S)
        self.state3.set_massflow()
        self.state3.set_rothalpy()


        #-----------------------------------------
        #      turbine state 4 (condensor inlet)
        # -----------------------------------------

        self.state4 = TurbineState(omega, r3, h3)
        self.state4.thermodynamic.static.P = P4
        self.state4.area = 10000
        c4mag = 0.01

        self.state4.thermodynamic.total.H = self.state3.thermodynamic.total.H
        self.state4.thermodynamic.static.H = self.state4.thermodynamic.total.H - .5*c4mag**2
        self.state4.thermodynamic.set_staticHP(self.state4.thermodynamic.static.H, self.state4.thermodynamic.static.P)
        self.state4.thermodynamic.set_totalHS(self.state4.thermodynamic.total.H, self.state4.thermodynamic.static.S)
        self.state4.thermodynamic.set_isentropic_staticHS(self.state4.thermodynamic.total.H, self.state4.thermodynamic.static.S)
        self.state4.thermodynamic.set_isentropic_totalHS(self.state4.thermodynamic.total.H, self.state4.thermodynamic.static.S)


        # -----------------------------------------
        #     All States
        # -----------------------------------------

        self.states.append(self.state1)
        self.states.append(self.state2)
        self.states.append(self.state3)
        self.states.append(self.state4)

        for i in self.states:
            i.kinematic.set_angles()
            i.thermodynamic.static.set_speedofsound()
            i.thermodynamic.isentropic_static.set_speedofsound()
            i.thermodynamic.isentropic_total.set_speedofsound()
            i.thermodynamic.isentropic_total.set_speedofsound()
            i.kinematic.set_mach_numbers(i.thermodynamic.static.A)

        self.calc_DOR_velocity()
        self.calc_DOR_enthalpy()
        self.set_work_enthalpy()
        self.set_work_velocity()


    def calc_cmag_given_massflow_alpha1(self,c1, alpha1):
        c1 = c1[0]
        self.state1.kinematic.set_state_alpha_cmag(alpha1, c1)
        self.state1.thermodynamic.set_static_total_c(c1)
        self.state1.set_massflow()
        return (self.desired_massflow - self.state1.massflow)**2

    # def calc_P2_given_massflow_alpha2_2(self,P2, alpha2, s2):
    #     p2 = P2[0]
    #     self.state2.thermodynamic.set_staticPS(p2,s2)
    #     cmag = self.state2.thermodynamic.get_c_static_total()
    #     self.state2.kinematic.set_state_alpha_cmag(alpha2, cmag)
    #     massflow = self.state2.area*self.state2.thermodynamic.static.D*self.state2.kinematic.c.r
    #     return (self.state2.massflow-massflow)**2
    #
    def calc_S3_given_w3P3(self,S3, P3, w3):
        s3 = S3[0]
        self.state3.thermodynamic.set_staticPS( P3,s3)
        w3n = np.sqrt(2 * (self.state2.rothalpy - self.state3.thermodynamic.static.H + (self.state3.kinematic.u.theta ** 2) / 2))
        return (w3-w3n)**2

    def calc_S2_given_c2P2(self,S2, P2, c2):
        s2 = S2[0]
        self.state2.thermodynamic.set_staticPS( P2,s2)
        ho2 = self.state2.thermodynamic.static.H + 0.5*c2**2

        return (self.state2.thermodynamic.isentropic_total.H-ho2)**2

    def calc_alpha2_given_massflow_P2(self,alpha2, P2, S2):
        a2 = alpha2[0]
        self.state2.thermodynamic.set_staticPS(P2,S2)
        cmag = self.state2.thermodynamic.get_c_static_total()
        self.state2.kinematic.set_state_alpha_cmag(a2, cmag)
        massflow = self.state2.area*self.state2.thermodynamic.static.D*self.state2.kinematic.c.r
        return (self.state2.massflow-massflow)**2

    def calc_P2_given_massflow_alpha2(self,P2, alpha2, s2):
        p2 = P2[0]
        self.state2.thermodynamic.set_staticPS(p2,s2)
        cmag = self.state2.thermodynamic.get_c_static_total()
        self.state2.kinematic.set_state_alpha_cmag(alpha2, cmag)
        massflow = self.state2.area*self.state2.thermodynamic.static.D*self.state2.kinematic.c.r
        return (self.state2.massflow-massflow)**2

    def set_work_enthalpy(self):
        self.work_enthalpy = self.state2.thermodynamic.total.H - self.state3.thermodynamic.total.H
        return

    def set_work_velocity(self):
        self.work_velocity = self.state2.kinematic.u.theta*self.state2.kinematic.c.theta - \
                self.state3.kinematic.u.theta*self.state3.kinematic.c.theta

    def calc_DOR_velocity(self):
        self.DOR_velocity = ( self.state3.kinematic.w.mag**2
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


    def print_stages(self):
        for index, state in enumerate(self.states):
            print("--------------------------------------")
            print("State "+ str(index))
            print("--------------------------------------")
            state.print_state()

    def print_turbine(self):
        print("--------------------------------------")
        print("DOR velocity2: " + str(    self.DOR_velocity))
        print("DOR enthalpy: " + str(    self.DOR_enthalpy))
        print("Rotational speed: " + str(self.omega))
        print("Work enthalpy: " + str(self.work_enthalpy))
        print("Work velocity: " + str(self.work_velocity))
        print("Should be H03: " +str(self.state3.thermodynamic.static.H+0.5*self.state3.kinematic.c.mag**2))
        print("--------------------------------------")



    def get_turbine_info(self):
        return self.flatten(dict({
               "optimized": self.optimized,
               "R": self.DOR_enthalpy,
               "omega": self.omega,
                "eta1": self.eta1,
                "eta2": self.eta2,
            "work": self.work_enthalpy,
               "state1": self.state1.get_turbinestate_info(),
               "state2": self.state2.get_turbinestate_info(),
               "state3": self.state3.get_turbinestate_info(),
               "state4": self.state4.get_turbinestate_info()}))


    def flatten(self,d, parent_key='', sep='_'):
        items = []
        for k, v in d.items():
            new_key = parent_key + sep + k if parent_key else k
            if isinstance(v, collections.MutableMapping):
                items.extend(self.flatten(v, new_key, sep=sep).items())
            else:
                items.append((new_key, v))
        return dict(items)


