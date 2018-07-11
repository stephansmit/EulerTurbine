from thermodynamic_state import ThermodynamicState
from kinematic_state import KinematicState
from math import pi
from scipy.optimize import minimize
import numpy as np

class TurbineState():
    def __init__(self, omega, r,h):
        self.thermodynamic = ThermodynamicState()
        self.kinematic = KinematicState(omega, r)
        self.area = 2*pi *r*h
        self.massflow = 0

    def set_massflow(self):
        self.massflow = self.area*self.thermodynamic.static.D*self.kinematic.c.r
    def set_previous(self, state_in):
        self.previous_state=state_in

    def set_first(self, state_first):
        self.first_state=state_first






    def set_state_with_etatotal_ptotal_alpha(self, etatt, ptotal_out, alpha):
        c0=300
        minimize(self.calc_error_etatotal, c0, args=(etatt, ptotal_out, alpha))
    def set_state_with_etastatic_pstatic_alpha(self, etats, pstatic_out, alpha):
        c0 = 300
        minimize(self.calc_error_etastatic, c0, args=(etats, pstatic_out, alpha))
    def set_state_with_etastatic_pstatic_R(self, etats, pstatic_out, R):
        alpha0 = 40
        minimize(self.calc_error_etastatic_alpha, alpha0, args=(etats, pstatic_out, R))



    #error functions for minization
    def calc_error_etastatic_alpha(self,alpha,etats, pstatic_out, R):
        print(alpha[0])
        h_static_out = (self.previous_state.thermodynamic.static.H - R*self.first_state.thermodynamic.static.H )/(1-R)
        self.thermodynamic.set_staticHP(h_static_out,pstatic_out)

        self.thermodynamic.set_htotal_with_etastatic_statetotal_pstatic(etats,
                                                                        self.previous_state.thermodynamic.total,
                                                                        pstatic_out)

        c_mag = self.thermodynamic.get_c_static_total()
        self.thermodynamic.set_total_static_c(c_mag)
        self.kinematic.set_state_alpha_cmag(alpha[0], c_mag)
        self.set_massflow()
        error = (self.massflow - self.previous_state.massflow)**2
        return error

    def calc_error_etastatic(self, c, etats, pstatic_out, alpha):

        self.thermodynamic.set_htotal_with_etastatic_statetotal_pstatic(etats,
                                                                        self.previous_state.thermodynamic.total,
                                                                        pstatic_out)
        self.thermodynamic.set_static_htotal_c_p(c[0], pstatic_out)
        self.kinematic.set_state_alpha_cmag(alpha, c[0])
        self.set_massflow()
        error = (self.massflow - self.previous_state.massflow)**2
        return error
    def calc_error_etatotal(self, c, etatt, ptotal_out, alpha):
        self.thermodynamic.set_total_with_etatotal_statetotal_ptotal(etatt,
                                                                     self.previous_state.thermodynamic.total,
                                                                     ptotal_out)
        self.kinematic.set_state_alpha_cmag(alpha, c[0])
        self.thermodynamic.set_static_total_c(c[0])
        self.set_massflow()
        error = (self.massflow - self.previous_state.massflow)**2
        return error

    def print_state(self):
        print("total Enthalpy: "     +str(self.thermodynamic.total.H))
        print("total Entropy: "      +str(self.thermodynamic.total.S))
        print("total Density: "      +str(self.thermodynamic.total.D))
        print("total Pressure: "     +str(self.thermodynamic.total.P))
        print("total Temperature: "  +str(self.thermodynamic.total.T))
        print("")
        print("static Enthalpy: "     +str(self.thermodynamic.static.H))
        print("static Entropy: "      +str(self.thermodynamic.static.S))
        print("static Density: "      +str(self.thermodynamic.static.D))
        print("static Pressure: "     +str(self.thermodynamic.static.P))
        print("static Temperature: "  +str(self.thermodynamic.static.T))
        print("")
        print("Mag Absolute Velocity: " + str(np.linalg.norm(self.kinematic.c.vec)))
        print("Absolute Velocity: " + str(self.kinematic.c.vec))
        print("Mag Relative Velocity: "  +str(np.linalg.norm(self.kinematic.w.vec)))
        print("Relative Velocity: "  +str(self.kinematic.w.vec))
        print("Absolute Mach: " + str(self.kinematic.c.mach))
        print("Relative Mach: " + str(self.kinematic.w.mach))
        print("Rotational Velocity: "+str(self.kinematic.u.vec))
        print("")
        print("Massflow: "+ str(self.massflow))
        print("Area: "+ str(self.area))

