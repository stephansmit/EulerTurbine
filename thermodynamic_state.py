from CoolProp.CoolProp import PropsSI
import numpy as np


class ThermodynamicState():
    def __init__(self,fluid="Toluene"):
        self.fluid = fluid
        self.total = ThermodynamicProperties()
        self.static = ThermodynamicProperties()

    #static conditions
    def set_staticPT(self, P,T):
        self.static.set_statePT(P,T)
    def set_staticHS(self, H, S ):
        self.static.set_stateHS(H,S)
    def set_staticHP(self, H, P ):
        self.static.set_stateHP(H,P)

    #total conditions
    def set_totalPT(self, P, T):
        self.total.set_statePT(P, T)
    def set_totalHS(self, H, S ):
        self.total.set_stateHS(H,S)
    def set_totalHP(self, H, P ):
        self.total.set_stateHP(H,P)

    def get_c_static_total(self):
        return np.sqrt(self.total.H - self.static.H)*2

    def set_total_static_c(self, c):
        htotal = self.static.H + (c**2)/2
        self.set_totalHS(htotal, self.static.S)

    def set_static_total_c(self, c):
        hstatic = self.total.H - (c**2)/2
        self.set_staticHS(hstatic, self.total.S)

    def set_static_htotal_c_p(self,c,P):
        hstatic = self.total.H - (c ** 2) / 2
        self.set_staticHP(hstatic, P)
        self.set_totalHS(self.total.H, self.static.S)

    def set_total_with_etatotal_statetotal_ptotal(self, etatotal, statetotal_in, ptotal_out):
        htotal_s_out = PropsSI('H', 'S', statetotal_in.S,'P', ptotal_out, self.fluid)
        htotal_out = statetotal_in.H - etatotal*(statetotal_in.H - htotal_s_out)
        self.set_totalHP(htotal_out, ptotal_out )

    def set_htotal_with_etastatic_statetotal_pstatic(self, etastatic, statetotal_in, pstatic_out):
        hstatic_s_out = PropsSI('H', 'S', statetotal_in.S, 'P', pstatic_out, self.fluid)
        htotal_out = statetotal_in.H - etastatic * (statetotal_in.H - hstatic_s_out)
        self.total.H = htotal_out



class ThermodynamicProperties():
    def __init__(self, fluid="Toluene", H=0, S=0, T=0, P=0, D=0, A=0):
        self.fluid = fluid
        self.H = H
        self.S = S
        self.T = T
        self.P = P
        self.D = D
        self.A = A
    def set_stateHP(self, H, P):
        self.P = P
        self.H = H
        self.T = PropsSI('T', 'H', self.H, 'P', self.P, self.fluid)
        self.D = PropsSI('D', 'H', self.H, 'P', self.P, self.fluid)
        self.S = PropsSI('S', 'H', self.H, 'P', self.P, self.fluid)
        self.A = PropsSI('A', 'H', self.H, 'P', self.P, self.fluid)
    def set_statePT(self, P, T):
        self.P = P
        self.T = T
        self.H = PropsSI('H', 'T', self.T, 'P', self.P, self.fluid)
        self.D = PropsSI('D', 'T', self.T, 'P', self.P, self.fluid)
        self.S = PropsSI('S', 'T', self.T, 'P', self.P, self.fluid)
        self.A = PropsSI('A', 'T', self.T, 'P', self.P, self.fluid)
    def set_stateHS(self, H,S):
        self.H=H
        self.S=S
        self.T = PropsSI('T', 'H', self.H, 'S', self.S, self.fluid)
        self.P = PropsSI('P', 'H', self.H, 'S', self.S, self.fluid)
        self.D = PropsSI('D', 'H', self.H, 'S', self.S, self.fluid)
        self.A = PropsSI('A', 'H', self.H, 'S', self.S, self.fluid)
