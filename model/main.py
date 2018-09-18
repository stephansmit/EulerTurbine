from turbine_state import TurbineState
from radial_turbine import RadialTurbine
from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt
from math import pi
import numpy as np



#rotationalspeed
omega = 520
etats = 0.5
#heights
h1 = 8.897e-3
h2 = 8.897e-3
h3 =10.e-3

#radial postions
r1=.1680
r2=0.115
r3=0.0873

#inlet conditions
P1=31.95e5
T1=587.65
C1=0.5
alpha1 = 0

#inlet/rotor interface conditions
alpha2=85.2
DOR = 0.295658
zeta3 = 0
#outlet conditions
P3 = 0.08926012e5
init=480

massflow = .18
optimized = False
turbine = RadialTurbine(omega, h1, h2, h3, r1, r2, r3, P1, T1, alpha1, alpha2, zeta3, P3, massflow, optimized)

#print(turbine.get_turbine_info())
turbine.print_stages()
turbine.print_turbine()


#turbine.draw_velocity_triangle()
#print(turbine.get_turbine_info())
