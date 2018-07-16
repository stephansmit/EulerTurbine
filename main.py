from turbine_state import TurbineState
from radial_turbine import RadialTurbine
from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt
from math import pi
import numpy as np



#rotationalspeed
omega = 2600
etats = 0.5
#heights
h1 = 0.002
h2 = 0.008
h3 = 0.01

#radial postions
r1=.210
r2=0.115
r3=0.090

#inlet conditions
P1=32e5
T1=525
C1=0.5
alpha1 = 45

#inlet/rotor interface conditions
alpha2=77
DOR = 0.05
#outlet conditions
P3 = .8e5


turbine = RadialTurbine(omega, h1, h2, h3, r1, r2, r3, P1, T1, C1, alpha1, alpha2, P3, DOR)


# turbine.print_stages()
# turbine.print_turbine()
# turbine.draw_velocity_triangle()
print(turbine.get_turbine_info())