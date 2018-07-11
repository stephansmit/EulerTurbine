from turbine_state import TurbineState
import matplotlib.pyplot as plt
from math import pi

P1=32e5
T1=525
C1=0.5
omega1=0
r1=.210
alpha1 = 45
height1 = 0.005
#P02=31.5e5
#etaTT12=0.99


etaTS12=0.15
P2 =0.8e5
alpha2=80
omega2 = 2600
r2=0.115
height2 = 0.009


etaTS23=0.65
r3=0.090
P3 = 0.2e5
alpha3 = 50


state1 = TurbineState( omega1, r1, height1)
state1.thermodynamic.set_staticPT(P1,T1)
state1.thermodynamic.set_total_static_c(C1)
state1.kinematic.set_state_alpha_cmag(alpha1, C1)
state1.set_massflow()
state1.kinematic.set_mach_numbers(state1.thermodynamic.static.A)


state2 = TurbineState( omega2, r2, height2)
state2.set_previous(state1)
state2.set_state_with_etastatic_pstatic_alpha(etaTS12, P2, alpha2)
state2.kinematic.set_mach_numbers(state2.thermodynamic.static.A)

# state2.set_state_with_etatotal_ptotal_alpha(etaTT12, P02, alpha2)

height3 = 0.05
state3 = TurbineState(omega2, r3,height3)
state3.set_previous(state2)
state3.set_first(state1)
R=0.0
state3.set_state_with_etastatic_pstatic_R(etaTS23, P3, R)
# state3.set_state_with_etastatic_pstatic_alpha(etaTS23, P3, alpha3)
state3.kinematic.set_mach_numbers(state3.thermodynamic.static.A)

DOR = (state2.thermodynamic.static.H - state3.thermodynamic.static.H)/ \
      (state1.thermodynamic.static.H- state3.thermodynamic.static.H)





print("")
print("STATE 1")
print("")
state1.print_state()

print("")
print("STATE 2")
print("")
state2.print_state()

print("")
print("STATE 3")
print("")
state3.print_state()
print("DOR: "+ str(DOR))
#
originstatorinlet = [0, 0]
originrotorinlet = [150, 0]
originrotoroutlet = [300, 0]
fig, ax = plt.subplots()
state1.kinematic.draw_with_c(ax, originstatorinlet)
state2.kinematic.draw_with_c(ax, originrotorinlet)
state3.kinematic.draw_with_w(ax, originrotoroutlet)
ax.set_xlim([0, 1000])
ax.set_ylim([-500, 500])
plt.show()
