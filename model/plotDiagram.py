import numpy as np
import CoolProp.CoolProp as CP
import CoolProp.Plots as CPP
import matplotlib.pyplot as plt
import pandas as pd
import warnings
from CoolProp.Plots.SimpleCycles import StateContainer, StatePoint
import CoolProp
from CoolProp.Plots import PropertyPlot
from CoolProp.CoolProp import PropsSI
from radial_turbine import RadialTurbine
from CoolProp.CoolProp import PhaseSI
plt.rc('text', usetex=True)
import pandas as pd
warnings.filterwarnings('ignore')

def create_Hs_diagram(H, S, H1, H01, S1,H2, H02, S2,H2s, S2s, H3,H03, S3,H3s, S3s,H4,H04, S4, P1, P2, P3,P4,dome1h, dome2h, dome1, dome2,  criticalisobar):
    #creating isobars
    stator_inlet = PropsSI("H", "P", P1, "S", S, 'Toluene')*1e-3
    stator_outlet = PropsSI("H", "P", P2, "S", S, 'Toluene')*1e-3
    rotor_outlet = PropsSI("H", "P", P3, "S", S, 'Toluene')*1e-3
    condenser_inlet = PropsSI("H", "P", P4, "S", S, 'Toluene')*1e-3


    #creating contour
    dftemp = pd.DataFrame(H, columns=['H'])
    dfentr = pd.DataFrame(S, columns=['S'])
    dftemp['key'] = 1
    dfentr['key'] = 1
    df = dftemp.merge(dfentr, how='outer')[['S', 'H']]
    df['Z'] = df.apply(lambda x: PropsSI('Z', 'H', x['H'], 'S', x['S'], 'Toluene'), axis=1)
    df['phase'] = df.apply(lambda x: PropsSI('Phase', 'H', x['H'], 'S', x['S'],'Toluene'), axis=1)
    gas = df['phase'] == CP.get_phase_index('phase_gas')
    supercritical_gas = df['phase'] == CP.get_phase_index('phase_supercritical_gas')
    df_filter = df[gas | supercritical_gas][['H', "S", "Z"]]
    df_filter['H']=df_filter['H']*1e-3
    df_filter['S']=df_filter['S']*1e-3
    hdfpivot = df_filter.pivot('H', 'S')
    X = hdfpivot.columns.levels[1].values
    Y = hdfpivot.index.values
    Z = hdfpivot.values
    Xi, Yi = np.meshgrid(X, Y)

    fig, ax = plt.subplots(figsize=(10, 7))

    # isobars
    ax.plot(S*1e-3, H01*np.ones(np.size(S))*1e-3, color='black', linestyle=':')
    fig.text(0.15, 0.77, "$\displaystyle h_{01}$",fontsize=20)
    ax.plot(S*1e-3, H03*np.ones(np.size(S))*1e-3, color='black', linestyle=':')
    fig.text(0.15, 0.5, "$\displaystyle h_{03}$",fontsize=20)


    #isobars
    ax.plot(S*1e-3, stator_inlet, color='black', linestyle='--')
    fig.text(0.3, 0.65, "{:4.2f}".format(P1 / 1e5) + ' [Bar]')
    ax.plot(S*1e-3, stator_outlet, color='black', linestyle='--')
    fig.text(0.3, 0.42, "{:4.2f}".format(P2 / 1e5) + ' [Bar]')
    ax.plot(S*1e-3, rotor_outlet, color='black', linestyle='--')
    fig.text(0.3, 0.29, "{:4.2f}".format(P3 / 1e5) + ' [Bar]')
    ax.plot(S*1e-3, condenser_inlet, color='black', linestyle='--')
    fig.text(0.3, 0.33, "{:4.2f}".format(P4 / 1e5) + ' [Bar]')
    # ax.plot(S, condenser_inlet, color='black', linestyle='--')

    # dome
    ax.plot(dome1*1e-3, dome1h*1e-3, color='black')
    ax.plot(dome2*1e-3, dome2h*1e-3, color='black')

    # process
    ax.plot([S1*1e-3,S2*1e-3, S3*1e-3,S4*1e-3], [H1*1e-3,H2*1e-3,H3*1e-3,H4*1e-3], color='black', marker="X")
    #isentropic process
    ax.plot([S1*1e-3,S2s*1e-3], [H1*1e-3,H2s*1e-3], color='black', marker="o")
    ax.plot([S2*1e-3,S3s*1e-3], [H2*1e-3,H3s*1e-3], color='black', marker="o")

    # compressibility
    fig = ax.contourf(Xi, Yi, Z, 25, alpha=0.7, cmap=plt.cm.jet)
    ax.set_xlim([.850, 1.5])
    ax.ticklabel_format(fontsize=20)
    ax.set_ylim([200, 700])
    ax.set_xlabel(r'Entropy [KJ/kgK]', fontsize=15)
    ax.set_ylabel(r'Enthalpy [KJ/kg]', fontsize=15)
    cb = plt.colorbar(fig)
    cb.set_label(r"Compressibility factor [-]", fontsize=15)
    return fig,ax


def create_Ts_diagram(T, S, T1, S1, T3, S3, P1, P2, P3, dome1, dome2, criticalisobar):
    stator_inlet = PropsSI("T", "P", P1, "S", S, 'Toluene')
    stator_outlet = PropsSI("T", "P", P2, "S", S, 'Toluene')
    rotor_outlet = PropsSI("T", "P", P3, "S", S, 'Toluene')
    dftemp = pd.DataFrame(T, columns=['T'])
    dfentr = pd.DataFrame(S, columns=['S'])
    dftemp['key'] = 1
    dfentr['key'] = 1
    df = dftemp.merge(dfentr, how='outer')[['S', 'T']]
    df['Z'] = df.apply(lambda x: PropsSI('Z', 'T', x['T'], 'S', x['S'], 'Toluene'), axis=1)
    df['phase'] = df.apply(lambda x: PropsSI('Phase', 'T', x['T'], 'S', x['S'], 'Toluene'), axis=1)

    gas = df['phase'] == CP.get_phase_index('phase_gas')
    supercritical_gas = df['phase'] == CP.get_phase_index('phase_supercritical_gas')
    df_filter = df[gas | supercritical_gas][['T', "S", "Z"]]
    hdfpivot = df_filter.pivot('T', 'S')
    X = hdfpivot.columns.levels[1].values
    Y = hdfpivot.index.values
    Z = hdfpivot.values
    Xi, Yi = np.meshgrid(X, Y)

    fig, ax = plt.subplots(figsize=(10, 7))

    # isobars
    ax.plot(S, stator_inlet, color='black', linestyle='--')
    fig.text(0.3, 0.80, "{:4.2f}".format(P1 / 1e5) + ' [Bar]')
    ax.plot(S, stator_outlet, color='black', linestyle='--')
    fig.text(0.3, 0.42, "{:4.2f}".format(P2 / 1e5) + ' [Bar]')
    ax.plot(S, rotor_outlet, color='black', linestyle='--')
    fig.text(0.3, 0.29, "{:4.2f}".format(P3 / 1e5) + ' [Bar]')
    # ax.plot(S, condenser_inlet, color='black', linestyle='--')

    # dome
    ax.plot(dome1, T, color='black')
    ax.plot(dome2, T, color='black')

    # process
    ax.plot([S1, S3], [T1, T3], color='black', marker="X")

    # compressibility
    fig = ax.contourf(Xi, Yi, Z, 25, alpha=0.7, cmap=plt.cm.jet)
    ax.set_xlim([850, 1500])
    ax.ticklabel_format(fontsize=20)
    ax.set_ylim([273.15, 600])
    ax.set_xlabel(r'Entropy [J/kgK]', fontsize=15)
    ax.set_ylabel(r'Temperature [K]', fontsize=15)
    cb = plt.colorbar(fig)
    cb.set_label(r"Compressibility factor [-]", fontsize=15)

    return fig, ax

if __name__ =="__main__":
    # rotationalspeed
    omega = 525
    etats = 0.5
    # heights
    h1 = 8.897e-3
    h2 = 8.897e-3
    h3 = 12.e-3

    # radial postions
    r1 = .1680
    r2 = 0.115
    r3 = 0.0873

    # inlet conditions
    P1 = 31.95e5
    P2 =1e5
    T1 = 587.65
    C1 = 0.5
    alpha1 = 0

    # inlet/rotor interface conditions
    alpha2 = 85.2
    DOR = 0.295658
    zeta3 = 0
    # outlet conditions
    P3 = 0.08926012e5
    init = 480

    P4 = 0.2e5
    massflow = .18
    optimized = True
    turbine = RadialTurbine(omega, h1, h2, h3, r1, r2, r3, P1, T1, alpha1, P2, zeta3, P3, P4, massflow, optimized)

    #stator inlet conditions
    P1 = turbine.state1.thermodynamic.static.P
    T1 = turbine.state1.thermodynamic.static.T
    S1 = turbine.state1.thermodynamic.static.S #PropsSI("S", "P", P1, "T", T1, 'Toluene')
    H1 = turbine.state1.thermodynamic.static.H # PropsSI("H", "P", P1, "T", T1, 'Toluene')
    H01 = turbine.state1.thermodynamic.total.H

    #stator outlet conditions
    P2 = turbine.state2.thermodynamic.static.P
    T2 = turbine.state2.thermodynamic.static.T
    S2 = turbine.state2.thermodynamic.static.S#PropsSI("S", "P", P2, "T", T2, 'Toluene')
    H2 = turbine.state2.thermodynamic.static.H#PropsSI("H", "P", P2, "T", T2, 'Toluene')
    H02 = turbine.state2.thermodynamic.total.H#PropsSI("H", "P", P2, "T", T2, 'Toluene')

    #stator outlet conditions
    P2s = turbine.state2.thermodynamic.isentropic_static.P
    T2s = turbine.state2.thermodynamic.isentropic_static.T
    S2s = turbine.state2.thermodynamic.isentropic_static.S#PropsSI("S", "P", P2, "T", T2, 'Toluene')
    H2s = turbine.state2.thermodynamic.isentropic_static.H#PropsSI("H", "P", P2, "T", T2, 'Toluene')
    H02s = turbine.state2.thermodynamic.isentropic_total.H#PropsSI("H", "P", P2, "T", T2, 'Toluene')

    #rotor outlet condition
    P3 = turbine.state3.thermodynamic.static.P#.1e5
    T3 = turbine.state3.thermodynamic.static.T#425
    S3 = turbine.state3.thermodynamic.static.S#PropsSI("S", "P", P3, "T", T3, 'Toluene')
    H3 = turbine.state3.thermodynamic.static.H#PropsSI("H", "P", P3, "T", T3, 'Toluene')
    #C3 = 300
    H03 = turbine.state3.thermodynamic.total.H
    #stator outlet conditions
    P3s = turbine.state3.thermodynamic.isentropic_static.P
    T3s = turbine.state3.thermodynamic.isentropic_static.T
    S3s = turbine.state3.thermodynamic.isentropic_static.S#PropsSI("S", "P", P2, "T", T2, 'Toluene')
    H3s = turbine.state3.thermodynamic.isentropic_static.H#PropsSI("H", "P", P2, "T", T2, 'Toluene')
    H03s = turbine.state3.thermodynamic.isentropic_total.H#PropsSI("H", "P", P2, "T", T2, 'Toluene')

    #condensor inlet
    P4 = turbine.state4.thermodynamic.static.P
    T4 = turbine.state4.thermodynamic.static.T
    S4 = turbine.state4.thermodynamic.static.S#PropsSI("S", "P", P4, "T", T4, 'Toluene')
    H4 = turbine.state4.thermodynamic.static.H#PropsSI("H", "P", P4, "T", T4, 'Toluene')
    H04 = turbine.state4.thermodynamic.total.H#PropsSI("H", "P", P4, "T", T4, 'Toluene')

    #domain
    nT = 100
    Tmin = 273.15
    Tmax = 650
    T = np.linspace(Tmin, Tmax, nT)

    nS = 100
    Smin = 850
    Smax = 1500
    S = np.linspace(Smin, Smax, nS)

    nH = nT
    Hmin = 200*1e3#PropsSI("H", "T", Tmin, "S", Smin, 'Toluene')
    Hmax = 700*1e3 #psSI("H", "T", Tmax, "S", Smax, 'Toluene')
    H = np.linspace(Hmin, Hmax, nH)

    print(Hmin, Hmax)



    condenser_inlet = PropsSI("T", "P", P4, "S", S, 'Toluene')

    # diagram bars
    criticalisobar = PropsSI("T", "P", CP.PropsSI("Pcrit", "Toluene"), "S", S, 'Toluene')
    dome1 = PropsSI("S", "T", T, "Q", 1., 'Toluene')
    dome2 = PropsSI("S", "T", T, "Q", 0., 'Toluene')

    dome1 = PropsSI("S", "T", T, "Q", 1., 'Toluene')
    dome2 = PropsSI("S", "T", T, "Q", 0., 'Toluene')

    dome1h=[]
    dome2h=[]
    for index, element in enumerate(dome1):
        try:
            dome1h.append(PropsSI("H", "T", T[index], "S", element, 'Toluene'))
        except:
            dome1h.append(float('Inf'))
    for index, element in enumerate(dome2):
        try:
            dome2h.append(PropsSI("H", "T", T[index], "S", element, 'Toluene'))
        except:
            dome2h.append(float('Inf'))
        #fig, ax = create_Ts_diagram(T, S, T1, S1, T3, S3, stator_inlet, stator_outlet, rotor_outlet, dome1, dome2, criticalisobar)
    dome1h = np.asarray(dome1h)
    dome2h = np.asarray(dome2h)

    fig, ax = create_Hs_diagram(H, S, H1,H01, S1,H2,H02,S2,H2s,S2s, H3,H03, S3,H3s,S3s, H4,H04, S4,P1, P2, P3, P4,dome1h, dome2h,dome1, dome2,criticalisobar)

    plt.show(fig)

    #
    # test = PropsSI("S", "T", temperatures, "Q", 1., 'Toluene')
    # test2 = PropsSI("S", "T", temperatures, "Q", 0., 'Toluene')
