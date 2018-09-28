from numpy import pi, sin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import pandas as pd
import math
import ipywidgets as widgets
from IPython.display import display

plt.rc('text', usetex=True)



def scale_variables(df):
    columns = {'work': 1e-3, "state2_area":1e6, "state3_area":1e6, "state3_height":1e3, "state2_height":1e3, "state2_thermodynamic_static_P":1e-5,"state3_thermodynamic_static_P":1e-5,"state3_thermodynamic_total_H":1e-3, "state3_thermodynamic_static_H":1e-3, "state2_thermodynamic_static_H":1e-3}
    for col in columns:
        df[col]=df[col]*columns[col]
    return df

if __name__=="__main__":
    # read info
    dfturbine = pd.read_csv("/home/azureuser/Documents/EulerTurbine/data/EulerTurbineData.txt", sep='\t')

    # create extra parameters
    dfturbine['tsefficiency'] = (dfturbine['state1_thermodynamic_total_H'] - dfturbine['state3_thermodynamic_total_H']) / \
                                (dfturbine['state1_thermodynamic_total_H'] - dfturbine['state3_thermodynamic_static_H'])
    dfturbine = scale_variables(dfturbine)
    dfturbine['omegaHz'] = dfturbine['omega'].apply(lambda x: x / (2 * math.pi))
    dfturbine['centrifugalterm'] = 0.5 * (dfturbine['state2_kinematic_u_mag'] * dfturbine['state2_kinematic_u_mag'] -
                                          dfturbine['state3_kinematic_u_mag'] * dfturbine['state3_kinematic_u_mag']) / 1000.
    dfturbine['state3_kinematic_ske'] = 0.5 * dfturbine['state3_thermodynamic_static_D'] * dfturbine[
        'state3_kinematic_c_mag'] * dfturbine['state3_kinematic_c_mag'] / 1000.0
    dfturbine['state2_kinematic_ske'] = 0.5 * dfturbine['state3_thermodynamic_static_D'] * dfturbine[
        'state2_kinematic_c_mag'] * dfturbine['state2_kinematic_c_mag'] / 1000.0
    dfturbine['ratio'] = dfturbine['centrifugalterm'] / (
    dfturbine['state2_thermodynamic_static_H'] - dfturbine['state3_thermodynamic_static_H'])

    # vectors used for filtering
    alpha1_vec = dfturbine['state1_kinematic_alpha'].unique()
    h2_vec = dfturbine['state2_height'].unique()
    alpha2_vec = dfturbine['state2_kinematic_alpha'].unique()
    omega_vec = dfturbine['omegaHz'].sort_values().unique()
    R_vec = dfturbine['R'].sort_values().unique()
    h2_vec = dfturbine['state2_height'].sort_values().unique()
    P2_vec = dfturbine['state2_thermodynamic_static_P'].sort_values().unique()
    a3_vec = dfturbine['state3_area'].sort_values().unique()
    P3_vec = dfturbine['state3_thermodynamic_static_P'].sort_values().unique()

    # optimized versus non-optimal
    dfnotopti = dfturbine[dfturbine['optimized'] == False]
    dfopti = dfturbine[dfturbine['optimized']]

    h3_vec = dfopti['state3_height'].unique()
    h3_vec2 = dfnotopti['state3_height'].sort_values().unique()
    zeta3_vec = dfopti['state3_kinematic_zeta'].unique()
    zeta3_vec2 = dfnotopti['state3_kinematic_zeta'].unique()

    h3 = h3_vec2[3]
    alpha2 = alpha2_vec[0]
    zeta3 = zeta3_vec2[0]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    axis_color = 'lightgoldenrodyellow'
    amp_0 = 5
    freq_0 = 3

    # Define an axes area and draw a slider in it
    amp_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axis_color)
    amp_slider = Slider(amp_slider_ax, '$\displaystyle P_3$', P3_vec[0], P3_vec[-1], valinit=P3_vec[len(P3_vec)/2])

    selection_slider = widgets.SelectionSlider(
        options=['scrambled', 'sunny side up', 'poached', 'over easy'],
        value='sunny side up',
        description='I like my eggs ...',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True
    )
    display(selection_slider)

    #amp_slider.on_changed(sliders_on_changed)
    plt.show()
