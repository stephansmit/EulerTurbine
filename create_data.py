from radial_turbine import RadialTurbine
import  matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import product
import numpy as np
from multiprocessing import cpu_count, Pool
import warnings
import sys

def getRadialTurbineInfo(x):
    omega=x['omega']
    h1=x['h1']
    h2=x['h2']
    h3=x['h3']
    r1=x['r1']
    r2=x['r2']
    r3=x['r3']
    P1=x['P1']
    T1=x['T1']
    C1=x['C1']
    alpha1=x['alpha1']
    alpha2=x['alpha2']
    P3=x['P3']
    DOR=x['DOR']
    try:
        data = RadialTurbine(omega, h1, h2, h3, r1, r2, r3, P1, T1, C1, alpha1, alpha2, P3, DOR).get_turbine_info()
    except:
        data = empty_dict
    return pd.Series(data)

def split(df, partitions):
    return np.array_split(df, partitions, axis=0)

def parallizable_getRadialTurbineInfo(df):
    result = df.apply(getRadialTurbineInfo,axis=1)
    return result

def parallelized_getRadialTurbineInfo(df, cores):
    partitions = cores
    dfs = split(df, partitions)
    pool = Pool(cores)
    data = pd.concat(pool.map(parallizable_getRadialTurbineInfo, dfs))
    pool.close()
    pool.join()
    return data

if __name__=="__main__":
    
    
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


    ncores = int(sys.argv[1])
    df_DOR = pd.DataFrame()
    df_h2 = pd.DataFrame()
    df_alpha2 = pd.DataFrame()
    df_omega = pd.DataFrame()
    df_h3 = pd.DataFrame()

    df_DOR['DOR'] = np.linspace(0,0.5,5)
    df_DOR['key'] = 1
    df_h2['h2'] = np.linspace(5e-3, 9e-3,5)
    df_h2['key']=1
    df_alpha2['alpha2'] = np.linspace(70, 80,5)
    df_alpha2['key']=1
    df_omega['omega']=np.linspace(2600,3600, 10)
    df_omega['key']=1
    df_h3['h3'] = np.linspace(1e-2,5e-2, 5)
    df_h3['key']=1

    df1 = pd.merge(df_DOR, df_omega, on='key')
    df2 = pd.merge(df_alpha2, df1, on='key')
    df3 = pd.merge(df_h2, df2, on='key')
    df = pd.merge(df_h3, df3, on='key')

    df['h1'] = 0.002
    df['r1']=.210
    df['r2']=0.115
    df['r3']=0.090
    df['P1']=32e5
    df['T1']=525
    df['C1']=0.5
    df['alpha1']=45
    df['P3'] = .2e5
    
    dfturbine = parallelized_getRadialTurbineInfo(df, ncores)
    dfplot = dfturbine.dropna()
    dfplot.to_csv("EulerTurbineData.txt", sep='\t', index=False)
