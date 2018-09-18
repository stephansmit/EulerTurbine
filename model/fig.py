import sys
from radial_turbine import RadialTurbine
import  matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import product
import numpy as np
from multiprocessing import cpu_count, Pool
import warnings
warnings.filterwarnings('ignore')
from matplotlib import rc
import os
rc('text', usetex=True)
from name_dict import namedictionary
import imageio
def filter_df_2cols(df, c1, f1,c2, f2):
    return df[
    (abs(df[c1]-f1)<1e-10) &
    (abs(df[c2]-f2)<1e-10) 
    ]

def filter_df_3cols(df, c1, f1,c2, f2,c3, f3):
    return df[
    (abs(df[c1]-f1)<1e-10) &
    (abs(df[c2]-f2)<1e-10) &
    (abs(df[c3]-f3)<1e-10) 
    ]

def scale_variables(df):
    columns = {'work': 1e-3, "state2_area":1e6, "state3_area":1e6, "state3_height":1e3, "state2_height":1e3, "state3_thermodynamic_static_P":1e-5, }
    for col in columns:
        df[col]=df[col]*columns[col]
    return df

def plot_contour(df, cx, cy, cz, c1, f1, c2, f2, c3, f3, figpath, figname, fontsize):
    df_filter = filter_df_3cols(df, c1, f1, c2, f2, c3, f3)
    hdf = df_filter.groupby([cx,cy])[[cz]].mean()
    hdfreset = hdf.reset_index()
    hdfreset.columns = [cx, cy, cz]
    hdfpivot=hdfreset.pivot(cx,cy)
    X=hdfpivot.columns.levels[1].values
    Y=hdfpivot.index.values
    Z=hdfpivot.values
    Xi,Yi = np.meshgrid(X, Y)
    fig, ax = plt.subplots(figsize=(12,9))
    CS = ax.contourf(Yi, Xi, Z, 20,alpha=0.7, cmap=plt.cm.jet)
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel(" ".join([namedictionary[cz]['symbol'],namedictionary[cz]['units']]), fontsize=fontsize)
    ax.set_title(namedictionary[cz]["longname"], fontsize=fontsize)
    ax.set_title(namedictionary[cz]["longname"]+'\n'+
                namedictionary[c1]['symbol']+' = '+ "{0:10.2f}".format(f1) + ' ' +namedictionary[c1]['units']+' , '+
                namedictionary[c2]['symbol']+' = '+ "{0:10.2f}".format(f2) + ' ' +namedictionary[c2]['units']+' , '+
                namedictionary[c3]['symbol']+' = '+ "{0:10.2f}".format(f3) + ' ' +namedictionary[c3]['units'], 
                fontsize=fontsize)
    ax.set_xlabel(" ".join([namedictionary[cx]['symbol'],namedictionary[cx]['units']]), fontsize=fontsize)
    ax.set_ylabel(" ".join([namedictionary[cy]['symbol'],namedictionary[cy]['units']]), fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)    
    suffix=".png"
    plt.tight_layout();plt.savefig(os.path.join(figpath,figname+suffix))
    plt.close()
    return 

def plot_contour_gif(df, cx, cy, cz, tc1, c1, f1, c2, f2, figpath, figname, fontsize):
    time_vec = df[tc1].sort_values().unique()
    filenames = []
    maxc = df
    df_filter_tmp = filter_df_2cols(df, c1, f1, c2, f2)
    cmax = df_filter_tmp[cz].max()
    cmin = df_filter_tmp[cz].min()
    xmax = df_filter_tmp[cx].max()
    xmin = df_filter_tmp[cx].min()
    ymax = df_filter_tmp[cy].max()
    ymin = df_filter_tmp[cy].min()

        
    for index, i in enumerate(time_vec):
        df_filter = filter_df_3cols(df, tc1, i, c1, f1, c2, f2)
        hdf = df_filter.groupby([cx,cy])[[cz]].mean()
        hdfreset = hdf.reset_index()
        hdfreset.columns = [cx, cy, cz]
        hdfpivot=hdfreset.pivot(cx,cy)
        X=hdfpivot.columns.levels[1].values
        Y=hdfpivot.index.values
        Z=hdfpivot.values
        Xi,Yi = np.meshgrid(X, Y)
        fig, ax = plt.subplots(figsize=(12,9))
        CS = ax.contourf(Yi, Xi, Z,alpha=0.7, cmap=plt.cm.jet, levels=np.linspace(cmin,cmax,40))
        cbar = plt.colorbar(CS)
        ax.annotate(" ".join([namedictionary[tc1]['symbol'],"=","{0:10.2f}".format(i),namedictionary[tc1]['units']]) ,
            xy=(675, 750), xycoords='figure pixels',fontsize=21)
        #print(" ".join([namedictionary[tc1]['symbol'],"=","{0:10.2f}".format(i),namedictionary[tc1]['units']])) 
        cbar.ax.set_ylabel(" ".join([namedictionary[cz]['symbol'],namedictionary[cz]['units']]), fontsize=fontsize)
        ax.set_title(namedictionary[cz]["longname"]+'\n'+
                namedictionary[c1]['symbol']+' = '+ "{0:10.2f}".format(f1) + ' ' +namedictionary[c1]['units']+' , '+
                namedictionary[c2]['symbol']+' = '+ "{0:10.2f}".format(f2) + ' ' +namedictionary[c2]['units'], 
                fontsize=fontsize)
        ax.set_xlabel(" ".join([namedictionary[cx]['symbol'],namedictionary[cx]['units']]), fontsize=fontsize)
        ax.set_ylabel(" ".join([namedictionary[cy]['symbol'],namedictionary[cy]['units']]), fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize)
        suffix=str(index)+".png"
        plt.ylim(ymin, ymax) 
        plt.xlim(xmin, xmax) 
        plt.tight_layout();plt.savefig(os.path.join(figpath,figname+suffix))
        plt.close()
        filenames.append(os.path.join(figpath,figname+suffix))
        images = []
        suffix2=".gif"
    for filename in filenames:
        images.append(imageio.imread(filename))
    imageio.mimsave(os.path.join(figpath,figname+suffix2), images,duration=0.5)
    #for i in filenames:
    #    os.system("rm -rf "+i)

        


if __name__=="__main__":
    figpath = "/home/azureuser/Documents/EulerTurbine/results/variableh3/figures/"
    moviepath = "/home/azureuser/Documents/EulerTurbine/results/variableh3/movies/"
    dfturbine = pd.read_csv("../data/EulerTurbineData2.txt", sep='\t')
    dfturbine = scale_variables(dfturbine)
    
    
    alpha1_vec = dfturbine['state1_kinematic_alpha'].unique()
    h2_vec = dfturbine['state2_height'].unique()
    
    alpha2_vec = dfturbine['state2_kinematic_alpha'].unique()
    omega_vec = dfturbine['omega'].sort_values().unique()
    R_vec = dfturbine['R'].sort_values().unique()
    h2_vec = dfturbine['state2_height'].sort_values().unique()

    a3_vec = dfturbine['state3_area'].sort_values().unique()
    alpha2_vec = dfturbine['state2_kinematic_alpha'].unique()
    omega_vec = dfturbine['omega'].sort_values().unique()
    zeta3_vec = dfturbine['state3_kinematic_zeta'].sort_values().unique()

    dfoptimized = dfturbine[dfturbine['optimized']]
    dfnotoptimized = dfturbine[dfturbine['optimized']==False]
    h3_vec = dfnotoptimized['state3_height'].sort_values().unique()
    #plot_contour_gif(dfoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "state3_height","state2_kinematic_alpha", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state3_kinematic_zeta", zeta3_vec[0],
    #        moviepath, "state3_height_alpha2", 21
    #        )
    #plot_contour_gif(dfoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "work","state2_kinematic_alpha", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state3_kinematic_zeta", zeta3_vec[0],
    #        moviepath, "work_alpha2", 21
    #        ) 
    #plot_contour_gif(dfoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "R","state2_kinematic_alpha", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state3_kinematic_zeta", zeta3_vec[0],
    #        moviepath, "R_alpha2", 21
    #        ) 
    #plot_contour_gif(dfoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "state3_height","state3_kinematic_zeta", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state2_kinematic_alpha", alpha2_vec[0],
    #        moviepath, "state3_height_zeta3", 21
    #        )
    #plot_contour_gif(dfoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "work","state3_kinematic_zeta", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state2_kinematic_alpha", alpha2_vec[0],
    #        moviepath, "work_zeta3", 21
    #        ) 
    #plot_contour_gif(dfoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "R","state3_kinematic_zeta", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state2_kinematic_alpha", alpha2_vec[0],
    #        moviepath, "R_zeta3", 21
    #        ) 
    #plot_contour_gif(dfoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "state3_kinematic_w_theta","state3_kinematic_zeta", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state2_kinematic_alpha", alpha2_vec[0],
    #        moviepath, "state3_kinematic_w_theta_zeta3", 21
    #        ) 
    #plot_contour_gif(dfoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "state3_kinematic_w_r","state3_kinematic_zeta", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state2_kinematic_alpha", alpha2_vec[0],
    #        moviepath, "state3_kinematic_w_r_zeta3", 21
    #        ) 
    #plot_contour_gif(dfnotoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "R","state3_kinematic_zeta", 
    #        "state3_height", h3_vec[0],
    #        "state2_kinematic_alpha", alpha2_vec[0],
    #        moviepath, "notoptimized_R_zeta3", 21
    #        ) 
    #plot_contour_gif(dfnotoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "state3_kinematic_beta","state3_height", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state3_kinematic_zeta", zeta3_vec[0],
    #        moviepath, "notoptimized_state3_kinematic_beta_height3", 21
    #        ) 
    #plot_contour_gif(dfnotoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "state3_kinematic_alpha","state3_height", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state3_kinematic_zeta", zeta3_vec[0],
    #        moviepath, "notoptimized_state3_kinematic_alpha_height3", 21
    #        ) 
    #plot_contour_gif(dfnotoptimized, 
    #        "omega", "state3_thermodynamic_static_P", "work","state3_height", 
    #        "state1_kinematic_alpha", alpha1_vec[0],
    #        "state3_kinematic_zeta", zeta3_vec[0],
    #        moviepath, "notoptimized_work_height3", 21
    #        ) 
    plot_contour_gif(dfnotoptimized, 
            "omega", "state3_thermodynamic_static_P", "R","state2_kinematic_alpha", 
            "state1_kinematic_alpha", alpha1_vec[0],
            "state3_kinematic_zeta", zeta3_vec[0],
            moviepath, "notoptimized_R_alpha2", 21
            ) 
    plot_contour_gif(dfnotoptimized, 
            "omega", "state3_thermodynamic_static_P", "state2_thermodynamic_static_H","state2_kinematic_alpha", 
            "state1_kinematic_alpha", alpha1_vec[0],
            "state3_kinematic_zeta", zeta3_vec[0],
            moviepath, "notoptimized_state2_thermodynamic_static_H_alpha2", 21
            ) 
    plot_contour_gif(dfnotoptimized, 
            "omega", "state3_thermodynamic_static_P", "state2_kinematic_c_mag","state2_kinematic_alpha", 
            "state1_kinematic_alpha", alpha1_vec[0],
            "state3_kinematic_zeta", zeta3_vec[0],
            moviepath, "notoptimized_state2_kinematic_c_mag_alpha2", 21
            ) 
    plot_contour_gif(dfnotoptimized, 
            "omega", "state3_thermodynamic_static_P", "state2_kinematic_c_r","state2_kinematic_alpha", 
            "state1_kinematic_alpha", alpha1_vec[0],
            "state3_kinematic_zeta", zeta3_vec[0],
            moviepath, "notoptimized_state2_kinematic_c_r_alpha2", 21
            ) 
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "work","state3_area", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             moviepath, "work_area3", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "state3_kinematic_w_mag","state3_area", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             moviepath, "wmag3_area3", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "state3_kinematic_w_r","state3_area", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             moviepath, "wr3_area3", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "state3_kinematic_c_mag","state3_area", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             moviepath, "cmag3_area3", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "state3_kinematic_alpha","state3_area", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             moviepath, "alpha3_area3", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "work","state3_thermodynamic_static_P", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             moviepath, "work_P3", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "state3_kinematic_alpha","state3_thermodynamic_static_P", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             moviepath, "alpha3_P3", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "work","state2_kinematic_alpha", 
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             "state3_height", h3_vec[-1],
    #             moviepath, "work_alpha2", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "state3_kinematic_alpha","state2_kinematic_alpha", 
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             "state3_height", h3_vec[-1],
    #             moviepath, "alpha3_alpha2", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "state2_area","state2_kinematic_alpha", 
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             "state3_height", h3_vec[-1],
    #             moviepath, "area2_alpha2", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "state2_kinematic_c_theta","state2_kinematic_alpha", 
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             "state3_height", h3_vec[-1],
    #             moviepath, "ctheta2_alpha2", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "state2_kinematic_c_mag","state2_kinematic_alpha", 
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             "state3_height", h3_vec[-1],
    #             moviepath, "cmag2_alpha2", 21
    #             )
    #plot_contour_gif(dfturbine, 
    #             "omega", "R", "state2_massflow","state2_kinematic_alpha", 
    #             "state3_thermodynamic_static_P", P3_vec[0],
    #             "state3_height", h3_vec[-1],
    #             moviepath, "m2_alpha2", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "work", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "work", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state2_kinematic_c_mag", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "c2mag", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state2_kinematic_c_r", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "c2r", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state2_kinematic_c_theta", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "c2theta", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state3_kinematic_c_theta", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "c3theta", 21
    #             )
    #
    #plot_contour(dfturbine, 
    #             "omega", "R", "state2_kinematic_beta", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "beta2", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state3_kinematic_c_mach", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "c3mach", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state3_kinematic_c_mag", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "c3mag", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state3_thermodynamic_static_A", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "a3", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state2_thermodynamic_static_A", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "a2", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state3_kinematic_c_theta", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "c3theta", 21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state3_kinematic_w_mag", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "w3mag",  21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state3_kinematic_alpha", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "alpha3",  21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state3_kinematic_beta", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "beta3",  21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state3_kinematic_w_mach", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "w3mach",  21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state3_kinematic_c_mach", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "c3mach",  21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state2_kinematic_w_mach", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "w2mach",  21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state2_kinematic_c_mach", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "c2mach",  21
    #             )
    #plot_contour(dfturbine, 
    #             "omega", "R", "state2_area", 
    #             "state2_kinematic_alpha", alpha2_vec[0],
    #             "state3_height", h3_vec[-1],
    #             "state3_thermodynamic_static_P", P3_vec[0], 
    #             figpath, "c2area",  21
    #             )
