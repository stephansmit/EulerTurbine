import math
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


def filter_df_4cols(df, c1, f1,c2, f2,c3, f3, c4, f4):
    return df[
    (abs(df[c1]-f1)<1e-10) &
    (abs(df[c2]-f2)<1e-10) &
    (abs(df[c3]-f3)<1e-10) &
    (abs(df[c4]-f4)<1e-10) 
    ]

def filter_df_5cols(df, c1, f1,c2, f2,c3, f3, c4, f4, c5, f5):
    return df[
    (abs(df[c1]-f1)<1e-10) &
    (abs(df[c2]-f2)<1e-10) &
    (abs(df[c3]-f3)<1e-10) &
    (abs(df[c4]-f4)<1e-10) & 
    (abs(df[c5]-f5)<1e-10) 
    ]
def filter_df_6cols(df, c1, f1,c2, f2,c3, f3, c4, f4, c5, f5, c6, f6):
    return df[
    (abs(df[c1]-f1)<1e-10) &
    (abs(df[c2]-f2)<1e-10) &
    (abs(df[c3]-f3)<1e-10) &
    (abs(df[c4]-f4)<1e-10) & 
    (abs(df[c5]-f5)<1e-10) & 
    (abs(df[c6]-f6)<1e-10) 
    ]

def plot_xy(df,cx, cy, c1, f1,c2, f2,c3, f3, c4, f4,c5,f5, figpath, figname, fontsize):
    df_filter = filter_df_5cols(df, c1, f1, c2, f2, c3, f3,c4,f4,c5,f5)

    fig, ax = plt.subplots(figsize=(12,9))
    df_filter.plot(x=cx, y=cy, ax=ax, legend=None)
    ax.set_xlabel(" ".join([namedictionary[cx]['symbol'],namedictionary[cx]['units']]), fontsize=fontsize)
    ax.set_ylabel(" ".join([namedictionary[cy]['symbol'],namedictionary[cy]['units']]), fontsize=fontsize)
    ax.set_title(namedictionary[cy]["longname"]+'\n'+
                namedictionary[c1]['symbol']+' = '+ "{0:10.2f}".format(f1) + ' ' +namedictionary[c1]['units']+' , '+
                namedictionary[c2]['symbol']+' = '+ "{0:10.2f}".format(f2) + ' ' +namedictionary[c2]['units']+' , '+
                namedictionary[c3]['symbol']+' = '+ "{0:10.2f}".format(f3) + ' ' +namedictionary[c3]['units']+' , '+ 
                namedictionary[c4]['symbol']+' = '+ "{0:10.2f}".format(f4) + ' ' +namedictionary[c4]['units']+' , '+ 
                namedictionary[c5]['symbol']+' = '+ "{0:10.2f}".format(f5) + ' ' +namedictionary[c5]['units'], 
                fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    suffix=".png"
    plt.tight_layout();plt.savefig(os.path.join(figpath,figname+suffix))
    plt.close()
    return

def plot_xy_group(df,cx, cy, cg1, c2, f2,c3, f3, c4, f4,c5,f5, figpath, figname, fontsize):
    df_filter = filter_df_4cols(df,  c2, f2, c3, f3,c4,f4,c5,f5)
    df_filter2 = df_filter.groupby(cg1)

    fig, ax = plt.subplots(figsize=(12,9))
    df_filter2.plot(x=cx, y=cy, ax=ax, legend=None)
    ax.set_xlabel(" ".join([namedictionary[cx]['symbol'],namedictionary[cx]['units']]), fontsize=fontsize)
    ax.set_ylabel(" ".join([namedictionary[cy]['symbol'],namedictionary[cy]['units']]), fontsize=fontsize)
    ax.set_title(namedictionary[cy]["longname"]+'\n'+
                namedictionary[c2]['symbol']+' = '+ "{0:10.2f}".format(f2) + ' ' +namedictionary[c2]['units']+' , '+
                namedictionary[c3]['symbol']+' = '+ "{0:10.2f}".format(f3) + ' ' +namedictionary[c3]['units']+' , '+ 
                namedictionary[c4]['symbol']+' = '+ "{0:10.2f}".format(f4) + ' ' +namedictionary[c4]['units']+' , '+ 
                namedictionary[c5]['symbol']+' = '+ "{0:10.2f}".format(f5) + ' ' +namedictionary[c5]['units'], 
                fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    suffix=".png"
    plt.tight_layout();plt.savefig(os.path.join(figpath,figname+suffix))
    plt.close()
    return

def scale_variables(df):
    columns = {'work': 1e-3, "state2_area":1e6, "state3_area":1e6, "state3_height":1e3, "state2_height":1e3, "state2_thermodynamic_static_P":1e-5,"state3_thermodynamic_static_P":1e-5,"state3_thermodynamic_total_H":1e-3, "state3_thermodynamic_static_H":1e-3, "state2_thermodynamic_static_H":1e-3}
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
    sys.path.append("/home/azureuser/Documents/EulerTurbine/model/")
    
    #set path to write figures
    figpath_notopti= "/home/azureuser/Documents/EulerTurbine/results/variableh/figures/notoptimized/"
    figpath_opti = "/home/azureuser/Documents/EulerTurbine/results/variableh/figures/optimized/"
    moviepath_opti = "/home/azureuser/Documents/EulerTurbine/results/variableh/movies/optimized/"
    moviepath_notopti = "/home/azureuser/Documents/EulerTurbine/results/variableh/movies/notoptimized/"

    
    #read info
    dfturbine = pd.read_csv("/home/azureuser/Documents/EulerTurbine/data/EulerTurbineData.txt", sep='\t')
    
    
    #create extra parameters
    dfturbine['tsefficiency']= (dfturbine['state1_thermodynamic_total_H'] - dfturbine['state3_thermodynamic_total_H']) / \
                                         (dfturbine['state1_thermodynamic_total_H']-dfturbine['state3_thermodynamic_static_H'])
    dfturbine = scale_variables(dfturbine)
    dfturbine['omegaHz']=dfturbine['omega'].apply(lambda x: x/(2*math.pi))
    dfturbine['centrifugalterm']=0.5*(dfturbine['state2_kinematic_u_mag']*dfturbine['state2_kinematic_u_mag'] - 
                                      dfturbine['state3_kinematic_u_mag']*dfturbine['state3_kinematic_u_mag'])/1000.
    dfturbine['state3_kinematic_ske']=0.5*dfturbine['state3_thermodynamic_static_D']*dfturbine['state3_kinematic_c_mag']*dfturbine['state3_kinematic_c_mag']/1000.0
    dfturbine['state2_kinematic_ske']=0.5*dfturbine['state3_thermodynamic_static_D']*dfturbine['state2_kinematic_c_mag']*dfturbine['state2_kinematic_c_mag']/1000.0
    dfturbine['ratio']=dfturbine['centrifugalterm']/(dfturbine['state2_thermodynamic_static_H']-dfturbine['state3_thermodynamic_static_H'])

    
    
    #vectors used for filtering
    alpha1_vec = dfturbine['state1_kinematic_alpha'].unique()
    h2_vec = dfturbine['state2_height'].unique()
    alpha2_vec = dfturbine['state2_kinematic_alpha'].unique()
    omega_vec = dfturbine['omegaHz'].sort_values().unique()
    R_vec = dfturbine['R'].sort_values().unique()
    h2_vec = dfturbine['state2_height'].sort_values().unique()
    P2_vec = dfturbine['state2_thermodynamic_static_P'].sort_values().unique()
    a3_vec = dfturbine['state3_area'].sort_values().unique()
    P3_vec = dfturbine['state3_thermodynamic_static_P'].sort_values().unique()



    #optimized versus non-optimal
    dfnotopti = dfturbine[dfturbine['optimized']==False]
    dfopti = dfturbine[dfturbine['optimized']]
    
    h3_vec = dfopti['state3_height'].unique()
    h3_vec2 = dfnotopti['state3_height'].sort_values().unique()
    zeta3_vec = dfopti['state3_kinematic_zeta'].unique()
    zeta3_vec2 = dfnotopti['state3_kinematic_zeta'].unique()

    h3=h3_vec2[3]
    alpha2=alpha2_vec[0]
    zeta3= zeta3_vec2[0]
    plotContour=True
    if plotContour:
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "tsefficiency",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_eta", 30
                      )
         plot_contour(dfnotopti,
                      "state3_thermodynamic_static_P", 
                      "state2_kinematic_alpha", 
                      "R",
                      "omegaHz", omega_vec[3],
                      "state3_height", h3_vec2[1],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      figpath_notopti,"P3alpha2_R", 30
                      )
         plot_contour(dfnotopti,
                      "state3_thermodynamic_static_P", 
                      "state2_kinematic_alpha", 
                      "state3_kinematic_beta",
                      "omegaHz", omega_vec[3],
                      "state3_height", h3_vec2[1],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      figpath_notopti,"P3alpha2_beta", 30
                      )
	 plot_contour(dfnotopti,
              "omegaHz", "state3_thermodynamic_static_P", "ratio",
              "state2_kinematic_alpha", alpha2,
              "state3_height", h3,
              "state3_kinematic_zeta", zeta3,
              figpath_notopti, "omegaP3_ratio", 30
              )
         plot_contour(dfnotopti,
                      "omegaHz", 
                      "state3_thermodynamic_static_P", 
                      "R",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_R", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_alpha",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_alpha3", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state2_thermodynamic_static_P", "state2_kinematic_w_mag",
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP2_w2mag", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state2_thermodynamic_static_P", "state2_kinematic_w_mach",
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP2_w2mach", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "work",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_work", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_beta",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_beta3", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_w_mag",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_w3mag", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_w_mach",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_w3mach", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state2_thermodynamic_static_P", "state2_kinematic_c_mag",
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP2_c2mag", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state2_thermodynamic_static_P", "state2_kinematic_c_mach",
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP2_c2mach", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "R",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_R", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_c_mag",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_c3mag", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_c_mach",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_c3mach", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_ske",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_ske3", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_thermodynamic_static_H",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_H3", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_thermodynamic_total_H",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_H03", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_thermodynamic_total_H",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_H03", 30
                      )
         plot_contour(dfnotopti,
                      "omegaHz", "state3_thermodynamic_static_P", "work",
                      "state2_kinematic_alpha", alpha2,
                      "state3_height", h3,
                      "state3_kinematic_zeta", zeta3,
                      figpath_notopti, "omegaP3_work", 30
                      )


         plot_contour(dfopti,
                  "omegaHz", "state3_thermodynamic_static_P", "work",
                  "state2_kinematic_alpha", alpha2_vec[0],
                  "state1_kinematic_alpha", alpha1_vec[0],
                  "state3_kinematic_zeta", zeta3_vec[0],
                  figpath_opti, "omegaP3_work", 30
                  )

         plot_contour(dfopti,
                      "omegaHz", "state2_kinematic_alpha", "work",
                      "state3_thermodynamic_static_P", P3_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaalpha2_work", 30
                     )
         plot_contour(dfopti,
                      "omegaHz", "state2_thermodynamic_static_P", "work",
                      "state3_thermodynamic_static_P", P3_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP2_work", 30
                     )
         plot_contour(dfopti,
                      "omegaHz", "state2_thermodynamic_static_P", "state2_kinematic_c_mag",
                      "state3_thermodynamic_static_P", P3_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP2_c2mag", 30
                     )
         plot_contour(dfopti,
                      "omegaHz", "state2_thermodynamic_static_P", "state2_kinematic_c_mach",
                      "state3_thermodynamic_static_P", P3_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP2_c2mach", 30
                     )
         plot_contour(dfopti,
                      "omegaHz", "state2_thermodynamic_static_P", "state2_kinematic_c_theta",
                      "state3_thermodynamic_static_P", P3_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP2_c2theta", 30
                     )
         plot_contour(dfopti,
                      "omegaHz", "state2_thermodynamic_static_P", "R",
                      "state3_thermodynamic_static_P", P3_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP2_R", 30
                     )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_height",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_height3", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_w_mag",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_w3mag", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_w_mach",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_w3mach", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_c_mag",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_c3mag", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_c_mach",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_c3mach", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_c_r",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_c3r", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_beta",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_beta3", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "R",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_R", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_alpha",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_alpha3", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_thermodynamic_static_D",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_D3", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_thermodynamic_static_H",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_H3", 30
                      )
         
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "centrifugalterm",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_cent", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_thermodynamic_total_H",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_H03", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "tsefficiency",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_eta", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_ske",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_ske3", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_c_mag",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_c3mag", 30
                      )
         plot_contour(dfopti,
                      "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_c_mach",
                      "state2_kinematic_alpha", alpha2_vec[0],
                      "state1_kinematic_alpha", alpha1_vec[0],
                      "state3_kinematic_zeta", zeta3_vec[0],
                      figpath_opti, "omegaP3_c3mach", 30
                      )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "work","state2_kinematic_alpha", 
                "state1_kinematic_alpha", alpha1_vec[0],
                "state3_height", h3_vec2[4],
                moviepath_notopti, "omegaP3_work_alpha2", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "R","state2_kinematic_alpha", 
                "state1_kinematic_alpha", alpha1_vec[0],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_R_alpha2", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_c_mag","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[0],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_c3mag-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_c_mag","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[0],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_c3mag-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "tsefficiency","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[0],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_efficiency-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "R","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[0],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_R-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "work","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[0],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_work-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "state3_thermodynamic_static_D","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[0],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_D3-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_w_mag","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[0],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_w3mag-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_alpha","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[0],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_alpha3-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "state3_thermodynamic_static_H","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[-1],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_H3-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_ske","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[-1],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_ske3-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "state3_thermodynamic_total_H","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[-1],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_H03-height", 21
                )
    plot_contour_gif(dfnotopti, 
                "omegaHz", "state3_thermodynamic_static_P", "state3_kinematic_w_mag","state3_height", 
                "state2_kinematic_alpha", alpha2_vec[-1],
                "state2_height", h2_vec[0],
                moviepath_notopti, "omegaP3_w3mag-height", 21
                )






    plot_contour_gif(dfopti, 
                "omegaHz", "state3_thermodynamic_static_P", "work","state2_kinematic_alpha", 
                "state1_kinematic_alpha", alpha1_vec[0],
                "state2_height", h2_vec[0],
                moviepath_opti, "omegaP3_work_alpha2", 21
                )
    plot_contour_gif(dfopti, 
                "omegaHz", "state3_thermodynamic_static_P", "R","state2_kinematic_alpha", 
                "state1_kinematic_alpha", alpha1_vec[0],
                "state2_height", h2_vec[0],
                moviepath_opti, "omegaP3_R_alpha2", 21
                )

