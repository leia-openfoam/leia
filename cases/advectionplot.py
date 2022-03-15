import pandas as pd
import numpy as np
pd.set_option("display.precision", 8)
from matplotlib import pyplot as plt
import os
from math import pi, log
from matplotlib import rcParams
rcParams["text.usetex"] = True
rcParams["figure.dpi"] = 200
rcParams["font.size"] = 18


def agglomerate_dframe(csv_filename="",study_pattern="",dframe_name=""):
    # Find all "csv_filename.csv" files in study_pattern*/ folders
    csv_files = [os.path.join(folder, csv_filename) 
                 for folder in os.listdir(os.curdir) 
                 if os.path.isfile(os.path.join(folder, csv_filename))
                 and study_pattern in folder]
    csv_files.sort()
    # Read all "csv_filename.csv" files into a pandas.DataFrame
    dframes = []
    for csv_file in csv_files:
        dframes.append(pd.read_csv(csv_file, header=0)) 
    final_dframe = pd.concat(dframes, ignore_index=True)
    final_dframe.to_csv(dframe_name + ".csv", index=False)
    return final_dframe

def plot_advection_errors(advection_dframe, R, study=""): 

    resolutions = advection_dframe["DELTA_X"].unique()
    Evmax = []
    for resolution in resolutions:
        advection_data = advection_dframe[advection_dframe["DELTA_X"] == resolution]
        Evmax.append(advection_data["E_VOL_ALPHA"].max())
        plt.plot(advection_data["TIME"], advection_data["E_VOL_ALPHA"], 
                 label="%d cells / radius" % (R * (1 / resolution)))

    # First and last h values for convergence-order computation.
    h_01 = [advection_dframe["DELTA_X"].iloc[0],advection_dframe["DELTA_X"].iloc[-1]]

    title = "%s $Ev$" % study
    plt.title(title)
    plt.xlabel("time in seconds")
    plt.ylabel("$E_v$")
    plt.savefig("%s-volume-conservation-evolution.pdf" % study.replace(" ", "") , bbox_inches='tight')
    plt.show()

    plt.title(title)
    plt.plot(resolutions, Evmax, 'x-')

    Ev_error2nd_01 = [Evmax[0], Evmax[0]*(h_01[1]/h_01[0])**2]
    Ev_error1st_01 = [Evmax[0], Evmax[0]*(h_01[1]/h_01[0])]

    plt.plot(h_01,Ev_error2nd_01,"k--",label="second-order")
    plt.plot(h_01,Ev_error1st_01,"r:",label="first-order")

    plt.semilogy()
    plt.ylabel("$\max(Ev)$")
    plt.xlabel("$h$")
    plt.savefig("%s-volume-conservation-convergence.pdf" % study.replace(" ", "") , bbox_inches='tight')
    plt.show()

    Eg = []
    for resolution in resolutions:
        advection_data = advection_dframe[advection_dframe["DELTA_X"] == resolution]
        Eg.append(advection_data["E_GEOM_ALPHA"].iloc[-1])
    
    plt.title("%s Eg" % study)
    plt.ylabel("$E_g$")
    plt.xlabel("$h$")
    plt.plot(resolutions, Eg, label="resolution 1 / %d" % (1 / resolution))

    Eg_error2nd_01 = [Eg[0], Eg[0]*(h_01[1]/h_01[0])**2]
    Eg_error1st_01 = [Eg[0], Eg[0]*(h_01[1]/h_01[0])]

    plt.plot(h_01,Eg_error2nd_01,"k--",label="second-order")
    plt.plot(h_01,Eg_error1st_01,"r:",label="first-order")

    plt.show()


