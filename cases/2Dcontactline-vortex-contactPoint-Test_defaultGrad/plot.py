#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def main():
    plot()

def reference(filename):
    df = pd.read_csv(filename)
    t = df.iloc[:,0]
    y = df.iloc[:,1]
    return (t,y)

def plot():
    df = pd.read_csv("contactPoint.csv", header=0)
    t = df.iloc[:,0]
    x = df.iloc[:,1]
    gamma = df.iloc[:,2]
    kappa = df.iloc[:,3]
    magGradPsi = df.iloc[:,4]
    
#    df2 = pd.read_csv("gradPsiError.csv", header=0)
#    t2 = df2.loc[:,"TIME"]
#    magGradPsi = df2.loc[:,"NARROW_MAX_MAG_GRAD_PSI"] 
    
    fig, axs = plt.subplots(2, 2, figsize=[13.2, 8.4])
    axs[0,0].plot(*reference('constant/reference/position.csv'), 'k--', label='reference')
    axs[0,0].plot(t,x, 'x', label='simulation')
    axs[0,0].set_xlabel("time")
    axs[0,0].set_ylabel("position")
    axs[0,0].legend()

    axs[0,1].plot(*reference('constant/reference/contactAngle.csv'), 'k--', label='reference')
    axs[0,1].plot(t,gamma, 'x', label='simulation')
    axs[0,1].set_xlabel("time")
    axs[0,1].set_ylabel("contactAngle")
    axs[0,1].legend()

    axs[1,0].plot(*reference('constant/reference/curvature.csv'), 'k--', label='reference')
    axs[1,0].plot(t,kappa, 'x', label='simulation')
    axs[1,0].set_xlabel("time")
    axs[1,0].set_ylabel("curvature")
    axs[1,0].legend()

    axs[1,1].plot(t,magGradPsi, 'x', label='simulation')
    axs[1,1].set_xlabel("time")
    axs[1,1].set_ylabel(r"|grad $\psi$|")
    axs[1,1].legend()
    plt.show()
#    plt.savefig("plot.png")
    

if __name__ == '__main__':
    main()
