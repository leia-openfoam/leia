import os
import os.path
from subprocess import run
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import PyFoam
import PyFoam.RunDictionary
import PyFoam.RunDictionary.SolutionFile
import PyFoam.RunDictionary.ParsedParameterFile
from PyFoam.Basics.DataStructures import Vector

import importlib.util
def import_module(path):
    name = os.path.basename(path).removesuffix('.py')
    # Create a spec object based on the module path
    spec = importlib.util.spec_from_file_location(name, path)
    # Import the module using the spec
    mod = importlib.util.module_from_spec(spec)
    # Load the module into the namespace
    spec.loader.exec_module(mod)
    return mod

module_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'postprocess_module.py')
mod = import_module(module_path)


def getInternalField(filename):
    
    file = PyFoam.RunDictionary.ParsedParameterFile.ParsedParameterFile(filename)
    field = file.getValueDict()["internalField"]
    return field.value()

def main():
    if not mod.isOpenFoamSourced:
        print("Source OpenFOAM before calling.")
        sys.exit(1)

    timedirs =  sorted(  
                        set(mod.getTimedirs()) 
                        & set(["0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8"])
    )

    for timedir in timedirs:
        # plot(x,Ux,psi,psi0,magGradPsi,timedir)
        plot(timedir)

# def plot(x,Ux,psi,psi0,magGradPsi,timedir):

def plot(timedir):

    U = getInternalField(os.path.join(timedir, "U"))
    Ux = list(map(lambda v: v[0], U))
    psi0 = getInternalField(os.path.join("0", "psi"))
    psi = getInternalField(os.path.join(timedir, "psi"))
    gradPsi = getInternalField(os.path.join(timedir, "gradPsi"))
    magGradPsi = list(map(abs, gradPsi))

    xStart, xEnd = 0, 1
    L = xEnd - xStart
    N = len(U)
    dx = L/N
    x = np.linspace(dx/2, L - dx/2, N)

    plt.figure(figsize=[8.8, 5.6])
    plt.plot(x, Ux, label=fr"$U=\sin(2 \pi x)$")
    psi_Line = plt.plot(x, psi, label=f"$\psi(t={timedir})$")[0]
    plt.plot(x, psi0, label=r"$\psi(t=0)$", color=psi_Line.get_color(), linestyle=':')
    plt.plot(x, magGradPsi, label=fr"$|\nabla \psi(t={timedir})|$")
    plt.xticks(np.linspace(0,1,11))
    plt.ylim((-2,12))
    legendtitle = f"source: {mod.getSourceType()}"
    if mod.getSourceType() == 'beta':
        alpha, beta = mod.getAlphaBeta()
        legendtitle += "\n" + rf'$\alpha={alpha}$, $\beta={beta}$'
        etaCutoff = mod.getEtaCutoff()
        legendtitle += "\n" + f'etaCutoff={etaCutoff}'
        psiCutoff = mod.getPsiCutoff()
        legendtitle += "\n" + f'psiCutoff={psiCutoff}'
    plt.legend(title=legendtitle)
    plt.grid()
    title = "1D problem with the compressible Level-Set Equation\n"\
        + r"$ \partial_t \psi + \nabla \cdot (U \psi) - (\nabla \cdot U) \psi = S(\psi, \nabla \psi; \alpha, \beta) $"\
        + "\nFields inside the domain"
    plt.title(title)
    plt.xlabel("position")
    plotfile = os.path.join(os.getcwd(), f"{os.path.basename(os.getcwd())}_fieldAt{timedir}Time.png")
    print(plotfile)
    # plt.show()
    plt.savefig(plotfile, bbox_inches='tight')
    

if __name__ == '__main__':
    main()
