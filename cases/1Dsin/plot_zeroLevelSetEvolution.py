import os
import os.path
import sys
from subprocess import run
import types
import matplotlib.pyplot as plt
import numpy as np
from PyFoam.RunDictionary.SolutionFile import SolutionFile
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


def writeCellCentres():
    run("writeMeshObj -time 0 > /dev/null", shell=True)

def readCellCentres():
    cc = []
    with open("meshCellCentres_0.obj") as f:
        for line in f:
            line = line.removeprefix('v ').removesuffix('\n')
            cc.append(Vector(*[float(fl_) for fl_ in line.split()]))
    return cc

def getCellCentres():
    if not os.path.isfile("meshCellCentres_0.obj"):
        writeCellCentres()
    return readCellCentres()


def getBetaSourcetermBounds():
    alpha, beta = mod.getAlphaBeta()
    return (beta-alpha, beta+alpha)


def plot(t,y):
    plt.figure(figsize=[8.8, 5.6])
    plt.plot(t, y, label=r"$| \nabla \psi |$")
    legendtitle = f"source: {mod.getSourceType()}"
    if mod.getSourceType() == 'beta':
        bounds = getBetaSourcetermBounds()
        plt.plot((t[0],t[-1]),bounds[1]*np.array([1,1]), 'k--')
        plt.plot((t[0],t[-1]),bounds[0]*np.array([1,1]), 'k--', label='bounds')

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
        + "\nZero level-set evolution"
    plt.title(title)
    plt.xlabel("time")
    plt.ylabel(r"$| \nabla \psi |$")
    plotfile = os.path.join(os.getcwd(), f"{os.path.basename(os.getcwd())}_zeroLevelSetEvolution.png")
    print(plotfile)
    # plt.show()
    plt.savefig(plotfile, bbox_inches='tight')

def main():
    if not mod.isOpenFoamSourced:
        print("Source OpenFOAM before calling.")
        sys.exit(1)

    timedirs = mod.getTimedirs()
    fieldname = 'magGradPsi'
    magGradPsi = []

    dc = Vector(0.5, 0, 0) # domain center
    centres = getCellCentres()
    dcidx = centres.index(dc)

    mgp = types.SimpleNamespace() # magGradPsi

    for tdir in timedirs:
        mgp.file = SolutionFile(tdir, fieldname)
        mgp.parsed = mgp.file.getContent() # -> return ParsedParameterFile
        mgp.dict = mgp.parsed.getValueDict()
        mgp.field = mgp.dict['internalField']
        magGradPsi.append(mgp.field[dcidx])

    times = np.asarray(timedirs, float)
    bounds = getBetaSourcetermBounds()
    plot(times, magGradPsi)

    

if __name__ == '__main__':
    main()
