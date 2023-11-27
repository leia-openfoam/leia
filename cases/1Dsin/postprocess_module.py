import os
from subprocess import run

def isOpenFoamSourced() -> bool:
    return "WM_PROJECT_DIR" in os.environ

def getSourceType() -> str:
    return run("foamDictionary -entry levelSet/sdplsSource/type -value system/fvSolution", 
        shell=True,
        capture_output=True,
        encoding='utf-8'
        ).stdout.removesuffix('\n')

def getAlphaBeta():
    if getSourceType() != 'beta':
        return (0,0)
    
    alpha = run("foamDictionary -entry levelSet/sdplsSource/alpha -value system/fvSolution", 
        shell=True,
        capture_output=True,
        encoding='utf-8'
        ).stdout.removesuffix('\n')
    alpha = float(alpha)
    beta = run("foamDictionary -entry levelSet/sdplsSource/beta -value system/fvSolution", 
        shell=True,
        capture_output=True,
        encoding='utf-8'
        ).stdout.removesuffix('\n')
    beta = float(beta)
    return (alpha, beta)

def getEtaCutoff() -> str:
    return run("foamDictionary -entry levelSet/sdplsSource/etaCutoff -value system/fvSolution", 
        shell=True,
        capture_output=True,
        encoding='utf-8'
        ).stdout.removesuffix('\n')

def getPsiCutoff() -> str:
    return run("foamDictionary -entry levelSet/sdplsSource/psiCutoff -value system/fvSolution", 
        shell=True,
        capture_output=True,
        encoding='utf-8'
        ).stdout.removesuffix('\n')

def getTimedirs(caseroot='.'):
    return run(f"foamListTimes -withZero -case {caseroot}",
    # return run(f"foamListTimes -case {caseroot}", 
            shell=True, capture_output=True, encoding='utf-8'
            ).stdout.splitlines() 