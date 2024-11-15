#!/usr/bin/env python

from optparse import OptionParser
import sys
from subprocess import call
import os

usage = """A wrapper for pyFoamRunParameterVariation.py that generates the
directory structure for a parameter study. 

Meshes are not generated and preprocessing is not done. 
Used to prepare for execution on a cluster.

Usage: ./create-study.py -c templateCase -p paramFile -s studyName"""

study_name=""

parser = OptionParser(usage=usage)

parser.add_option("-c", "--case", dest="casedir",
                  help="Template case directory.", 
                  metavar="CASEDIR")

parser.add_option("-p", "--parameter-file", dest="paramfile", 
                  help="PyFoam parameter file used by pyFoamRunParameterVariation.py.", 
                  metavar="PARAMFILE")

parser.add_option("-s", "--study-name", dest="studyname", 
                  help="Name of the parameter study.", 
                  metavar="STUDYNAME")
if __name__ == "__main__":

    (options, args) = parser.parse_args()

    if ((options.casedir == None) or  
        (options.paramfile == None) or 
        (options.studyname == None)): 
        print ("Error: case or parameter option not used. Use --help option for more information.") 
        sys.exit(1)

    (options, args) = parser.parse_args()

    numberOfVariations = [0,1,2,3,4,5,6,7,8,9,10,11] 
 
    study_name = options.studyname
    
    call(["pyFoamRunParameterVariation.py", 
          "--no-execute-solver", 
          "--no-server-process", 
          "--no-mesh-create", 
          "--no-case-setup", 
          "--cloned-case-prefix=%s" % options.studyname, 
          "--every-variant-one-case-execution", 
          "--create-database", 
          options.casedir, 
          options.paramfile])
