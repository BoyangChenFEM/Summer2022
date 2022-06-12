"""
Authors: T.H.E. Gulikers (Dec 2018), Boyang Chen (June 2022)
Main Script for Design of Experiments (DoE) loop
Description: 
    imports the CompositesOpenHole module;
    perform DoE loop on the Open Hole model_generator;
    DoE loop gives combinations of (eps1, eps2, eps12) as strain loadings, 
    creates a subfolder "abaqus-models" and stores inputfiles (and jobfiles if
    runJob parameter is set to True).
NOTE: 
    eps12 is the tensor strain which is half of gamma12 (engineering shear strain). 
    The DoE loop input (eps1, eps2, eps12) to the model_generator
"""

# custom imports
import os
import shutil
import time
import numpy as np
import itertools
from Module_CompositesOpenHole import model_generator

# ===================================================================================
# Set some parameters:
runJob = False   # True to directly submit the job; False for only writing input file
ncpus  = 1      # abaqus student edition only supports single cpu job

# Design of Experiments on strain loadings
# loading ranges
limitsEps1 = limitsEps2 = limitsEps12 = [-1.0e-2, 1.0e-2]
# no. of grid points in each loading range
ngrid = 5
# generate the list of strains in each direction
Eps1_DoE = np.linspace(limitsEps1[0], limitsEps1[1], ngrid)
Eps2_DoE = np.linspace(limitsEps2[0], limitsEps2[1], ngrid)
Eps12_DoE = np.linspace(limitsEps12[0], limitsEps12[1], ngrid)
# ===================================================================================

# set a reference time to calculate analysis time afterwards
tick = time.time()

# get the current working directory
main_directory = os.getcwd()

# Run models and generate job files in a sub directory
new_directory = main_directory + r'\abaqus-models\\'
if not os.path.exists(new_directory): os.makedirs(new_directory)
os.chdir(new_directory)

# See what jobs are already there
currentModels = [file for file in os.listdir(new_directory) if file.endswith(".odb")]

print "Start DoE loop at time: ", time.time() - tick

for Eps1, Eps2, Eps12 in itertools.product(Eps1_DoE, Eps2_DoE, Eps12_DoE):

    # skip directly to next iteration if none of the strain is at limit
    if all([Eps1 not in limitsEps1, Eps2 not in limitsEps2, Eps12 not in limitsEps12]): continue
    
    # form modelName based on micro-strain values
    modelName = 'Ea' +  str(int(Eps1*1e6)) + '_Eb' + str(int(Eps2*1e6)) + '_Es' + str(int(Eps12*1e6))
    
    # if this model has already been run, then skip to the next one
    if modelName+'.odb' in currentModels: continue
    
    # record time before starting the job
    toc = time.time()
    
    # run the model generator (and the job if runJob=True)
    model_generator(modelName, Eps1, Eps2, Eps12, runJob, ncpus)
    
    # record duration of the job
    print "Finished job: ", modelName, "; Time taken: ", time.time() - toc

# change back to main directory
os.chdir(main_directory)

print "Full script ended at time: ", time.time() - tick