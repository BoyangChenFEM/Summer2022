"""
Authors: Karthik Venkatesan (Aug 2019), Boyang Chen (June 2022)
Script for extracting images from odb files
Description: 
    imports the CompositesOpenHole module;
    extract damage pattern images and store them in a subfolder;
Inputs:
    damageLabels (list or tuple of strings): Abaqus damage variables to be plotted on 
    Plies (list or tuple of + integers): ply indices in Abaqus composite layup for plotting
    Nfiles (+ integer): no. of images to be extracted per damage variable per ply per job
Outputs:
    subfolder "damagePatterns" with all the extracted images of damage variable plots
"""

# custom imports
import os
import shutil
import time
import numpy as np
from Module_CompositesOpenHole import extractor_DamagePatterns

# ===================================================================================
# Set some parameters (see docstring of extractor_DamagePatterns for details):
damageLabels = ('DAMAGEFT',) # fibre tensile damage variable
Plies = (4,)
Nfiles = 10
# ===================================================================================

# set a reference time to calculate analysis time afterwards
tick = time.time()

# get the current working directory
main_directory = os.getcwd()

# create new subfolder to store images of damagePatterns
image_folder = main_directory + r'\damagePatterns\\'
if not os.path.exists(image_folder):os.makedirs(image_folder)

# Go to abaqus-models subfolder for data processing
job_folder = main_directory + r'\abaqus-models\\'
os.chdir(job_folder)

print "Start reading ODB files at time: ", time.time() - tick

# Loop over all ODB files in the job folder
for jobfile in os.listdir(job_folder): 
    if not jobfile.endswith(".odb"): continue
    
    # record time before starting the job
    toc = time.time()
    
    # form modelName
    modelName = jobfile.strip(".odb")

    # run the data extractor, print the requested images, get the list of image file names
    imageList=extractor_DamagePatterns(modelName, damageLabels, Plies, Nfiles)
    
    # move extracted images to image folder
    for imageName in imageList: shutil.move(job_folder+imageName, image_folder+'\\'+imageName)
    
    # record duration of the extraction
    print "Finished extraction: ", modelName, "; Time taken: ", time.time() - toc


# change back to main directory
os.chdir(main_directory)

print "Full script ended at time: ", time.time() - tick