"""
Authors: T.H.E. Gulikers (Dec 2018), Boyang Chen (June 2022)
Script for extracting strain and stress data from odb job files
Description: 
    imports the CompositesOpenHole module's data_extractor function to 
    extracts the strains & stresses from all odb files in the abaqus-models folder 
    and print the extracted data of each odb file to a txt file named after the odb.
NOTE: 
    eps12 is the tensor strain which is half of gamma12 (the engineering shear strain). 
    the data_extractor outputs (eps1, eps2, gamma12, sig1, sig2, sig12)
"""

# custom imports
import os
import shutil
import time
import numpy as np
from Module_CompositesOpenHole import extractor_StrainStress

# ===================================================================================
# Set some parameters:

# ===================================================================================

# set a reference time to calculate analysis time afterwards
tick = time.time()

# get the current working directory
main_directory = os.getcwd()

# Run models and generate job files in a sub directory
data_folder = main_directory + r'\StrainStressData\\'
if not os.path.exists(data_folder): os.makedirs(data_folder)
os.chdir(data_folder)

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
    
    # Open a text file to write the extracted stress-strain data
    reportFileName = modelName+'.txt'
    reportFile = open(reportFileName, 'w')
    reportTitle = 'Eps1, Eps2, 2*Eps12, Sigma_1, Sigma_2, Sigma_12'
    reportFile.write(reportTitle+'\n')
    dataArray = extractor_StrainStress(modelName)
    np.savetxt(reportFile, np.array(dataArray), delimiter = ',')
    reportFile.close()
    
    # move extracted data to data folder
    shutil.move(job_folder+reportFileName, data_folder+'\\'+reportFileName)
    
    # record duration of the extraction
    print "Finished extraction: ", modelName, "; Time taken: ", time.time() - toc

# change back to main directory
os.chdir(main_directory)

print "Full script ended at time: ", time.time() - tick