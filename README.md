# Summer2022
Training ANNs with FEM data on Open Hole Composite plate - for summer school 2022 in Delft
Based on the former Master theses work of Tom Gulikers (https://github.com/tgulikers/ABAQUS_ANN_constitutive_model) and Karthik Ventakesan (https://github.com/k-venkatesan/composite-damage-prediction). For more in-depth details, see their repos and theses.

## Pre-requisites
  Abaqus (full or student edition) installed; Student edition only runs on Windows, can be downloaded from https://edu.3ds.com/en/software/abaqus-learning-edition
  Google Drive and Colab notebook (or locally installed Tensorflow Keras and jupyter notebook)

## General description
  A square composite plate with an open hole is a recurring substructure in bigger composite structures. The FE modelling of this substructure requires fine meshes, particularly for damage analysis. The objective of this project is to see if we can build a Neural Network surrogate model for this substructure under uniform planar strain loadings (epsilon 1, epsilon 2 and gamma 12). Abaqus FE models of this substructure (quasi-isotropic layup [45/90/-45/0]s) under different combinations of strain loading are built and run using Python scripting. Hashin criterion with cohesive damage propagation is used to model the composite damage. No delamination is considered for simplicity. After the jobs are finished, data (strain/stress) are extracted from job files using Python scripts and stored in subfolders. A Google Colab notebook is then written to read the data and train a dense NN using Tensorflow Keras to quickly estimate the average stresses given the strain loadings.

## Description of files

*Module_CompsitesOpenHole.py*
  A module written in Abaqus python scripts which defines the following functions:
  - model_generator: generate the FE model of the substructure (geometry, material, layup, mesh, element type, step settings) given the strain loadings. It run the models if the runJob argument is passed in as True
  - extractor_StrainStress: extract strain-stress data history
  - extractor_DamagePatterns: extract a user-defined number of damage patterns of a request ply for a requested damage variable 
  The functions in this module needs to be called by a script (see below). Currently, only the strain loadings are considered variable. Other parameters such as the ply angles may also be added to the list of arguments for sampling.

*Module_MatDatabse.py*:
  A module to store the materials, to be used by other modules such as the Module_CompositesOpenHole

*Script-generateModels.py*:
  This is the Abaqus python script to generate the input files of the substructure. It will directly run all the jobs if the runJob parameter in the script is set to True. Its function is to sample the strain loadings through a (static) Design of Experiment, then calls the Module_CompsitesOpenHole module_generator function to generate/run the models.
  
*Script-extractStrainStress.py*:
  imports the CompositesOpenHole module's extractor_StrainStress function to extracts the strains & stresses from all odb files in the abaqus-models folder and print the extracted data of each odb file to a txt file named after the odb.
NOTE: eps12 is the tensor strain which is half of gamma12 (the engineering shear strain). the extractor_StrainStress outputs (eps1, eps2, gamma12, sig1, sig2, sig12)

*Script-extractDamagePatterns.py*:
  imports the CompositesOpenHole module; extract damage pattern images and store them in a subfolder;
  Inputs:
    damageLabels (list or tuple of strings): Abaqus damage variables to be plotted on 
    Plies (list or tuple of + integers): ply indices in Abaqus composite layup for plotting
    Nfiles (+ integer): no. of images to be extracted per damage variable per ply per job
  Outputs:
    subfolder "damagePatterns" with all the extracted images of damage variable plots

*abaqus-models.zip*: input files created by running Script-generateModels.py

*StrainStressData.zip*: strain-stress data (txt files) of all the jobs from running the above input files

*damagePatterns.zip*: damage images of 0-dgree ply (ply 4 of the layup definition) fiber tensile damage ('DAMAGEFT' variable in Abaqus)

*ANN-StrainStress.ipynb*: jupyter notebook (run in Colab) to read StrainStressData from Google Drive, then train a Keras Tensorflow NN model. Before running this notebook, first download the StrainStressData.zip and unzip it in your Colab Notebooks folder. The notebook will ask to mount your Google Drive, which will require your authorization (remember to allow pop-ups from Colabo website).
