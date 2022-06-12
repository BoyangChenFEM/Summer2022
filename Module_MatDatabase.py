# Database of materials for Abaqus models
# inputs: 
# - current model (mdb.Model object), 
# - name of material (string),
# - user-defined material parameters (optional)
# output: updated current model with material defined

from abaqus import *
from abaqusConstants import *

import material

def formCompositeMaterial(currentModel, materialChoice, eta_vr=1.e-5):
# eta_vr - viscous parameter for damage regularization (unit: time), default=1.e-5
    if materialChoice == 'IM7-8552':
        currentModel.Material(name=materialChoice)
        currentModel.materials[materialChoice].Elastic(table=((161000.0, 11000.0, 11000.0, 0.32, 0.32, 0.45,
                                                        5170.0, 5170.0, 3980.0),), type=ENGINEERING_CONSTANTS)
        currentModel.materials[materialChoice].HashinDamageInitiation(table=((2800.0, 1700.0, 60.0, 125.0, 90.0, 90.0),))
        currentModel.materials[materialChoice].hashinDamageInitiation.DamageEvolution(table=((50.0, 50.0, 0.22, 0.72),),
            type=ENERGY)
        currentModel.materials[materialChoice].hashinDamageInitiation.DamageStabilization(fiberCompressiveCoeff=eta_vr,
            fiberTensileCoeff=eta_vr, matrixCompressiveCoeff=eta_vr, matrixTensileCoeff=eta_vr)
# end formMaterial function    

# can be expanded to form other types of materials, difference would be in the optional user-defined parameters