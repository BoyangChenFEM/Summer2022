"""
Authors: 
    T.H.E. Gulikers (Dec 2018), Karthik Venkatesan (Aug 2019), Boyang Chen* (June 2022)
    @ Aerospace Engineering, TU Delft, b.chen-2@tudelft.nl
Description: 
    Model generator & data_extractor for Composites Open Hole plate 
    under biaxial&shear loading with Periodic boundary conditions
Currently fixed parameters:
    - Geometry: hole diameter = 6mm; width(X)=length(Y)=5*hole diameter
    - Material: composite, IM7-8552, viscous damage parameter eta_vr= 1.e-5
    - ply thickness: 0.125 mm
    - Layup: (45, 90, -45, 0)_sym
    - Mesh: 10 nodes on each edge, 40 nodes around hole, radial mesh
    - Element type: S4 (conventional quad shell, no reduced integration)
    - Step settings: 
      - initialInc=maxInc=0.1; minInc=0.01*eta_vr; num_inc=1/eta_vr
      - I_a (max number of increment size cutbacks)=15 (abaqus default=5)
      - nlgeom = OFF (i.e., geometrically linear analysis)
      - extrapolation = NONE (i.e., full Newton-Raphson solver)
      - no general damping/viscous regularization in step
Functions:
    - model_generator: build model, write input file or run the job
    - data_extractor: extract some selected data from the job files
"""

# default abaqus imports
from abaqus import *
from abaqusConstants import *
import __main__

#================================================================================
#------------------------ START PARAMETERS BLOCK --------------------------------
#================================================================================
# geometrical dimensions of the problem based on Green et al., Comp. Part A, 2007
diameter = 6. # mm
length = width = 5*diameter # square section from (0,0) to (width, length)
position_ellipse = (width/2., length/2.)     # [mm] width, length indicating ellipse centre
dimension_ellipse = (diameter/2., diameter/2.)    # [mm] half of the hole size axis along [x, y] axis.

# material name (must be in the Module_MatDatabase.py)
materialChoice = 'IM7-8552'

# composite layup parameters
t_ply = 0.125 # mm
theta = 45. # degrees
layup = (theta, 90., -theta, 0.)
makeSymLayup = True            # creates symmetric layup if True
depth = t_ply*len(layup)
if makeSymLayup: depth *= 2

# mesh densities
seedsOnEdge = 10            # amount of seeds per outer edge. 
seedsOnHole = 4*seedsOnEdge # The elliptic hole edge has 4x seedsOnEdge to generate radial mesh
maxNodeID = 10*seedsOnHole**2 # ref nodes' starting index, must be new nodes, hence a large starting index.

# time incrementation (assumed total time period = 1)
max_dt  = 0.1  # initial & max time increment
eta_vr  = 1.e-5 # viscous regularization parameter for damage 
min_dt  = 0.01*eta_vr # min. time increment, must be smaller than eta_vr
num_inc = int(1./eta_vr) # estimate of the total inc. num.
I_a = 15.0 # no. of increment size cutbacks allowed, used on step's solver control settings

#================================================================================
# ----------------------- END PARAMETERS BLOCK ----------------------------------
#================================================================================


def model_generator(modelName, epsilon_1, epsilon_2, epsilon_12, runJob=False, ncpus=1):
    """
    Inputs:
    - modelName  : string for name of the model
    - epsilon_1  : normal strain loading along X
    - epsilon_2  : normal strain loading along Y
    - epsilon_12 : shear  strain loading (note: epsilon_12=1/2*gamma_12, where gamma_12 is the engineering shear strain)
    - runJob     : default is .False. to generate the inputfile only; Set it to .True. to run the job directly
    - ncpus      : default is 1, set higher for parallel computing (note: student version does not allow parallel)
    (Note: other parameters can be added as inputs if needed by the Design of Experiment, e.g.: layup, hole diameter)
    outputs (in the current work directory, jobName = modelName):
    - if runJob is False, then only the inputfile
    - if runJob is True, then the full set of abaqus results
    """
    # default abaqus imports
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    # import visualization
    # import xyPlot
    # import displayGroupOdbToolset as dgo
    # import connectorBehavior
    
    # custom imports
    import numpy as np
    from Module_MatDatabase import formCompositeMaterial


    # create the model
    mdb.Model(name=modelName)
    currentModel = mdb.models[modelName]

    # create the sketch
    currentSketch = currentModel.ConstrainedSketch(name='2D_sketch', sheetSize=float(np.max((width, length, depth))))
    # draw rectangle
    currentSketch.rectangle((0., 0.), (width, length))
    # draw hole
    currentSketch.EllipseByCenterPerimeter(center= (width/2., length/2.), axisPoint1 =(width/2., length/2.-dimension_ellipse[1]),
                                                        axisPoint2 = (width/2.-dimension_ellipse[0], length/2.))
                                                        
    # create part from sketch
    currentPart = currentModel.Part(name='part_3D', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    currentPart.BaseShell(sketch=currentSketch)
    del currentModel.sketches['2D_sketch']

    # coordinate system
    CSYS_base = currentPart.DatumCsysByThreePoints(coordSysType= CARTESIAN, line1=(1.0, 0.0, 0.0),
                                                line2=(0.0, 1.0, 0.0), name='csys_base', origin=(0.0, 0.0, 0.0))

    # Create materials using function from Module_MatDatabase
    formCompositeMaterial(currentModel, materialChoice, eta_vr)
    
    #region create composite section
    # Note: the outer brackets are needed, otherwise the region object under CompositePly does not accept it.
    sectionRegion = (currentPart.faces.findAt((width*0.99, length*0.99, 0.),),) # location on the plate

    # composite layup object initiation
    currentSection = currentPart.CompositeLayup(description='', elementType=SHELL, name='QI_section',
                offsetType=MIDDLE_SURFACE, symmetric=makeSymLayup, thicknessAssignment=FROM_SECTION)
    currentSection.Section(integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, thicknessType=UNIFORM,
                           useDensity=OFF)
    currentSection.ReferenceOrientation(additionalRotationType=ROTATION_NONE, angle=0.0, axis=AXIS_3, fieldName='',
        localCsys=None, orientationType=GLOBAL)

    # create plies
    for i, angle_ply in enumerate(layup, start=1):
        currentSection.CompositePly(thickness=t_ply, angle=angle_ply, axis=AXIS_3, material= materialChoice,
                                    orientationType=CSYS, orientation=currentPart.datums[CSYS_base.id],
                                    plyName='Ply-'+str(i), region=sectionRegion, thicknessType=SPECIFY_THICKNESS)

    #region create assembly and define edges
    currentAssembly = currentModel.rootAssembly
    currentAssembly.Instance(name='Part Instance', part=currentPart, dependent=ON)
    currentInstance = currentAssembly.instances['Part Instance']

    # bottom, left, top, right, edge objects
    edgefind = currentInstance.edges.findAt
    edges_b = edgefind((( width/2., 0., 0.),),)
    edges_l = edgefind(((0., length/2., 0.),),)
    edges_t = edgefind(((width/2., length, 0.),),)
    edges_r = edgefind(((width, length/2., 0.),),)
    edges_ellipse = edgefind(((width/2.-dimension_ellipse[0], length/2.,0.),),)

    # define sets from edges
    currentAssembly.Set(name='edge_bottom', edges=edges_b)
    currentAssembly.Set(name='edge_left',   edges=edges_l)
    currentAssembly.Set(name='edge_top',    edges=edges_t)
    currentAssembly.Set(name='edge_right',  edges=edges_r)
    currentAssembly.Set(name='edge_circle',  edges=edges_ellipse)

    # define mesh controls (plane stress element)
    currentMeshRegion = sectionRegion
    meshElemType = mesh.ElemType(elemCode=S4, elemLibrary=STANDARD)
    currentPart.setMeshControls(algorithm=MEDIAL_AXIS, elemShape=QUAD, regions=currentMeshRegion)
    currentPart.setElementType(regions=currentMeshRegion, elemTypes=(meshElemType,))

    # define mesh seeds
    currentPart.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=seedsOnEdge)
    currentPart.seedEdgeByNumber(constraint=FINER, edges=edges_ellipse, number=seedsOnHole)
    for edge in (edges_l, edges_r):  currentPart.seedEdgeByNumber(constraint=FINER, edges=edge, number=seedsOnEdge)
    for edge in (edges_b, edges_t):  currentPart.seedEdgeByNumber(constraint=FIXED, edges=edge, number=seedsOnEdge)

    # generate mesh
    currentPart.generateMesh()
    currentAssembly.regenerate()

    # define the analysis steps
    currentModel.StaticStep(name='ApplyLoad1', previous='Initial', description='Load is applied in this step',
                            initialInc=max_dt, maxInc=max_dt, maxNumInc=num_inc, minInc=min_dt, extrapolation=NONE)
                            
    # set the step's solver control parameters for the incrementation scheme
    # non-default is in symbol with value defined above
    currentModel.steps['ApplyLoad1'].control.setValues(allowPropagation=OFF, resetDefaultValues=OFF,
                                                       timeIncrementation=(4.0, 8.0, 9.0, 16.0, 10.0, 4.0,
                                                                           12.0, I_a, 6.0, 3.0, 50.0))

    # Field output requests
    currentModel.fieldOutputRequests.changeKey(fromName='F-Output-1', toName='FieldOutputGlobal')
    currentModel.fieldOutputRequests['FieldOutputGlobal'].setValues(variables=('S', 'E', 'U', 'RF', 'CF'))
    currentModel.FieldOutputRequest(createStepName='ApplyLoad1', layupLocationMethod=SPECIFIED, layupNames=(
        'Part Instance.QI_section', ), name='FieldOutputComposite', outputAtPlyBottom=False,
        outputAtPlyMid=True, outputAtPlyTop=False, rebar=EXCLUDE, variables=('S', 'E',
        'DAMAGEFT','DAMAGEFC', 'DAMAGEMT', 'DAMAGEMC', 'DAMAGESHR', 'DMICRT'))

    # History output requests
    currentModel.historyOutputRequests.changeKey(fromName='H-Output-1', toName='HistoryOutputs')
    currentModel.historyOutputRequests['HistoryOutputs'].setValues(variables=PRESELECT)


    # =======================================================================================
    # ==================== START BOUNDARY CONDITIONS AND LOADS ==============================
    # =======================================================================================
    #
    # ---------------------------------------------------------------------------------------
    # Start region apply (periodic) boundary conditions
    # ---------------------------------------------------------------------------------------
    # find matching node pairs on model boundaries
    pairs1, pairs2 = [], []                           # matching nodes for PBC in DoF 1 and DoF 2 (global coordinate system)
    edgelabels = []

    for set_name in ['edge_left', 'edge_right']:     # loop over direction 1
        for i_node in range(len(currentAssembly.sets[set_name].nodes)):         # loop over nodes in set
            label = currentAssembly.sets[set_name].nodes[i_node].label      # save node number and coordinates in list
            coords = currentAssembly.sets[set_name].nodes[i_node].coordinates
            pairs1.append([label]+list(coords)) # list of lists, each sublist contains node label and coordiates

            # define set (required for equation constraint)
            currentAssembly.SetFromNodeLabels(name='Node_' + str(label), nodeLabels=(('Part Instance', (label,)),))
            edgelabels.append(label)

    for set_name in ['edge_bottom', 'edge_top']:     # loop over direction 2
        for i_node in range(len(currentAssembly.sets[set_name].nodes)):         # loop over nodes in set
            label = currentAssembly.sets[set_name].nodes[i_node].label      # save node number and coordinates in list
            coords = currentAssembly.sets[set_name].nodes[i_node].coordinates
            pairs2.append([label]+list(coords))

            # define set (required for equation constraint)
            currentAssembly.SetFromNodeLabels(name='Node_' + str(label), nodeLabels=(('Part Instance', (label,)),))
            edgelabels.append(label)

    # convert lists of pairs to numpy arrays
    # i is a list containing node lable and coordinates, it is converted to tuple
    # then all these tuples are grouped into a list and converted to np array
    # the nodelst_dtype give the names and types of 4 columns in the np array
    nodelst_dtype = [('label', 'i4'), ('x', 'f4'), ('y', 'f4'), ('z', 'f4')]  # data type specification needed for sorting
    pairs1, pairs2 = np.array([tuple(i) for i in pairs1], dtype=nodelst_dtype), np.array([tuple(i) for i in pairs2], dtype=nodelst_dtype)

    # sort pair-lists based on coordinates such that the first half of the list is edge A and the second half is edge B
    pairs1.sort(axis=0, order = ['x', 'y'])     # first x to divide the sides, then y to sort sequence on edge
    pairs2.sort(axis=0, order = ['y', 'x'])     # first y to divide the sides, then x to sort sequence on edge
    mp1, mp2 = len(pairs1)/2, len(pairs2)/2     # indices of mid-points, which is the division between edges

    # define reference points to be used for loading and define a set for each
    refNnormal = maxNodeID+1 # define them to be larger than the max node ID in the mesh to avoid conflict
    refNshear = maxNodeID+2
    RP_normal_node = currentPart.Node(coordinates=(width*1.1, length/2., 0.)).setValues(label=refNnormal)
    RP_shear_node = currentPart.Node(coordinates=(width*-0.1, length/2., 0.)).setValues(label=refNshear)
    currentAssembly.SetFromNodeLabels(name='RP_normal', nodeLabels=(('Part Instance', (refNnormal,)),))
    currentAssembly.SetFromNodeLabels(name='RP_shear', nodeLabels=(('Part Instance', (refNshear,)),))

    # equation constraints
    for i_pair in range(0, int(np.shape(pairs1)[0])/2):      # constraints on pairs in x-direction
        equationName = 'eq_'+str(pairs1[i_pair][0])+'_'+str(pairs1[i_pair+mp1][0])   # combination of node numbers

        # name, sequence of (Float, String, Int): (coefficient, Set name, DoF)
        # right edge U2 - left edge U2 = RP_shear U1
        # U2 of left corner nodes coupled to U2 of right corner nodes
        currentModel.Equation('S'+equationName, ((1, 'Node_' + str(pairs1[i_pair][0]), 2),
                (-1, 'Node_' + str(pairs1[i_pair+mp1][0]), 2), (1, 'RP_shear', 1)))   # shear, coupled to RP_shear DoF1
        
        # right edge U1 - left edge U1 = RP_normal U1
        # U1 of bottom corner nodes are excluded to avoid overconstraint conflict
        # they are constraint to the U1 of the top corner nodes (see below, first constraint in the next for loop)
        # while the U1 of the two top corner nodes will be coupled to each other here, hence it implies that
        # U1 of the bottom corner nodes will be automatically coupled to each other.
        if i_pair != 0:# and i_pair != int(np.shape(pairs1)[0]-2): 
            currentModel.Equation('N'+equationName, ((1, 'Node_' + str(pairs1[i_pair][0]), 1),
                (-1, 'Node_' + str(pairs1[i_pair+mp1][0]), 1), (1, 'RP_normal', 1)))  # normal, coupled to RP_normal DoF1

    for i_pair in range(0, int(np.shape(pairs2)[0])/2):      # constraints on pairs in y-direction
        equationName = 'eq_'+str(pairs2[i_pair][0])+'_'+str(pairs2[i_pair+mp2][0])   # combination of node numbers

        # name, sequence of (Float, String, Int): (coefficient, Set name, DoF)
        # top edge U1 - bottom edge U1 = RP_shear U2
        # U1 of top corner nodes coupled to U1 of bottom corner nodes
        currentModel.Equation('S'+equationName, ((1, 'Node_' + str(pairs2[i_pair][0]), 1),
          (-1, 'Node_' + str(pairs2[i_pair+mp2][0]), 1), (1, 'RP_shear', 2)))  # shear, coupled to RP_shear DoF2

        # top edge U2 - bottom edge U2 = RP_normal U2
        # U2 of left corner nodes are excluded to avoid overconstraint conflict
        # U2 of left corner nodes are coupled to U2 of right corner nodes in the first equation of the above for loop
        # while U2 of the right corner nodes are coupled to each other here. hence no need to specify the coupling of U2
        # of the left corner nodes.
        if i_pair != 0:# and i_pair != int(np.shape(pairs2)[0]-2): # exclude left pair to avoid overconstraint
            currentModel.Equation('N'+equationName, ((1, 'Node_' + str(pairs2[i_pair][0]), 2),
                    (-1, 'Node_' + str(pairs2[i_pair+mp2][0]), 2), (1, 'RP_normal', 2)))  # normal, coupled to RP_normal DoF2
    # ---------------------------------------------------------------------------------------
    # end region apply (periodic) boundary conditions
    # ---------------------------------------------------------------------------------------
    #
    # ---------------------------------------------------------------------------------------
    # start region apply fixed boundary conditions and loads
    # ---------------------------------------------------------------------------------------
    # Region object from regionToolset module is needed repeatedly in this section, extract it here for brevity
    Region = regionToolset.Region

    # fix out of plane displacement of the whole instance
    # use the findAt tool to find the whole region of the current instance
    # sidenote: attempted to use the sectionRegion object formed earlier as my Region, it led to exception error
    currentFaces = currentInstance.faces.findAt(((width*0.99, length*0.99, 0.),),)
    myRegion = Region(faces=currentFaces)
    currentModel.DisplacementBC('fixU3', createStepName='Initial',
                                region=myRegion, u3=0.)

    # find a node that is not on the edge and pin it such that the model doesn't fly away
    for i in range(1,maxNodeID):
        if i not in edgelabels:
            break
    currentModel.DisplacementBC('pin', createStepName='Initial', u1 = 0., u2 = 0.,
                                 region=Region(nodes=currentInstance.nodes.sequenceFromLabels((i,))))

    # create meshNodeObjects from reference nodes (required for Region function)
    meshNodeObjNormal = currentInstance.nodes.sequenceFromLabels((refNnormal,))
    meshNodeObjShear = currentInstance.nodes.sequenceFromLabels((refNshear,))

    # apply normal loads
    # u1 of RP_normal represents the normal load applied on pairs1 (left - right edges), hence it is epsilon_1*width
    # u2 of RP_normal represents the normal load applied on pairs2 (top - bottom edges), hence it is epsilon_2*length
    currentModel.DisplacementBC('load_RP_normal', createStepName='Initial', region=Region(nodes=meshNodeObjNormal), u1=0.,
                                u2=0.)
    currentModel.boundaryConditions['load_RP_normal'].setValuesInStep('ApplyLoad1', u1 = epsilon_1*width,
                                                                      u2 =epsilon_2*length)

    # apply shear load
    # u1 of RP_shear represents the shear load applied on pairs1 (left - right edges), hence it is epsilon_12*width
    # u2 of RP_shear represents the shear load applied on pairs2 (top - bottom edges), hence it is epsilon_12*length
    currentModel.DisplacementBC('load_RP_shear', createStepName='Initial', region=Region(nodes=meshNodeObjShear), u1=0.,
                                u2=0.)
    currentModel.boundaryConditions['load_RP_shear'].setValuesInStep('ApplyLoad1', u1 = epsilon_12*width,
                                                                     u2 = epsilon_12*length)
    # ---------------------------------------------------------------------------------------
    # endregion apply fixed boundary conditions and loads
    # ---------------------------------------------------------------------------------------
    #
    # =======================================================================================
    # ==================== END BOUNDARY CONDITIONS AND LOADS ================================
    # =======================================================================================


    # Create job (ncpus =1 by default, but is an optinal user-input)
    # NOTE: jobname is the same as modelName
    mdb.Job(name=modelName, model=modelName, description='Analysis of '+modelName,
            type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None,
            memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
            scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=ncpus,
            numGPUs=0, numDomains=ncpus)
    if runJob:
        mdb.jobs[modelName].submit()
        mdb.jobs[modelName].waitForCompletion()
    else:
        mdb.jobs[modelName].writeInput(consistencyChecking=OFF)
#End of function model_generator (add this comment line so that code folding ends here)



def extractor_StrainStress(modelName):
    """
    function to extract strain&stress data from odb files of Composites Open Hole models.
    It assumes that the job file exists in the current work directory.
    It uses some module parameters and is highly tailored to the model_generator function
    Input: modelName (string, = jobName),
    Output: 
        dataArray (list of lists), each row list is [eps1, eps2, gamma12, sig1, sig2, sig12]
        corresponding to a certain time frame of the analysis step. Hence, dataArray contains 
        the strain and stress data of all the time frames in the job, excluding the initial
        frame where all data are zero.
    """
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    
    # open .odb object
    odbName = modelName+'.odb'
    odbObject = session.openOdb(name=odbName) #session is imported from abaqus
    
    # label relevant node sets
    # setBottom   = odbObject.rootAssembly.nodeSets['EDGE_BOTTOM']
    # setLeft     = odbObject.rootAssembly.nodeSets['EDGE_LEFT']
    # setTop      = odbObject.rootAssembly.nodeSets['EDGE_TOP']
    # setRight    = odbObject.rootAssembly.nodeSets['EDGE_RIGHT']
    setNormal   = odbObject.rootAssembly.nodeSets['RP_NORMAL']
    setShear    = odbObject.rootAssembly.nodeSets['RP_SHEAR']
    
    # # find corner nodes
    # nodeBottom_L = [i for i in range(len(setBottom.nodes[0])) if setBottom.nodes[0][i].coordinates[0] == 0.][0]
    # nodeBottom_R = [i for i in range(len(setBottom.nodes[0])) if setBottom.nodes[0][i].coordinates[0] == width][0]
    # nodeLeft_B   = [i for i in range(len(setLeft.nodes[0])) if setLeft.nodes[0][i].coordinates[1] == 0.][0]
    # nodeLeft_T   = [i for i in range(len(setLeft.nodes[0])) if setLeft.nodes[0][i].coordinates[1] == length][0]
    
    dataArray = [] # initialize data array for printing in file
    
    step = odbObject.steps.values()[0] # list of step 1 frames
    
    for frame in step.frames:
        
        # store reaction-force and displacement field
        RFfield = frame.fieldOutputs['RF']
        Ufield  = frame.fieldOutputs['U']

        # extract lists of forces and displacements at edges
        RF_1, RF_2, RF_12 = RFfield.getSubset(region=setNormal).values[0].data[0],\
                            RFfield.getSubset(region=setNormal).values[0].data[1],\
                            [RFfield.getSubset(region=setShear).values[0].data[0],
                            RFfield.getSubset(region=setShear).values[0].data[1]]
        U_1,  U_2,  U_12  = Ufield.getSubset(region=setNormal).values[0].data[0], \
                            Ufield.getSubset(region=setNormal).values[0].data[1], \
                            [Ufield.getSubset(region=setShear).values[0].data[0],
                            Ufield.getSubset(region=setShear).values[0].data[1]]
        
        # compute current strains
        eps_1  = U_1/width
        eps_2  = U_2/length
        eps_12 = 0.5*(U_12[0]/width + U_12[1]/length)
        gamma_12 = 2*eps_12
        # NOTE: 1st DoF of RP_shear is for shearing the 1st pair of edges (Left-Right), hence the corresponding
        # shear strain should be U_12[0]/width; similarly, the 2nd DoF is for Top-Bottom pair, hence the 
        # corresponding strain should be U_12[1]/length. The two components shall be averaged to get eps12
        # or summed to get gamma_12 (engineering strain)

        # compute current stresses 
        # NOTE: use Hill-Mandel principle of energy balance to derive stresses
        # strain energy = 1/2*[eps]*[sig]*V = sum(1/2 * RF_i * U_i) 
        sig_1 = RF_1/length/depth
        sig_2  = RF_2/width/depth
        if abs(gamma_12) > 0.:#avoid division by 0
            sig_12 = (RF_12[0]*U_12[0]+RF_12[1]*U_12[1])/width/length/depth/gamma_12
        else:
            sig_12 = 0.
        
        #append to data array
        dataArray.append([eps_1, eps_2, gamma_12, sig_1, sig_2, sig_12])
    
    # return dataArray excluding the first all-zero row for the initial time frame
    return dataArray[1:]
#End of function extractor_StrainStress (add this comment line so that code folding ends here)



def extractor_DamagePatterns(modelName, damageLabels, Plies, Nfiles=10):
    """
    function to extract damage images from odb files of Composites Open Hole models.
    It assumes that the job file exists in the current work directory.
    It uses some module parameters and is highly tailored to the model_generator function
    Inputs: 
        modelName (string, = jobName)
        damageLabels (list or tuple): ('HSNFTCRT',) or any one(s) of the following:
            HSNFTCRT:   Fiber tensile initiation criterion.
            HSNFCCRT:   Fiber compressive initiation criterion.
            HSNMTCRT:   Matrix tensile initiation criterion.
            HSNMCCRT:   Matrix compressive initiation criterion.
            DAMAGEFT:   Fiber tensile damage variable.
            DAMAGEFC:   Fiber compressive damage variable.
            DAMAGEMT:   Matrix tensile damage variable.
            DAMAGEMC:   Matrix compressive damage variable.
            DAMAGESHR:  Shear damage variable.
        Plies (list or tuple): Abaqus ply indices to be printed (start indexing from 1):
            Since layup is [45,90,-45,0]_Sym, Plies = (4,) would request 0-degree ply
        Nfiles (positive int): number of files for output; default = 10; sampled across frames [1:] 
    Output: 
        - images (PNG) of the request damage labels of the requested plies. A uniform scale 
        factor of 1.0 is used to plot the deformed shape. Images are dimensions 256x256 pixels.
        The name of the file contains all the necessary information for postprocessing, e.g.:
        damage label, ply number, total microstrains (Ea,Eb,Es), and current time*1e6;
        - Filelist (list): a list of names (string) of the images printed
    """
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import numpy as np
    
    # ============= Set some parameters =======================
    # set interval Values to better plot damage variable range
    intervalValues = (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, )
    numIntervals = len(intervalValues)-1
    
    # set image dimensions for printing:
    ImageX = 256
    ImageY = 256
    # =========================================================

    # open .odb object
    myViewport = session.viewports['Viewport: 1']
    # myViewport = session.Viewport(name='Print damage pattern', origin=(10, 10),
    # width=150, height=100)
    
    odbName = modelName+'.odb'
    odbObject = session.openOdb(name=odbName) #session is imported from abaqus
    
    myViewport.setValues(displayedObject=odbObject)
    
    myViewport.view.setValues(session.views['Front'])

    # Contour Options
    # myViewport.odbDisplay.commonOptions.setValues(
        # renderStyle=FILLED, visibleEdges=FEATURE)
    myViewport.odbDisplay.commonOptions.setValues(
        renderStyle=SHADED, visibleEdges=NONE, deformationScaling=UNIFORM, 
        uniformScaleFactor=1.0)
        
    myViewport.odbDisplay.contourOptions.setValues(
        intervalType=USER_DEFINED, intervalValues=intervalValues)
        
    myViewport.odbDisplay.contourOptions.setValues(numIntervals=numIntervals)

    myViewport.viewportAnnotationOptions.setValues(triad=OFF,
            legend=OFF, title=OFF, state=OFF, annotations=OFF, compass=OFF)
    
    # set print options
    session.pngOptions.setValues(imageSize=(ImageX, ImageY))
    session.printOptions.setValues(vpDecorations=OFF)
    
    # loop over step 1 frames to extract image data at intervals of Nfiles
    step = odbObject.steps.values()[0]
    
    numFrames = len(step.frames)
    
    StartFrame = 1 # starting from 0 is not interesting
    EndFrame = numFrames-1 # note the -1
    k_list = np.linspace(StartFrame, EndFrame, Nfiles).astype(int) # form list of frames to print
    
    currentTime = 0.0
    
    Filelist = []
    
    for k in k_list:
        
        frame = step.frames[k]
    
        # update the current time if the above line passes
        currentTime = frame.frameValue
        
        myViewport.odbDisplay.setFrame(step=0, frame=k)
        myViewport.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
        
        for label in damageLabels:

            # Select damage label
            myViewport.odbDisplay.setPrimaryVariable(variableLabel=label, outputPosition=INTEGRATION_POINT, )

            # Print for selected plies
            # the name of the file must contain all the necessary information for postprocessing, e.g.:
            # damage label, ply number, total strains (from modelName), and current increment (currentTime)
            for i in Plies:

                # Select layer i
                myViewport.odbDisplay.setPrimarySectionPoint(activePly="PLY-" + str(i))
                myViewport.odbDisplay.basicOptions.setValues(sectionPointScheme=PLY_BASED)

                # Print layer i
                fileName=label+'_Ply' + str(i) + '_' + modelName + '_T' + str(int(currentTime*1e6))
                session.printToFile(fileName=fileName, format=PNG, canvasObjects=(myViewport, ))
                Filelist.append(fileName+'.png')
    
    return Filelist

#End of function extractor_DamagePatterns