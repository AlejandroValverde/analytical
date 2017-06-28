import math
import numpy as np
import os
import platform

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

from moduleCommon import *

### GO TO ABAQUS FOLDER
cwd = os.getcwd()

globalChangeDir(cwd, '-..-..-workfolder')

############################
# FUNCTIONS

class structtype():
    pass

def mergeInstances(model, instancesToMerge, newName):

    model.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
        instances=instancesToMerge, name=newName, originalInstances=SUPPRESS)

def buildRib(model, ribDesign):
    model.ConstrainedSketch(name='__profile__', sheetSize=1.5 * max(H, B))

    #Create shape
    #1 -> 2
    model.sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
                ribDesign.rib1, 0.0))
    model.sketches['__profile__'].HorizontalConstraint(
                addUndoState=False, entity=
                model.sketches['__profile__'].geometry[2])
    #2 -> 3
    model.sketches['__profile__'].Line(point1=(ribDesign.rib1, 0.0), point2=(
                ribDesign.rib1, ribDesign.rib2))
    model.sketches['__profile__'].VerticalConstraint(addUndoState=
                False, entity=model.sketches['__profile__'].geometry[3])
    model.sketches['__profile__'].PerpendicularConstraint(
                addUndoState=False, entity1=
                model.sketches['__profile__'].geometry[2], entity2=
                model.sketches['__profile__'].geometry[3])

    #3 -> 4
    model.sketches['__profile__'].Line(point1=(ribDesign.rib1, ribDesign.rib2), point2=(
                0.0, ribDesign.rib2))
    model.sketches['__profile__'].HorizontalConstraint(
                addUndoState=False, entity=
                model.sketches['__profile__'].geometry[4])

    #4 -> 1
    model.sketches['__profile__'].Line(point1=(0.0, ribDesign.rib2), point2=(
                0.0, 0.0))
    model.sketches['__profile__'].VerticalConstraint(addUndoState=
                False, entity=model.sketches['__profile__'].geometry[5])
    model.sketches['__profile__'].PerpendicularConstraint(
                addUndoState=False, entity1=
                model.sketches['__profile__'].geometry[4], entity2=
                model.sketches['__profile__'].geometry[5])

    #5 -> 6
    model.sketches['__profile__'].Line(point1=(ribDesign.a, ribDesign.a), point2=(
                ribDesign.rib1 - ribDesign.a, ribDesign.a))
    model.sketches['__profile__'].HorizontalConstraint(
                addUndoState=False, entity=
                model.sketches['__profile__'].geometry[6])

    #6 -> 7
    model.sketches['__profile__'].Line(point1=(ribDesign.rib1 - ribDesign.a, ribDesign.a), point2=(
                ribDesign.rib1 - ribDesign.a, ribDesign.rib2 - ribDesign.a))
    model.sketches['__profile__'].VerticalConstraint(addUndoState=
                False, entity=model.sketches['__profile__'].geometry[7])
    model.sketches['__profile__'].PerpendicularConstraint(
                addUndoState=False, entity1=
                model.sketches['__profile__'].geometry[6], entity2=
                model.sketches['__profile__'].geometry[7])

    #7 -> 8
    model.sketches['__profile__'].Line(point1=(ribDesign.rib1 - ribDesign.a, ribDesign.rib2 - ribDesign.a), point2=(
                ribDesign.a, ribDesign.rib2 - ribDesign.a))
    model.sketches['__profile__'].HorizontalConstraint(
                addUndoState=False, entity=
                model.sketches['__profile__'].geometry[8])

    #8 -> 5
    model.sketches['__profile__'].Line(point1=(ribDesign.a, ribDesign.rib2 - ribDesign.a), point2=(
                ribDesign.a, ribDesign.a))
    model.sketches['__profile__'].VerticalConstraint(addUndoState=
                False, entity=model.sketches['__profile__'].geometry[9])
    model.sketches['__profile__'].PerpendicularConstraint(
                addUndoState=False, entity1=
                model.sketches['__profile__'].geometry[8], entity2=
                model.sketches['__profile__'].geometry[9])

    #Shell creation
    model.Part(dimensionality=THREE_D, name='Rib', type=
        DEFORMABLE_BODY)
    model.parts['Rib'].BaseShell(sketch=
        model.sketches['__profile__'])
    del model.sketches['__profile__']

    #General set "all" for the part
    model.parts['Rib'].Set(faces=
        model.parts['Rib'].faces.findAt(((ribDesign.a/2,ribDesign.a/2,0.0),),)
        , name='rib-all') 

    #Section
    model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material='Mat_rib', name='Section-Rib', numIntPts=5, 
        poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
        thickness=ribDesign.ribt, thicknessField='', thicknessModulus=None, thicknessType=
        UNIFORM, useDensity=OFF)

    #Section assignment, assign material to part. This part is not a shell, therefore, no thickness is defined when specifying the section.
    model.parts['Rib'].SectionAssignment(offset=0.0, offsetField=''
        , offsetType=TOP_SURFACE, region=
        model.parts['Rib'].sets['rib-all'], sectionName=
        'Section-Rib', thicknessAssignment=FROM_SECTION)

    #Instance operations
    stepInX = ribDesign.L / (ribDesign.numberOfRibs)
    ribDesign.innerRibsXpos = [ribDesign.L - ((i-1)*stepInX) for i in range(1, int(ribDesign.numberOfRibs) + 1)]

    instances_ribs = ()
    for r in range(int(ribDesign.numberOfRibs)):
        instanceRib = model.rootAssembly.Instance(dependent=ON, name='Rib-inst'+str(r), part=
            model.parts['Rib']) #Create instance
        model.rootAssembly.translate(instanceList=('Rib-inst'+str(r), ), 
            vector=(0.0, 0.0, ribDesign.innerRibsXpos[r]))

        instances_ribs += (instanceRib, )

    return instances_ribs

def loadParameters(paraRead, fileName):

    file = open(fileName, 'r')

    lines = file.readlines()

    for i in range(int(int((len(lines))/2.0))):

        nameParater = lines[(i*2)]
        valueParater = lines[(2*i)+1]

        if platform.system() == 'Linux':

            valueParater = valueParater.replace('\r\n','')
            nameParater = nameParater.replace('\r\n','')            

        elif platform.system() == 'Windows':

            valueParater = valueParater.replace('\n','')
            nameParater = nameParater.replace('\n','')

        else:

            valueParater = valueParater.replace('\r\n','')
            nameParater = nameParater.replace('\r\n','')

            print('OS not recognized, assumed unix based')

        setattr(paraRead, nameParater, valueParater)

    file.close()

    return paraRead

#############################

paraRead = structtype()

# [paraRead] = loadParameters(paraRead, 'inputAbaqusAnalytical.txt')

file = open('inputAbaqusAnalytical.txt', 'r')

lines = file.readlines()

for i in range(int(int((len(lines))/2.0))):

    nameParater = lines[(i*2)]
    valueParater = lines[(2*i)+1]

    if platform.system() == 'Linux':

        valueParater = valueParater.replace('\r\n','')
        nameParater = nameParater.replace('\r\n','')            

    elif platform.system() == 'Windows':

        valueParater = valueParater.replace('\n','')
        nameParater = nameParater.replace('\n','')

    else:

        valueParater = valueParater.replace('\r\n','')
        nameParater = nameParater.replace('\r\n','')

        print('OS not recognized, assumed unix based')

    setattr(paraRead, nameParater, valueParater)

file.close()  

#Build simple wing box.
L = 200 #mm
B = 50 #mm
H = 30 #mm
t1 = 2 #mm
t2 = 2 #mm

E1 = 69000.0 #N/mm^2
v1 = 0.3269 # steel: 0.28 #Poisson's ratio for aluminum G = E / (2*(nu + 1) )
E2 = E1/float(paraRead.E1overE2) #69000.0 #N/mm^2
v2 = 0.3269 #Poisson's ratio for aluminum 0.3269 G = E / (2*(nu + 1) )

#Real materials - Aluminum for section 1 and Steel for section 2
# E1 = 200000.0 #N/mm^2
# v1 = 0.28 # steel: 0.28 0.3269 #Poisson's ratio for aluminum G = E / (2*(nu + 1) )
# # E2 = E1/float(paraRead.E1overE2) #69000.0 #N/mm^2
# E2 = 69000.0 #float(paraRead.E2) #N/mm^2
# v2 = 0.3269 #Poisson's ratio for aluminum 0.3269 G = E / (2*(nu + 1) )

meshSize = 4

forceMagnitude = -4000 #Total force


numberOfPointsWhereMeasureU2 = 30;

forceRelXPos = 0.5

x_load = -((B/2) - (forceRelXPos*B)) #The abaqus and the analytical model have opposite signs for the axes

x_sc = -float(paraRead.y_sc)

x_moment = x_load - x_sc

#Rib
ribDesign = structtype
ribDesign.ribt = t1
ribDesign.E = E1 * 100
ribDesign.v = v1
ribDesign.a = 10
ribDesign.numberOfRibs = int(paraRead.numberOfRibs)

#Load
numberOfNodesToApplyLoad = 1;
typeOfLoad = 'moment' #options: 'externalForceAndCoupling', 'distributedLoad', 'distributedLoadOnRib', 'moment'

############################

model = mdb.models['Model-1']

#Materials
model.Material(name='Mat1')
model.materials['Mat1'].Elastic(table=((E1, v1), ))
model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, material='Mat1', name='Section-Mat1', numIntPts=5, 
    poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
    thickness=t1, thicknessField='', thicknessModulus=None, thicknessType=
    UNIFORM, useDensity=OFF)

model.Material(name='Mat2')
model.materials['Mat2'].Elastic(table=((E2, v2), ))
model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, material='Mat2', name='Section-Mat2', numIntPts=5, 
    poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
    thickness=t2, thicknessField='', thicknessModulus=None, thicknessType=
    UNIFORM, useDensity=OFF)

model.Material(name='Mat_rib')
model.materials['Mat_rib'].Elastic(table=((ribDesign.E, ribDesign.v), ))

#Sketch for beam section
model.ConstrainedSketch(name='__profile__', sheetSize=1.5 * max(H, B))

#Create rectangular shape
model.sketches['__profile__'].Line(point1=(0.0, 0.0), 
    point2=(B, 0.0))
model.sketches['__profile__'].HorizontalConstraint(
    addUndoState=False, entity=
    model.sketches['__profile__'].geometry[2])
model.sketches['__profile__'].Line(point1=(B, 0.0), 
    point2=(B, H))
model.sketches['__profile__'].VerticalConstraint(addUndoState=
    False, entity=model.sketches['__profile__'].geometry[3])
model.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    model.sketches['__profile__'].geometry[2], entity2=
    model.sketches['__profile__'].geometry[3])
model.sketches['__profile__'].Line(point1=(B, H), 
    point2=(0.0, H))
model.sketches['__profile__'].HorizontalConstraint(
    addUndoState=False, entity=
    model.sketches['__profile__'].geometry[4])
model.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    model.sketches['__profile__'].geometry[3], entity2=
    model.sketches['__profile__'].geometry[4])
model.sketches['__profile__'].Line(point1=(0.0, H), 
    point2=(0.0, 0.0))
model.sketches['__profile__'].VerticalConstraint(addUndoState=
    False, entity=model.sketches['__profile__'].geometry[5])
model.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    model.sketches['__profile__'].geometry[4], entity2=
    model.sketches['__profile__'].geometry[5])

#Shell extrude
model.Part(dimensionality=THREE_D, name='Wing-box', type=
    DEFORMABLE_BODY)
model.parts['Wing-box'].BaseShellExtrude(depth=L, sketch=
    model.sketches['__profile__'])
del model.sketches['__profile__']

#Set for part with properties #1
model.parts['Wing-box'].Set(faces=
    model.parts['Wing-box'].faces.findAt(((0.0, H/2, L/2),), 
    ((B/2, 0.0, L/2),), 
    ((B/2, H, L/2),), 
    ), name='part_1')

#Set for part with properties #2
model.parts['Wing-box'].Set(faces=
    model.parts['Wing-box'].faces.findAt(((B, H/2, L/2),),), name='part_2')

#Assign section
model.parts['Wing-box'].SectionAssignment(offset=0.0, offsetField=
    '', offsetType=MIDDLE_SURFACE, region=
    model.parts['Wing-box'].sets['part_1'], sectionName=
    'Section-Mat1', thicknessAssignment=FROM_SECTION)

#Assign section
model.parts['Wing-box'].SectionAssignment(offset=0.0, offsetField=
    '', offsetType=MIDDLE_SURFACE, region=
    model.parts['Wing-box'].sets['part_2'], sectionName=
    'Section-Mat2', thicknessAssignment=FROM_SECTION)

#Create instance
model.rootAssembly.Instance(dependent=ON, name='Wing-box_instance', part=model.parts['Wing-box'])

################################################
#Rib
if ribDesign.numberOfRibs != 0:
    ribDesign.rib1 = B
    ribDesign.rib2 = H
    ribDesign.L = L

    instances_ribs = buildRib(model, ribDesign)

    mergeInstances(model, (model.rootAssembly.instances['Wing-box_instance'], ) + instances_ribs, 'BoxPlusRib')

    partForTheNext = 'BoxPlusRib'
    instanceForTheNext = 'BoxPlusRib-1'

else:
    partForTheNext = 'Wing-box'
    instanceForTheNext = 'Wing-box_instance'
################################################

#Meshing
#Global Mesh
model.parts[partForTheNext].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=meshSize) #fine mesh
model.parts[partForTheNext].generateMesh()

### Define BC
##############################
instanceToApplyLoadAndBC = model.rootAssembly.instances[instanceForTheNext]
rf = model.rootAssembly.ReferencePoint(point=(B/2, H/2, 0.0))
model.rootAssembly.Set(name='referencePoint', referencePoints=(model.rootAssembly.referencePoints[rf.id], ))

#Edges to apply boundary condition
model.rootAssembly.Set(edges = instanceToApplyLoadAndBC.edges.findAt(((0.0, H/2, 0.0),), 
    ((B/2, 0.0, 0.0),), 
    ((B/2, H, 0.0),),
    ((B, H/2, 0.0),), 
    ), name='fixed')

#Enable coupling condition
model.Coupling(controlPoint= model.rootAssembly.sets['referencePoint'], couplingType=
    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
    'Constraint', surface= model.rootAssembly.sets['fixed'], u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

#Fix reference point
model.DisplacementBC(name='fixed', createStepName='Initial', 
    region=model.rootAssembly.sets['referencePoint'], u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

###############################

### Apply load
###############################
#Create step for load
model.StaticStep(description=
    'Step for load, linear', name='load', previous=
    'Initial')

if typeOfLoad == 'distributedLoad' or typeOfLoad == 'distributedLoadOnRib':

    #Create a set where a the displacement is imposed 
    mdb.models['Model-1'].rootAssembly.regenerate()
    # 
    #Create set
    tupleOfNodes = ()

    xPosForLoad = np.linspace(L, 0, numberOfNodesToApplyLoad, endpoint = False)

    for x in xPosForLoad:

        FlagSearch = True
        radius = 1
        while FlagSearch:

            if typeOfLoad == 'distributedLoadOnRib':
                node = instanceToApplyLoadAndBC.nodes.getByBoundingSphere((forceRelXPos*B, H - (ribDesign.a/2), x), radius)
            else:
                node = instanceToApplyLoadAndBC.nodes.getByBoundingSphere((forceRelXPos*B, H, x), radius)
        
            if node: #Continue if node found

               #Apply force just to the first node found
               nodeTuple = node.getMask()
               maskStr0 = nodeTuple[0].split('#')
               maskStr = '[#'+str(maskStr0[1])+'#'+str(maskStr0[2]) + ']'
               tupleOfNodes += (instanceToApplyLoadAndBC.nodes.getSequenceFromMask(mask=(maskStr,),), )

               FlagSearch = False

            else:
               print('Node to apply force not found, radius increased')
               radius += 1

            if radius > 10:

               raise ValueError('Nodes where force should be applied not found')

        print('Final radius utilized for finding the node: '+str(radius))

    #Create set
    model.rootAssembly.Set(name='setForDistributedLoad', nodes=tupleOfNodes)

    #Define displacement condition
    model.ConcentratedForce(cf2=forceMagnitude/len(xPosForLoad), createStepName='load', 
        distributionType=UNIFORM, field='', localCsys=None, name='distributedLoad', 
        region=model.rootAssembly.sets['setForDistributedLoad'])

elif typeOfLoad == 'moment':
    ##############################

    #Moment load
    model.rootAssembly.Set(edges = instanceToApplyLoadAndBC.edges.findAt(((0.0, H/2, L),), 
        ((B/2, 0.0, L),), 
        ((B/2, H, L),),
        ((B, H/2, L),), 
        ), name='momentEdgesOnBox')

    model.rootAssembly.Set(edges = instanceToApplyLoadAndBC.edges.findAt(((ribDesign.a, H/2, L),), 
        ((B/2, ribDesign.a, L),), 
        ((B/2, H - ribDesign.a, L),),
        ((B - ribDesign.a, H/2, L),), 
        ), name='momentEdgesOnRib')

    rf2 = model.rootAssembly.ReferencePoint(point=(B/2, H/2, L))
    model.rootAssembly.Set(name='referencePointMoment', referencePoints=(model.rootAssembly.referencePoints[rf2.id], ))

    if ribDesign.numberOfRibs != 0:
        model.Coupling(controlPoint=
            model.rootAssembly.sets['referencePointMoment'], 
            couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, 
            name='Constraint-2', surface=
            model.rootAssembly.sets['momentEdgesOnRib'], u1=OFF, u2=OFF, u3=
            OFF, ur1=OFF, ur2=OFF, ur3=ON)
    else:
        model.Coupling(controlPoint=
            model.rootAssembly.sets['referencePointMoment'], 
            couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, 
            name='Constraint-2', surface=
            model.rootAssembly.sets['momentEdgesOnBox'], u1=OFF, u2=OFF, u3=
            OFF, ur1=OFF, ur2=OFF, ur3=ON)

    model.Moment(cm3=x_moment * forceMagnitude, createStepName='load', 
        distributionType=UNIFORM, field='', localCsys=None, name='moment', region=
        model.rootAssembly.sets['referencePointMoment'])

    #################################

elif typeOfLoad == 'externalForceAndCoupling':
    #Moment from the edge, reference and coupling
    displacementExternal = B/2

    moment = y_load * forceMagnitude

    forceExternal = -moment / displacementExternal

    rf_external = model.rootAssembly.ReferencePoint(point=(-displacementExternal, H/2, L))
    model.rootAssembly.Set(name='referencePointExternal', referencePoints=(model.rootAssembly.referencePoints[rf_external.id], ))

    #Slave
    model.rootAssembly.Set(edges = instanceToApplyLoadAndBC.edges.findAt(((0.0, H/2, L),), ), name='edgesForExternalMoment')

    mdb.models['Model-1'].Coupling(controlPoint=
        mdb.models['Model-1'].rootAssembly.sets['referencePointExternal'], 
        couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, 
        name='Constraint-2', surface=
        mdb.models['Model-1'].rootAssembly.sets['edgesForExternalMoment'], u1=OFF, u2=ON, u3=
        OFF, ur1=OFF, ur2=OFF, ur3=OFF)

    model.ConcentratedForce(cf2=forceExternal, createStepName='load', 
        distributionType=UNIFORM, field='', localCsys=None, name='externalLoadForMoment', 
        region=model.rootAssembly.sets['referencePointExternal']) #GIVES ERROR

else:

    raise ValueError('ERROR: Not correct load case selected')

######################################
#Job operations
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='Job_simpleBeam', nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

#Submit job
model.rootAssembly.regenerate()
mdb.jobs['Job_simpleBeam'].submit(consistencyChecking=OFF)
mdb.jobs['Job_simpleBeam'].waitForCompletion()

#Post processing
#Get current folder
cwd = os.getcwd()

#Check if postProc folder already exists
globalCreateDir(cwd, '-postProc')

#Create folder for simulation results
globalCreateDir(cwd, '-postProc-'+paraRead.Iter)
#Get current folder
cwd = os.getcwd()

o3 = session.openOdb(name='Job_simpleBeam.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
a = mdb.models['Model-1'].rootAssembly #in case we want to see the assembly
odb = session.odbs[cwd + '\\Job_simpleBeam.odb']

xPosForLoad = np.linspace(1, L, numberOfPointsWhereMeasureU2)

#Obtain data from path - U2
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))

#Move to simulation results folder
globalChangeDir(cwd, '-postProc-'+paraRead.Iter)

#Path to analyze rotation of points
for x in xPosForLoad:

    session.Path(name='pth_u2_'+str(x), type=POINT_LIST, expression=((0, H, x), (B, H, x)))
    pth_u2 = session.paths['pth_u2_'+str(x)]

    session.XYDataFromPath(name='XYData_u2_'+str(x), path=pth_u2, includeIntersections=False, 
        projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=30, 
        projectionTolerance=0, shape=UNDEFORMED, labelType=NORM_DISTANCE)

    #Write u2 data to XY file
    x0_u2 = session.xyDataObjects['XYData_u2_'+str(x)]
    session.writeXYReport(fileName='XYData_u2_'+str(x)+'.rpt', xyData=(x0_u2, ), appendMode=OFF)

#Obtain data from path - UR3
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UR', outputPosition=NODAL, refinement=(COMPONENT, 'UR3'))

session.Path(name='pth_ur3', type=POINT_LIST, expression=((B/2, H, 0), (B/2, H, L)))
pth_ur3 = session.paths['pth_ur3']
session.XYDataFromPath(name='XYData_ur3', path=pth_ur3, includeIntersections=False, 
        projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=30, 
        projectionTolerance=0, shape=UNDEFORMED, labelType=NORM_DISTANCE)

#Write ur3 data to XY file
x0_ur3 = session.xyDataObjects['XYData_ur3']

session.writeXYReport(fileName='XYData_ur3'+'.rpt', xyData=(x0_ur3, ), appendMode=OFF)

#Return to original working folder
globalChangeDir(cwd, '.')