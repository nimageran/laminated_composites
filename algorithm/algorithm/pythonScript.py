# -*- coding: mbcs -*-
"""A Python/MATLAB-based computational framework for optimization of various laminated composite plate models in Abaqus by A. Kaveh, A. Dadras, N. Geran Malek, R. Ansari."""
from __future__ import division 
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
import math
import numpy
import sys
from abaqus import *
from abaqusConstants import *
import __main__
import regionToolset
import displayGroupMdbToolset as dgm
import xyPlot
import displayGroupOdbToolset as dgo
import re
import mesh
import os 
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
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import fileinput


jj=0
################################################
# Counting the number of files in the directory (nFiles)
################################################
# Open a file
path = "D:/optLaminatedComp/JOB"
dirs = os.listdir( path )

# This would print all the files and directories
nFiles=0 # Number of files existed in the directory path
for file in dirs:
        print (file)
        nFiles=nFiles+1
print(nFiles)
################################################
# Making folder for the current Model: 
################################################
os.mkdir('D:\\optLaminatedComp\\JOB\\%d'%(jj+1+nFiles))

################################################
# Changing the default working directory to an arbitrary working directory :
################################################
os.chdir(r"D:\\optLaminatedComp\\JOB\\%d"%(jj+1+nFiles))

##################################### Regenerating  : 
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
referenceRepresentation=OFF)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, 
optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
meshTechnique=ON)
############################################################
##################################### Creating Part
############################################################
s = mdb.models['Model-%d'%(jj+1)].ConstrainedSketch(name='__profile__', sheetSize=2)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(decimalPlaces=3)
s.setPrimaryObject(option=STANDALONE)
L=0.508# a [m]
W=0.254 # b [m]
h=0.025 # Total thickness of the composite plate= 64*0.127*0.001
############################################################
#################################### Units used in Abaqus [N] and [m]
############################################################
s.rectangle(point1=(0.0, 0.0), point2=(L, W))
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
############################################################
###################### Regenerating the created meshes in the assembly : 
############################################################
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
############################################################
##################################### Creating Assembly :
############################################################
a1 = mdb.models['Model-1'].rootAssembly
a1.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a1.Instance(name='Part-1-1', part=p, dependent=OFF)
############################################################
##################################### Creating Step(Buckling):
############################################################
mdb.models['Model-1'].BuckleStep(name='Step-1', previous='Initial', numEigen=5, vectors=10)
# Changing number of Iterations in order to rectify the Error
mdb.models['Model-1'].steps['Step-1'].setValues(maxIterations=10000)
############################################################
##################################### Material Definition :
############################################################
nPartition_horizontal=1
nPartition_vertical=1
nCell=nPartition_vertical*nPartition_horizontal
for l in range(nCell):
    mdb.models['Model-1'].Material(name='Material_%d'%(l+1))
    mdb.models['Model-1'].materials['Material_%d'%(l+1)].Elastic(type=LAMINA, table=((127553009900, 13031091280, 0.3, 6412124283,6412124283, 6412124283), )) # [mpduli]=GPa
############################################################
##################################### Asigning Mesh Controls :
############################################################
nPartition_horizontal=1
nPartition_vertical=1
counter=0
for j in range(nPartition_horizontal):
		for k in range(nPartition_vertical):
			counter=k+j*(nPartition_vertical)
			print(counter)
			x=(2*k+1)*(L/nPartition_vertical)/2
			y=(2*j+1)*(W/nPartition_horizontal)/2
			p = mdb.models['Model-1'].parts['Part-1']
                        f = p.faces
			faces = f.findAt(((x,y,0), ))
                        a = mdb.models['Model-1'].rootAssembly
                        pickedRegions = faces
                        a.setMeshControls(regions=pickedRegions, elemShape=QUAD, technique=FREE)
############################################################
##################################### Meshing
############################################################
for k in range(nPartition_vertical):
        
        x=(2*k+1)*(L/nPartition_vertical)/2
        y=W
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
        pickedEdges = e1.findAt(((x,y,0), ))
	NN=20*(L/W)
	NN=int(NN)
        a.seedEdgeByNumber(edges=pickedEdges, number=NN, constraint=FINER)
        

for k in range(nPartition_vertical):
        
        x=(2*k+1)*(L/nPartition_vertical)/2
        y=0
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
        pickedEdges = e1.findAt(((x,y,0), ))
        a.seedEdgeByNumber(edges=pickedEdges, number=NN, constraint=FINER)
        
for j in range(nPartition_vertical):
        
        x=L
        y=(2*j+1)*(W/nPartition_horizontal)/2
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
        pickedEdges = e1.findAt(((x,y,0), ))
        a.seedEdgeByNumber(edges=pickedEdges, number=20, constraint=FINER)

for j in range(nPartition_vertical):
        
        x=0
        y=(2*j+1)*(W/nPartition_horizontal)/2
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
        pickedEdges = e1.findAt(((x,y,0), ))
        a.seedEdgeByNumber(edges=pickedEdges, number=20, constraint=FINER)
        
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['Part-1-1'], )
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a.generateMesh(regions=partInstances)

###############################################
####################################### Layup :
###############################################
counter=0
for j in range(nPartition_horizontal):
        for k in range(nPartition_vertical):
                counter=k+j*(nPartition_vertical)
		compositeLayup = mdb.models['Model-1'].parts['Part-1'].CompositeLayup(
                        name='CompositeLayup_%d'%(counter+1), description='', elementType=SHELL, 
                        offsetType=MIDDLE_SURFACE, symmetric=True, 
                        thicknessAssignment=FROM_SECTION)
		compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
                        useDensity=OFF)
		compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
                        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, axis=AXIS_3)
		print(counter+1)
                x=(2*k+1)*(L/nPartition_vertical)/2
                y=(2*j+1)*(W/nPartition_horizontal)/2
                print(x)
                print(y)
                f = p.faces                
                faces = f.findAt(((x,y,0), ))
                layupOrientation = None
                p = mdb.models['Model-1'].parts['Part-1']
                region1=regionToolset.Region(faces=faces)
                region2=regionToolset.Region(faces=faces)
                region3=regionToolset.Region(faces=faces)

                ##############################################################
                ##############################################################
                f=open('D:\\optLaminatedComp\\TEXT2.txt')
                OrientationsList = []
                for l in f:
                        
                        row = l.split()
                        print(row)
                        OrientationsList=row
                OrientationsList = list(map(float, OrientationsList))
                print(OrientationsList)
                myOrientationsList=OrientationsList

                for iii in range(len(myOrientationsList)):

                        compositeLayup.CompositePly(suppressed=False, plyName='Ply-%d'%(iii+1), 
                        region=region1, material='Material_%d'%(counter+1), thicknessType=SPECIFY_THICKNESS, 
                        thickness=0.127E-3, orientationType=SPECIFY_ORIENT, orientationValue=myOrientationsList[iii], 
                        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                        axis=AXIS_3, angle=0.0, numIntPoints=3)
                
##########################################################################
##################################### Creating a Boundary Conditions :
##########################################################################
# Right
for j in range(nPartition_vertical):
        
        x=L
        y=(2*j+1)*(W/nPartition_horizontal)/2
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
	edges1 = e1.findAt(((x,y,0), ))
	region = regionToolset.Region(edges=edges1)
	mdb.models['Model-1'].DisplacementBC(name='BC-%d'%(j+6), createStepName='Initial', region=region, u1=UNSET, u2=UNSET, u3=0.0, ur1=0, ur2=UNSET, ur3=0, 
		amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
# Left
for j in range(nPartition_vertical):
        
        x=0
        y=(2*j+1)*(W/nPartition_horizontal)/2
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
	edges1 = e1.findAt(((x,y,0), ))
	region = regionToolset.Region(edges=edges1)
	mdb.models['Model-1'].DisplacementBC(name='BC-%d'%(j+20), createStepName='Initial', region=region, u1=0, u2=UNSET, u3=0, ur1=0, ur2=UNSET, ur3=0, 
		amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
# Top
for k in range(nPartition_vertical):
        
        x=(2*k+1)*(L/nPartition_vertical)/2
        y=W
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
        pickedEdges = e1.findAt(((x,y,0), ))
	edges1 = e1.findAt(((x,y,0), ))
	region = regionToolset.Region(edges=edges1)
	mdb.models['Model-1'].DisplacementBC(name='BC-%d'%(k+40), createStepName='Initial', region=region, u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=0, ur3=0, 
		amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
# Bottom
for k in range(nPartition_vertical):
        
        x=(2*k+1)*(L/nPartition_vertical)/2
        y=0
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
        pickedEdges = e1.findAt(((x,y,0), ))
	edges1 = e1.findAt(((x,y,0), ))
	region = regionToolset.Region(edges=edges1)
	mdb.models['Model-1'].DisplacementBC(name='BC-%d'%(k+60), createStepName='Initial', region=region, u1=UNSET, u2=0, u3=0.0, ur1=UNSET, ur2=0, ur3=0, 
		amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
############################################################
######## Creating an element set for whole model :
############################################################
a = mdb.models['Model-1'].rootAssembly
els_1 = a.instances['Part-1-1'].elements
Coh_els = els_1.getByBoundingBox(-3 ,-3 ,-3 ,3 ,3 ,3) 
a.Set(elements=Coh_els, name='myElements')

############################################################
############################## Applying shell edge load :
############################################################
    
## on Right Edge
for j in range(nPartition_vertical):
        
        x=L
        y=(2*j+1)*(W/nPartition_horizontal)/2
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
	edges1 = e1.findAt(((x,y,0), ))
        region = a.Surface(side1Edges=edges1, name='SurfBCRight-%d'%(j+1))
	mdb.models['Model-1'].ShellEdgeLoad(name='loadBCRight-%d'%(j+1), createStepName='Step-1', region=region, magnitude=1,
                                            distributionType=UNIFORM, field='', localCsys=None)

## on Left Edge
for j in range(nPartition_vertical):
        
		
	x=0
        y=(2*j+1)*(W/nPartition_horizontal)/2
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
	edges1 = e1.findAt(((x,y,0), ))
        region = a.Surface(side1Edges=edges1, name='SurfBCLeft-%d'%(j+1))
	mdb.models['Model-1'].ShellEdgeLoad(name='loadBCLeft-%d'%(j+1), createStepName='Step-1', region=region, magnitude=1,
                                            distributionType=UNIFORM, field='', localCsys=None)

## on Top Edge
for j in range(nPartition_vertical):
        
		
	x=(2*j+1)*(L/nPartition_horizontal)/2
        y=W
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
	edges1 = e1.findAt(((x,y,0), ))
        region = a.Surface(side1Edges=edges1, name='SurfBCTop-%d'%(j+1))
	mdb.models['Model-1'].ShellEdgeLoad(name='SurfBCTop-%d'%(j+1), createStepName='Step-1', region=region, magnitude=1,
                                            distributionType=UNIFORM, field='', localCsys=None)

## on Bottom Edge
for j in range(nPartition_vertical):
        
		
	x=(2*j+1)*(L/nPartition_horizontal)/2
        y=0
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['Part-1-1'].edges
	edges1 = e1.findAt(((x,y,0), ))
        region = a.Surface(side1Edges=edges1, name='SurfBCBottom-%d'%(j+1))
	mdb.models['Model-1'].ShellEdgeLoad(name='SurfBCBottom-%d'%(j+1), createStepName='Step-1', region=region, magnitude=1,
                                            distributionType=UNIFORM, field='', localCsys=None)

###############################################################
#################################### Applied equation constraint on top edge :
###############################################################
### Defining Equation Constraint :
# Creating a Set for The Top Edge :
mdb.models['Model-1'].rootAssembly.Set(edges=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.findAt(((L/2,W,0), )), name='TOPSET')
# Creating a Reference Point for The Top Edge :
a = mdb.models['Model-1'].rootAssembly
myRefPoint=a.ReferencePoint(point=(L/2, W, 0.0))

# Creating a Set for The Reference Point of The Top Edge :
RP_id=myRefPoint.id # # Creating a Set for Reference Point Needs the id of Reference Point, so we Create it with ((id)) syntax .
refpoint1= (mdb.models['Model-1'].rootAssembly.referencePoints[RP_id], ) # Creating a Set for Reference Point Needs a Tuple, So We Create a Tuple
a = mdb.models['Model-1'].rootAssembly
a.Set(referencePoints=refpoint1, name='REFPOINTSET_TOP')
# Creating Top Edge Equation Constraint :
mdb.models['Model-1'].Equation(name='EquationConst_TOP', terms=((1.0, 'TOPSET', 2), (-1.0, 'REFPOINTSET_TOP', 2)))

print('TOP Equation Constraint Has Been Created')


###############################################################
#################################### Applied equation constraint on right edge :
###############################################################
### Defining Equation Constraint :
# Creating a Set for The Right Edge :
mdb.models['Model-1'].rootAssembly.Set(edges=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.findAt(((L,W/2,0), )), name='RIGHTSET')
# Creating a Reference Point for The Right Edge :
a = mdb.models['Model-1'].rootAssembly
myRefPoint=a.ReferencePoint(point=(L, W/2, 0.0))

# Creating a Set for The Reference Point of The Right Edge :
RP_id=myRefPoint.id # # Creating a Set for Reference Point Needs the id of Reference Point, so we Create it with ((id)) syntax .
refpoint1= (mdb.models['Model-1'].rootAssembly.referencePoints[RP_id], ) # Creating a Set for Reference Point Needs a Tuple, So We Create a Tuple
a = mdb.models['Model-1'].rootAssembly
a.Set(referencePoints=refpoint1, name='REFPOINTSET_RIGHT')
# Creating Right Edge Equation Constraint :
mdb.models['Model-1'].Equation(name='EquationConst_RIGHT', terms=((1.0, 'RIGHTSET', 1), (-1.0, 'REFPOINTSET_RIGHT', 1)))

print('RIGHT Equation Constraint Has Been Created')


jj=0
##################################### Creating and submiting the current job 
mdb.Job(name='Job-%d'%(jj+1+nFiles), model='Model-%d'%(jj+1), description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=False, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)
#################################################################
############# Saving the completed job in arbitrary folder path :
#################################################################

mdb.saveAs(pathName='D:/optLaminatedComp/JOB/%d'%(jj+1+nFiles) + '/Job-%d'%(jj+1+nFiles))
##################################################
############# Submiting the current job :
##################################################
mdb.jobs['Job-%d'%(jj+1+nFiles)].submit(consistencyChecking=OFF)
##################################################
############# Waiting for completion of the current job
##################################################
mdb.jobs['Job-%d'%(jj+1+nFiles)].waitForCompletion()

##################################################
###################   Result extraction : Buckling Loads 
##################################################
import odbAccess
a1=odbAccess.openOdb('Job-%d.odb'%(jj+1+nFiles))
a1.steps['Step-1'].frames[1].mode
a1.steps['Step-1'].frames[1].description
f1=a1.steps['Step-1'].frames[1].description
First_Buckling_Load=float(f1[28:48])


########################################################################
########################################################################
#non_Dimentional_First_Buckling_Load=First_Buckling_Load*((W**2)/(2.1e10*(h**3)))
########################################################################
########################################################################
file_id1=open('D:/optLaminatedComp/JOB/%d'%(jj+1+nFiles) + '/loadFactor.txt','w')
file_id1.write(str(First_Buckling_Load))
file_id1.close()
########################################################################
########################################################################
#Writing Ply Orientations for every simulations into the associated folder of that Job 
########################################################################
########################################################################
file_id1=open('D:/optLaminatedComp/JOB/%d'%(jj+1+nFiles) + '/thisJobOrientations.txt','w')
file_id1.write(str(OrientationsList))
file_id1.close()

