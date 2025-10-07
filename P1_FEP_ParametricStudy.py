"""
=======================================================================
 Title:       Flush End-Plate Beam-to-Column Joints – Seismic and Fire Performance
 File:        P1_FEP_ParametricStudy.py
 Version:     1.0
 Date:        2025-10-07
=======================================================================

 Authors:
     Md. Ibrahim Kholil
         Postgraduate Student, Department of Civil Engineering,
         Khulna University of Engineering & Technology (KUET), Bangladesh
         Email: engikholil@gmail.com
         ORCID: https://orcid.org/0000-0002-3349-8496

     Muhammad Harunur Rashid
         Professor, Department of Civil Engineering,
         Khulna University of Engineering & Technology (KUET), Bangladesh
         ORCID: https://orcid.org/0000-0002-3656-5025

     Aziz Ahmed
         Senior Lecturer, School of Civil, Mining, Environmental and Architectural Engineering,
         University of Wollongong, Australia
         ORCID: https://orcid.org/0000-0001-9707-2606

=======================================================================
 Description:
     This script automates the creation, execution, and post-processing of
     finite element (FE) models in Abaqus CAE for simulating steel
     Flush End-Plate (FEP) beam-to-column joints under:
         • Monotonic loading (static performance)
         • Cyclic loading (seismic simulation)
         • Thermal-mechanical coupling (fire exposure)

     It represents the numerical backbone of the study:
         "Numerical Investigation of Monotonic Seismic and Fire Performance
          of Flush End-Plate Beam-to-Column Joints" (Kholil et al., 2025).

=======================================================================
 Key Features:
     • Automated geometry generation, meshing, and material assignment
     • Partitioning of components for stress-strain monitoring
     • Bolt, plate, and column parameter variability (t, D, spacing)
     • Contact, tie, and constraint definitions (rigid and deformable)
     • Application of thermal and mechanical boundary conditions
     • Batch job creation and output summary (ultimate load, rotation, etc.)

=======================================================================
 Usage:
     Run from Abaqus Command Line Interface (CLI):

         abaqus cae noGUI=P1_FEP_ParametricStudy.py

     Configure model parameters in the “User Parameter Section”
     inside this script before running.

=======================================================================
 Outputs:
     • Abaqus model files (.cae, .inp)
     • Job logs and output databases (.odb)
     • CSV summary file containing key response metrics:
         – Ultimate load
         – Rotation capacity
         – Maximum temperature
         – Failure mode and stiffness degradation

=======================================================================
 Compatibility:
     • Abaqus/CAE 2021–2024
     • Python 3.x (via Abaqus scripting interface)
=======================================================================
 License:
     MIT License (© 2025 Md. Ibrahim Kholil, Muhammad Harunur Rashid, and Aziz Ahmed)
     See LICENSE file for full license text.
=======================================================================
"""


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
from abaqus import *
from abaqusConstants import *
from odbAccess import *
import mesh   # Import the mesh module
from abaqus import *
from abaqusConstants import *
from assembly import *
from regionToolset import *
import numpy as np
import math

#import os
#os.chdir(r"E:\Final Research Value & model\Result data new")

#Column Parameter
myC_FlangeTop_W = 120    #Column width
myC_FlangeTop_T = 12     #Column Thickness
myC_FlangeBottom_W = 120    
myC_FlangeBotom_T = 12
myC_Depth = 240
myC_Web_H = myC_Depth - (myC_FlangeTop_T + myC_FlangeBotom_T)  #just web height not in total
myC_Web_T = 10
myC_H = 1500     #Column Height

#Beam Paraeter
myB_FlangeTop_W = 120    #Beam width
myB_FlangeTop_T = 12     #Beam flange Thickness
myB_FlangeBottom_W = 120    
myB_FlangeBotom_T = 12
myB_Depth = 240
myB_Web_H = myB_Depth-(myB_FlangeTop_T + myB_FlangeBotom_T)    #Just web height not in total
myB_Web_H_cc = (myB_Web_H+(myB_FlangeTop_T+myB_FlangeBotom_T)/2)  ##Just top to bottom flange CC distance
myB_Web_T = 10
myB_H = 1500     #Beam Height


myBolt_k = 10.0 #outer thickness
MyBolt_S = 24.0 #outer dia of T & B
MyBolt_D = 16.0
MyBoltClear = 2.0

#Bolt Parameter
myBolt_T_Dia = MyBolt_S
myBolt_T_T = myBolt_k           #top thickness
myBolt_M_Dia = MyBolt_D
myBolt_B_Dia = MyBolt_S
myBolt_B_T = myBolt_k              #Bottom Thickness

#End Plate Parameter

##################################### For Z Varies##################################################
##################################### For Z Varies##################################################

myEndPlate_W = 120    #width ep_w
myEndPlate_H = 260     #Height ep_h
myEndPlate_T = 8      #Thickness ept
myEndPlate_H_CC_T = 70   # Horigontal dstanc of two bolt hole center to certer in top ep_cc_t
myEndPlate_H_CC_B = 70   # Horigontal dstanc of two bolt hole center to certer in Bottom
Z_origonal = 179
New_Z = 179
Z_vary = New_Z-Z_origonal
myEndPlate_T_C = 65-Z_vary       # Plate top to bolt hole center ep_tc
myEndPlate_B_C = 65       # Plate Bottom to bolt hole center ep_tc_1
myBoltHoleDia = myBolt_M_Dia+MyBoltClear     # bh_d
myBolt_M_T =  myEndPlate_T +  myC_FlangeTop_T    #Middle Thickness
myEP_V_D_second_Row = 65+Z_vary  #1st row to second row ep_vccb_1
myEP_V_D_Third_Row = 130+Z_vary #1st row to Third row ep_vccb_2
Cc_V = (myEP_V_D_Third_Row-myEP_V_D_second_Row)   #2nd to 3rd row 

##################################### For Z Varies##################################################
##################################### For Z Varies##################################################

#Loading
myLoad_D = 1470    #with loading point to column edge surface 

#Material
myE = 200000
myPoiratio = 0.3
myDensity = 8.5e-09
myFy = 355
myFu = 470   
mySry = 0.0
mySru = 0.18



myJobmodelname = "Column_Trial_5"
myPart_1 = "Steel Column"
myPart_2 = "Steel Beam"
myPart_3 = "End Plate"
myPart_4 = "Bolt"

myMaterial_1 = "Flange"
myMaterial_2 = "Web"
myMaterial_3 = "End Plate"
myMaterial_4 = "Bolt"

myCS_1_1 = "Column Flange CS" 
myCS_1_2 = "Column Web CS" 
myCS_2 = "Beam Web CS"
myCS_2_1 = "Beam Flange_CS"
myCS_3 = "End Plate CS"
myCS_4 = "Bolt CS"

myInstance_1 = "Steel Column"
myInstance_2 = "Steel Beam"
myInstance_3 = "End Plate"
myInstance_4 = "Bolt"


myString = myJobmodelname
mdb.Model(name=myString)
myJobName= myJobmodelname

#My Material 
myDensity = 7.85e-9
MyBoltEM = 191500
MyFlangeEM = 200000
MyWebEM = 200000
MyEPEM = 198000

MyBoltPlastic = (588.480037, 0.0), (620.221932, 0.001969588), (639.92125, 0.003119448), (660.066409, 0.004869265), (680.884675, 0.007497904), (702.708338, 0.011397498), (726.018148, 0.017111364), (751.50272, 0.025380009), (780.139135, 0.037193653), (813.301382, 0.053846562), (852.905183, 0.076983358), (901.6, 0.108620591), (1610.0, 0.68473987)
MyFlangePlastic = (324.883797, 0.0), (342.263405, 0.001986836), (372.640806, 0.004168264), (403.927543, 0.008251454), (436.871142, 0.015530915), (472.714356, 0.027949365), (513.453487, 0.048276771), (562.200436, 0.080231161), (623.679429, 0.128434998), (704.895798, 0.198070686), (816.021472, 0.294140858), (971.55, 0.420409985), (1270.0, 0.686797181)
MyWebPlastic = (324.883797, 0.0), (342.263405, 0.001986836), (372.640806, 0.004168264), (403.927543, 0.008251454), (436.871142, 0.015530915), (472.714356, 0.027949365), (513.453487, 0.048276771), (562.200436, 0.080231161), (623.679429, 0.128434998), (704.895798, 0.198070686), (816.021472, 0.294140858), (971.55, 0.420409985), (1270.0, 0.686797181)
MyEPPlastic = (319.160569, 0.0), (336.231125, 0.001987108), (366.703443, 0.004193395), (398.090211, 0.00833728), (431.147856, 0.015741831), (467.135196, 0.028390912), (508.07507, 0.049107916), (557.120339, 0.081671038), (619.054958, 0.130755713), (700.967089, 0.201566892), (813.138202, 0.299072753), (970.2, 0.426931416), (1260.0, 0.686847181)

def Create_Material(model,mats,density,elastic,plastic):
    mdb.models[model].Material(name=mats)
    mdb.models[model].materials[mats].Density(table=((density, ), ))
    mdb.models[model].materials[mats].Elastic(table=((elastic, 0.3), ))
    mdb.models[model].materials[mats].Plastic(scaleStress=None, table=(plastic))

Create_Material(myString,myMaterial_1,myDensity,MyFlangeEM,MyFlangePlastic)
Create_Material(myString,myMaterial_2,myDensity,MyWebEM,MyWebPlastic)
Create_Material(myString,myMaterial_3,myDensity,MyEPEM,MyEPPlastic)
Create_Material(myString,myMaterial_4,myDensity,MyBoltEM,MyBoltPlastic)



def Create_Column(model,part,flangetw,fwdepth,webh,webt,flangebw,height):
    s1 = mdb.models[model].ConstrainedSketch(name='__profile__',  sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.Spot(point=(-flangetw/2, fwdepth/2))
    s1.Spot(point=(flangetw/2, fwdepth/2))
    s1.Spot(point=(-flangetw/2, webh/2))
    s1.Spot(point=(-webt/2, webh/2))
    s1.Spot(point=(webt/2, webh/2))
    s1.Spot(point=(flangetw/2, webh/2))
    s1.Spot(point=(-flangebw/2, -webh/2))
    s1.Spot(point=(-webt/2, -webh/2))
    s1.Spot(point=(webt/2, -webh/2))
    s1.Spot(point=(flangebw/2, -webh/2))
    s1.Spot(point=(-flangebw/2, -fwdepth/2))
    s1.Spot(point=(flangebw/2, -fwdepth/2))
    s1.Line(point1=(-flangetw/2, fwdepth/2), point2=(flangetw/2, fwdepth/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(flangetw/2, fwdepth/2), point2=(flangetw/2, webh/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(webt/2, webh/2), point2=(flangetw/2, webh/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(webt/2, webh/2), point2=(webt/2, -webh/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(flangebw/2, -webh/2), point2=(webt/2, -webh/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(flangebw/2, -webh/2), point2=(flangebw/2, -fwdepth/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(-flangebw/2, -fwdepth/2), point2=(flangebw/2, -fwdepth/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(-flangebw/2, -fwdepth/2), point2=(-flangebw/2, -webh/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(-webt/2, -webh/2), point2=(-flangebw/2, -webh/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(-webt/2, -webh/2), point2=(-webt/2, webh/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(-flangetw/2, fwdepth/2), point2=(-flangetw/2, webh/2))
    s1.VerticalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(-flangetw/2, webh/2), point2=(-webt/2, webh/2))
    s1.VerticalConstraint(entity=g[2], addUndoState=False)
    p = mdb.models[myString].Part(name=part, dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
    p = mdb.models[myString].parts[part]
    p.BaseSolidExtrude(sketch=s1, depth=height)
    s1.unsetPrimaryObject()


Create_Column(myString,myPart_1,myC_FlangeTop_W,myC_Depth,myC_Web_H,myC_Web_T,myC_FlangeBottom_W,myC_H)
#Create_Column(myString,myPart_2,myB_FlangeTop_W,myB_Depth,myB_Web_H,myB_Web_T,myB_FlangeBottom_W,myB_H)


def Create_Beam(model,part,flanget_w,web_h,flanget_t,flangeb_w,flangeb_t,beam_h):
    s1 = mdb.models[model].ConstrainedSketch(name='__profile__', 
        sheetSize=300.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.Spot(point=(-flanget_w/2, (web_h+flanget_t)/2))
    s1.Spot(point=(0.0, (web_h+flanget_t)/2))
    s1.Spot(point=(flanget_w/2, (web_h+flanget_t)/2))
    s1.Spot(point=(-flangeb_w/2, -(web_h+flangeb_t)/2))
    s1.Spot(point=(0.0, -(web_h+flangeb_t)/2))
    s1.Spot(point=(flangeb_w/2, -(web_h+flangeb_t)/2))
    s1.Line(point1=(-flanget_w/2, (web_h+flanget_t)/2), point2=(0.0, (web_h+flanget_t)/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(0.0, (web_h+flanget_t)/2), point2=(flanget_w/2, (web_h+flanget_t)/2))
    s1.HorizontalConstraint(entity=g[3], addUndoState=False)
    s1.ParallelConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
    s1.Line(point1=(0.0, (web_h+flanget_t)/2), point2=(0.0, -(web_h+flangeb_t)/2))
    s1.VerticalConstraint(entity=g[4], addUndoState=False)
    s1.Line(point1=(-flangeb_w/2, -(web_h+flangeb_t)/2), point2=(0.0, -(web_h+flangeb_t)/2))
    s1.HorizontalConstraint(entity=g[5], addUndoState=False)
    s1.Line(point1=(0.0, -(web_h+flangeb_t)/2), point2=(flangeb_w/2, -(web_h+flangeb_t)/2))
    s1.HorizontalConstraint(entity=g[6], addUndoState=False)
    p = mdb.models[model].Part(name=part, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models[model].parts[part]
    p.BaseShellExtrude(sketch=s1, depth=beam_h)
    s1.unsetPrimaryObject()


Create_Beam(myString,myPart_2,myB_FlangeTop_W,myB_Web_H,myB_FlangeTop_T,myB_FlangeBottom_W,myB_FlangeBotom_T,myB_H)

#del mdb.models['Model-1']



def Create_End_Plate(model,part,ep_w,ep_h,ep_cc_t,ep_tc,bh_d,ep_vccb_1,ep_vccb_2,ept,ep_tc_1):
    s1 = mdb.models[model].ConstrainedSketch(name='__profile__', 
        sheetSize=1000.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.Spot(point=(-ep_w/2, ep_h/2))
    s1.Spot(point=(ep_w/2, ep_h/2))
    s1.Spot(point=(-ep_w/2, -ep_h/2))
    s1.Spot(point=(ep_w/2, -ep_h/2))
    s1.Line(point1=(-ep_w/2, ep_h/2), point2=(ep_w/2, ep_h/2))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    s1.Line(point1=(ep_w/2, ep_h/2), point2=(ep_w/2, -ep_h/2))
    s1.VerticalConstraint(entity=g[3], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
    s1.Line(point1=(-ep_w/2, -ep_h/2), point2=(ep_w/2, -ep_h/2))
    s1.HorizontalConstraint(entity=g[4], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    s1.Line(point1=(-ep_w/2, -ep_h/2), point2=(-ep_w/2, ep_h/2))
    s1.VerticalConstraint(entity=g[5], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
    s1.CircleByCenterPerimeter(center=(-ep_cc_t/2, (ep_h/2)-ep_tc), point1=(-(ep_cc_t/2)+bh_d/2, (ep_h/2)-ep_tc))
    s1.CircleByCenterPerimeter(center=(ep_cc_t/2, (ep_h/2)-ep_tc), point1=((ep_cc_t/2)+bh_d/2, (ep_h/2)-ep_tc))
    s1.linearPattern(geomList=(g[6], g[7]), vertexList=(), number1=1, 
        spacing1=100.0, angle1=0.0, number2=2, spacing2=ep_vccb_1, angle2=270.0)
    s1.linearPattern(geomList=(g[6], g[7]), vertexList=(), number1=1, 
        spacing1=100.0, angle1=0.0, number2=2, spacing2=ep_vccb_2, angle2=270.0)
    p = mdb.models[model].Part(name=part, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models[model].parts[part]
    p.BaseSolidExtrude(sketch=s1, depth=ept)
    s1.unsetPrimaryObject()


Create_End_Plate(myString,myPart_3,myEndPlate_W,myEndPlate_H,myEndPlate_H_CC_T,myEndPlate_T_C,myBoltHoleDia,myEP_V_D_second_Row,myEP_V_D_Third_Row,myEndPlate_T,myEndPlate_B_C)



def Create_Bolt(model,part,m_d,m_t,bt_d,bb_d,b_t_t,b_b_t):
    s1 = mdb.models[model].ConstrainedSketch(name='__profile__', 
        sheetSize=300.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(m_d/2, 0.0))
    p = mdb.models[model].Part(name=part, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models[model].parts[part]
    p.BaseSolidExtrude(sketch=s1, depth=m_t)
    s1.unsetPrimaryObject()
    f, e = p.faces, p.edges
    t = p.MakeSketchTransform(sketchPlane=f[1], sketchUpEdge=e[0], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, myBolt_M_T))
    s1 = mdb.models[model].ConstrainedSketch(name='__profile__', 
        sheetSize=103.86, gridSpacing=2.59, transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(bt_d/2, 0.0))
    f1, e1 = p.faces, p.edges
    p.SolidExtrude(sketchPlane=f1[1], sketchUpEdge=e1[0], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s1, depth=b_t_t, flipExtrudeDirection=OFF)
    s1.unsetPrimaryObject()
    t = p.MakeSketchTransform(sketchPlane=f[4], sketchUpEdge=e[3], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0))
    s1 = mdb.models[model].ConstrainedSketch(name='__profile__', 
        sheetSize=103.86, gridSpacing=2.59, transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(bb_d/2, 0.0))
    p.SolidExtrude(sketchPlane=f1[4], sketchUpEdge=e1[3], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s1, depth=b_b_t, flipExtrudeDirection=OFF)
    s1.unsetPrimaryObject()


Create_Bolt(myString,myPart_4,myBolt_M_Dia,myBolt_M_T,myBolt_T_Dia,myBolt_B_Dia,myBolt_T_T,myBolt_B_T)


#------------------------------------------------------------------------------

# Create_Datum_Plane
def Create_Datum_Plane(type_plane,part,model,offset_plane):
    p = mdb.models[model].parts[part]
    myPlane = p.DatumPlaneByPrincipalPlane(principalPlane=type_plane, offset=offset_plane)
    myID = myPlane.id
    return myID

#------------------------------------------------------------------------------


def Create_Partion(model,part,id_plane):
    p = mdb.models[model].parts[part]
    c = p.cells[:]
    d = p.datums
    p.PartitionCellByDatumPlane(datumPlane=d[id_plane], cells=c)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

myID_0 = Create_Datum_Plane(XYPLANE,myPart_4,myString,0.0)
myID_1 = Create_Datum_Plane(XYPLANE,myPart_4,myString,myBolt_M_T)
myID_2 = Create_Datum_Plane(YZPLANE,myPart_4,myString,0.0)
myID_3 = Create_Datum_Plane(XZPLANE,myPart_4,myString,0.0)
myID_4 = Create_Datum_Plane(XZPLANE,myPart_1,myString,((myC_Web_H/2)))
myID_5 = Create_Datum_Plane(XZPLANE,myPart_1,myString,-((myC_Web_H/2)))
myID_6 = Create_Datum_Plane(XYPLANE,myPart_1,myString,myC_H/2)
myID_7 = Create_Datum_Plane(XYPLANE,myPart_1,myString,((myC_H/2-myEP_V_D_second_Row)))
myID_8 = Create_Datum_Plane(YZPLANE,myPart_3,myString,myEndPlate_H_CC_B/2)
myID_9 = Create_Datum_Plane(YZPLANE,myPart_3,myString,-myEndPlate_H_CC_B/2)
myID_10 = Create_Datum_Plane(XZPLANE,myPart_3,myString,myEP_V_D_second_Row)
myID_11 = Create_Datum_Plane(XZPLANE,myPart_3,myString,-Cc_V)
myID_12 = Create_Datum_Plane(YZPLANE,myPart_3,myString,0.0)
myID_13 = Create_Datum_Plane(XZPLANE,myPart_3,myString,0.0)
myID_14 = Create_Datum_Plane(XYPLANE,myPart_1,myString,((myC_H/2+Cc_V)))
myID_15 = Create_Datum_Plane(XYPLANE,myPart_1,myString,(((myC_H/2)+(myEndPlate_H/2))))
myID_16 = Create_Datum_Plane(XYPLANE,myPart_1,myString,(((myC_H/2)-(myEndPlate_H/2))))
myID_17 = Create_Datum_Plane(YZPLANE,myPart_1,myString,((myEndPlate_H_CC_B/2)))
myID_18 = Create_Datum_Plane(YZPLANE,myPart_1,myString,(-(myEndPlate_H_CC_B/2)))
myID_19 = Create_Datum_Plane(XYPLANE,myPart_2,myString,((myB_H/6)))
myID_20 = Create_Datum_Plane(XYPLANE,myPart_2,myString,((myLoad_D-myEndPlate_T)))
#myID_12 = Create_Datum_Plane(XZPLANE,myPart_3,myString,myEndPlate_Forth_Row/2)
#myID_13 = Create_Datum_Plane(XZPLANE,myPart_3,myString,-myEndPlate_Forth_Row/2)


#------------------------------------------------------------------------------

#Create_Partion
Create_Partion(myString,myPart_4,myID_0)
Create_Partion(myString,myPart_4,myID_1)
Create_Partion(myString,myPart_4,myID_2)
Create_Partion(myString,myPart_4,myID_3)
Create_Partion(myString,myPart_1,myID_4)
Create_Partion(myString,myPart_1,myID_5)
Create_Partion(myString,myPart_1,myID_6)
Create_Partion(myString,myPart_1,myID_7)
Create_Partion(myString,myPart_3,myID_8)
Create_Partion(myString,myPart_3,myID_9)
Create_Partion(myString,myPart_3,myID_10)
Create_Partion(myString,myPart_3,myID_11)
Create_Partion(myString,myPart_3,myID_12)
Create_Partion(myString,myPart_3,myID_13)
Create_Partion(myString,myPart_1,myID_14)
Create_Partion(myString,myPart_1,myID_15)
Create_Partion(myString,myPart_1,myID_16)
Create_Partion(myString,myPart_1,myID_17)
Create_Partion(myString,myPart_1,myID_18)
#Create_Partion(myString,myPart_3,myID_12)
#Create_Partion(myString,myPart_3,myID_13)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#def Cut_Extrude_Column(model,part,ep_w,c_h,eph_cc_b,bh_d,cft_t,cc_v)
def Cut_Extrude_Column(model,part,ep_w,c_h,eph_cc_b,bh_d,cft_t,cc_v,cc_v_1):
    p = mdb.models[model].parts[part]
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    f, e = p.faces, p.edges
    f1, e1 = p.faces, p.edges
    p.HoleBlindFromEdges(plane=f[167], edge1=e[163], edge2=e[81], planeSide=SIDE1, 
        diameter=bh_d, distance1=c_h/2, distance2=(ep_w-eph_cc_b)/2, depth=cft_t+10)
    p.HoleBlindFromEdges(plane=f1[173], edge1=e1[179], edge2=e1[98], 
        planeSide=SIDE1, diameter=bh_d, distance1=c_h/2, distance2=eph_cc_b+(ep_w-eph_cc_b)/2, 
        depth=cft_t+10)
    p.HoleBlindFromEdges(plane=f[179], edge1=e[197], edge2=e[119], planeSide=SIDE1, 
        diameter=bh_d, distance1=c_h/2-cc_v_1, distance2=(ep_w-eph_cc_b)/2, depth=cft_t+10)
    p.HoleBlindFromEdges(plane=f1[183], edge1=e1[211], edge2=e1[133], 
        planeSide=SIDE1, diameter=bh_d, distance1=c_h/2-cc_v_1, distance2=eph_cc_b+(ep_w-eph_cc_b)/2, 
        depth=cft_t+10)
    p.HoleBlindFromEdges(plane=f[188], edge1=e[226], edge2=e[151], planeSide=SIDE1, 
        diameter=bh_d, distance1=c_h/2+cc_v, distance2=(ep_w-eph_cc_b)/2, depth=cft_t+10)
    p.HoleBlindFromEdges(plane=f1[193], edge1=e1[243], edge2=e1[169], 
        planeSide=SIDE1, diameter=bh_d, distance1=c_h/2+cc_v, distance2=eph_cc_b+(ep_w-eph_cc_b)/2, 
        depth=cft_t+10)

Cut_Extrude_Column(myString,myPart_1,myEndPlate_W,myC_H,myEndPlate_H_CC_T,myBoltHoleDia,myC_FlangeTop_T,Cc_V,myEP_V_D_second_Row)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def Create_Shell_Beam_Partition(model,part):
    p = mdb.models[model].parts[part]
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#1f ]', ), )
    d1 = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d1[3], faces=pickedFaces)
    pickedFaces = f.getSequenceFromMask(mask=('[#3d0 ]', ), )
    d = p.datums
    p.PartitionFaceByDatumPlane(datumPlane=d[2], faces=pickedFaces)

Create_Shell_Beam_Partition(myString,myPart_2)
#------------------------------------------------------------------------------


def Create_Section(model,crosssection,material):
    mdb.models[model].HomogeneousSolidSection(name=crosssection, material=material, thickness=None)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#Create Section
Create_Section(myString,myCS_1_1,myMaterial_1)
Create_Section(myString,myCS_1_2,myMaterial_2)
Create_Section(myString,myCS_3,myMaterial_3)
Create_Section(myString,myCS_4,myMaterial_4)

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def Create_Rebar_And_Prism_Material(model,materialname,fy,fu,sry,sru,em,pr,density):
    mdb.models[model].Material(name=materialname)
    mdb.models[model].materials[materialname].Density(table=((
        density, ), ))
    mdb.models[model].materials[materialname].Elastic(table=((
        em, pr), ))
    mdb.models[model].materials[materialname].Plastic(
        scaleStress=None, table=((fy, sry), (fu, sru)))

#Create_Rebar_And_Prism_Material(myString,myMaterial_1,myFy,myFu,mySry,mySru,myE,myPoiratio,myDensity)
#Create_Rebar_And_Prism_Material(myString,myMaterial_2,myFy,myFu,mySry,mySru,myE,myPoiratio,myDensity)
#Create_Rebar_And_Prism_Material(myString,myMaterial_3,myFy,myFu,mySry,mySru,myE,myPoiratio,myDensity)
#Create_Rebar_And_Prism_Material(myString,myMaterial_4,myFy,myFu,mySry,mySru,myE,myPoiratio,myDensity)
#Create_Rebar_And_Prism_Material(myString,myMaterial_5,myFy,myFu,mySry,mySru,myE,myPoiratio,myDensity)
#Create_Rebar_And_Prism_Material(myString,myMaterial_6,myFy,myFu,mySry,mySru,myE,myPoiratio,myDensity)
#------------------------------------------------------------------------------

def Create_Shell_CS_Beam(model,cs_name,materials,thick_t):
    mdb.models[model].HomogeneousShellSection(name=cs_name, 
        preIntegrate=OFF, material=materials, thicknessType=UNIFORM, thickness=thick_t, 
        thicknessField='', nodalThicknessField='', idealization=NO_IDEALIZATION, 
        poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
        useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)


Create_Shell_CS_Beam(myString,myCS_2,myMaterial_2,myB_Web_T)
Create_Shell_CS_Beam(myString,myCS_2_1,myMaterial_1,myB_FlangeTop_T)

#------------------------------------------------------------------------------

def Section_Assignment(model,part,set_name,section_name):
    p = mdb.models[model].parts[part]
    c = p.cells[:]
    region = p.Set(cells=c, name=set_name)
    p.SectionAssignment(region=region, sectionName=section_name, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Section Assignment
Section_Assignment(myString,myPart_3,"E Plate",myCS_3)
Section_Assignment(myString,myPart_4,"Bolt",myCS_4)

#------------------------------------------------------------------------------

def Web_Assignment(model,part,set_name,web_section):
    p = mdb.models[model].parts[part]
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#908 ]', ), )
    region = p.Set(faces=faces, name=set_name)
    p.SectionAssignment(region=region, sectionName=web_section, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

Web_Assignment(myString,myPart_2,"Web",myCS_2)
#------------------------------------------------------------------------------
def Flange_Assignment(model,part,set_name,flange_section):
    p = mdb.models[model].parts[part]
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#76f7 ]', ), )
    region = p.Set(faces=faces, name=set_name)
    p.SectionAssignment(region=region, sectionName=flange_section, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

Flange_Assignment(myString,myPart_2,"Flange",myCS_2_1)
#------------------------------------------------------------------------------
def Create_Column_Web_Flange_Assignment(model,part,column_flange, column_web, column_fl_cs,column_web_cs):
    p = mdb.models[model].parts[part]
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#ffffffff #27 ]', ), )
    p.Set(cells=cells, name=column_flange)
    cells = c.getSequenceFromMask(mask=('[#0 #3d8 ]', ), )
    p.Set(cells=cells, name=column_web)
    region = p.sets[column_flange]
    p.SectionAssignment(region=region, sectionName=column_fl_cs, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    region = p.sets[column_web]
    p.SectionAssignment(region=region, sectionName=column_web_cs, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

Create_Column_Web_Flange_Assignment(myString,myPart_1,'Column Flange','Column Web', myCS_1_1,myCS_1_1)
#------------------------------------------------------------------------------
from abaqus import *
from abaqusConstants import *
def Create_Csys(model):
# Create a new datum coordinate system
    csys = mdb.models[model].rootAssembly.DatumCsysByThreePoints(
        coordSysType=CARTESIAN, name='csys-1',
        origin=(0.0, 0.0, 0.0),
        point1=(1.0, 0.0, 0.0),
        point2=(0.0, 1.0, 0.0))

Create_Csys(myString)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

def Assemply(model,part,instance,x,y,z):
    a = mdb.models[model].rootAssembly
    p = mdb.models[model].parts[part]
    a.Instance(name=instance, part=p, dependent=ON)
    p =a.instances[instance]
    p.translate(vector=(x,y,z))


# Assemply
Assemply(myString,myPart_1,myInstance_1,0,0,0)
Assemply(myString,myPart_2,myInstance_2,0,0,0)
Assemply(myString,myPart_3,myInstance_3,0,0,0)
Assemply(myString,myPart_4,myInstance_4,0,0,0)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def Create_Beam_Rotation(model,instance,angles,Cx, Cy, Cz):
    a1 = mdb.models[model].rootAssembly
    a1.rotate(instanceList=(instance, ), axisPoint=(0, 0, 0), 
        axisDirection=(Cx, Cy, Cz), angle=angles)

Create_Beam_Rotation(myString,myInstance_1,90.0,myC_Web_T/3, 0.0,0.0)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def Translet_And_Setup(model,instance,Tx, Ty, Tz):
    a = mdb.models[model].rootAssembly
    a.translate(instanceList=(instance, ), vector=(Tx, Ty, Tz))

Translet_And_Setup(myString,myInstance_1,0.0,myC_H,0.0)
Translet_And_Setup(myString,myInstance_2,0.0,myC_H/2,(myC_Depth/2)+myEndPlate_T)
Translet_And_Setup(myString,myInstance_3,0.0,myC_H/2,(myC_Depth/2))
Translet_And_Setup(myString,myInstance_4,myEndPlate_H_CC_B/2,myC_H/2-Cc_V,myC_Web_H/2)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def Create_Bar_Instance_By_Lenear_Pattern(model,instance,spacing,dx,dy, dz,num):
    a1 = mdb.models[model].rootAssembly
    a1.LinearInstancePattern(instanceList=(instance , ), direction1=(dx,dy, dz), direction2=(0.0, 1.0, 0.0), number1=num, number2=1, spacing1=spacing, spacing2=1.0)


Create_Bar_Instance_By_Lenear_Pattern(myString,myInstance_4,myEndPlate_H_CC_T,-1.0,0.0,0.0, 2)
#Create_Bar_Instance_By_Lenear_Pattern(myString,myInstance_4,myEndPlate_Forth_Row,0.0,1.0,0.0, 2)
#Create_Bar_Instance_By_Lenear_Pattern(myString,"Bolt-lin-2-1",myEndPlate_Forth_Row,0.0,1.0,0.0, 2)
Create_Bar_Instance_By_Lenear_Pattern(myString,myInstance_4,Cc_V,0.0,1.0,0.0, 2)
Create_Bar_Instance_By_Lenear_Pattern(myString,"Bolt-lin-2-1",Cc_V,0.0,1.0,0.0, 2)
Create_Bar_Instance_By_Lenear_Pattern(myString,myInstance_4,myEP_V_D_Third_Row ,0.0,1.0,0.0, 2)
Create_Bar_Instance_By_Lenear_Pattern(myString,"Bolt-lin-2-1",myEP_V_D_Third_Row,0.0,1.0,0.0, 2)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------


session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(
    datumAxes=OFF, datumPlanes=OFF)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
myFieldOutName = "F-Output-1"
myStepName_1 = 'Loading'

def Create_Step(model,step_name,fieldname,ini,max_in,min_in,total_t,pre_step):
    mdb.models[model].StaticStep(name=step_name, 
        previous=pre_step, timePeriod=total_t, maxNumInc=100000, initialInc=ini, minInc=min_in, 
        maxInc=max_in, nlgeom=ON)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step=step_name)
    mdb.models[model].FieldOutputRequest(name=fieldname, 
        createStepName=step_name, variables=('S', 'PE', 'PEEQ', 'U', 'RF', 'CF',  
    'EVOL', 'STATUS'), frequency=1)

Create_Step(myString, myStepName_1,myFieldOutName,0.01, 0.1, 1e-15,1.0,'Initial')
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------


def Contact_Property(model,contactprop,frictionfactor):
    mdb.models[model].ContactProperty(contactprop)
    mdb.models[model].interactionProperties[contactprop].NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)
    mdb.models[model].interactionProperties[contactprop].TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((frictionfactor, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005, elasticSlipStiffness=None)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#Contact Property
Contact_Property(myString,"Intprop-1",0.35)

#------------------------------------------------------------------------------
#For Surface
#------------------------------------------------------------------------------

def Create_Surface(model, part, points, surface_name):
    a = mdb.models[model].rootAssembly
    s1 = a.instances[part].faces
    side1Faces = s1.findAt(*points)  # Use multiple points here
    a.Surface(side1Faces=side1Faces, name=surface_name)

# Example call with multiple points:
Create_Surface(myString, myPart_3, (((-(myEndPlate_W/2-2.0), myC_H/2+(myEndPlate_H/2-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),),((myEndPlate_W/2-2.0, myC_H/2+(myEndPlate_H/2-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),),((-(myEndPlate_W/2-2.0), myC_H/2-(myEndPlate_H/2-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),), ((myEndPlate_W/2-2.0, myC_H/2-(myEndPlate_H/2-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),), ((-(myEndPlate_W/4-2.0), myC_H/2+(myEndPlate_H/2-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),), ((myEndPlate_W/4-2.0, myC_H/2+(myEndPlate_H/2-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),),((-(myEndPlate_W/4-2.0), myC_H/2-(myEndPlate_H/2-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),), ((myEndPlate_W/4-2.0, myC_H/2-(myEndPlate_H/2-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),), ((1.0, myC_H/2+(myEndPlate_H/4-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),), ((-1.0, myC_H/2+(myEndPlate_H/4-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),),((1.0, myC_H/2-(myEndPlate_H/4-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),), ((-1.0, myC_H/2-(myEndPlate_H/4-2.0), myC_Web_H/2+myC_FlangeBotom_T+myEndPlate_T),)), 'EPlate_Surface_For_Beam')
#------------------------------------------------------------------------------
#For surface
#------------------------------------------------------------------------------
def Create_Surface_Set(model, part, points, set_name):
    a = mdb.models[model].rootAssembly
    s1 = a.instances[part].edges  # Access the edges for the part
    # Select edges based on the provided points
    side1Edges = s1.findAt(*points) 
    # Create the surface using the selected edges
    a.Surface(side1Edges=side1Edges, name=set_name)

# Modify your points and call the function to create a surface set
Create_Surface_Set(myString, myPart_2, points=(((0.0, myC_H/2, myC_Web_H/2 + myC_FlangeBotom_T + myEndPlate_T),), ((myB_FlangeTop_W/3, myC_H/2 + myB_Web_H/2 + myB_FlangeTop_T/2, myC_Web_H/2 + myC_FlangeBotom_T + myEndPlate_T),), ((myB_FlangeTop_W/3, myC_H/2 - myB_Web_H/2 - myB_FlangeBotom_T/2, myC_Web_H/2 + myC_FlangeBotom_T + myEndPlate_T),), ((-myB_FlangeTop_W/3, myC_H/2 + myB_Web_H/2 + myB_FlangeTop_T/2, myC_Web_H/2 + myC_FlangeBotom_T + myEndPlate_T),), ((-myB_FlangeTop_W/3, myC_H/2 - myB_Web_H/2 - myB_FlangeBotom_T/2, myC_Web_H/2 + myC_FlangeBotom_T + myEndPlate_T),)), set_name='Beam_surf_EP')
#------------------------------------------------------------------------------
def Create_Tie_EP_To_Beam(model,ep_surf,beam_surf,tie_name):
    a = mdb.models[model].rootAssembly
    region1=a.surfaces[ep_surf]
    region2=a.surfaces[beam_surf]
    mdb.models[model].Tie(name=tie_name, main=region1, 
        secondary=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
        tieRotations=ON, thickness=ON)

Create_Tie_EP_To_Beam(myString,'EPlate_Surface_For_Beam','Beam_surf_EP','EPlate_to_Beam_Tie')
#------------------------------------------------------------------------------
#: The interaction "SElf_Contact" has been created.

def Self_Contact(model,set_name, con_prop):
    a = mdb.models[model].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
    mdb.models[model].ContactStd(name=set_name, 
        createStepName='Initial')
    mdb.models[model].interactions[set_name].includedPairs.setValuesInStep(
        stepName='Initial', useAllstar=ON)
    mdb.models[model].interactions[set_name].contactPropertyAssignments.appendInStep(
        stepName='Initial', assignments=((GLOBAL, SELF, con_prop), ))

Self_Contact(myString,'Self Contact',"Intprop-1")
#------------------------------------------------------------------------------
def Create_Reference_Point(x,y,z,model,setname):
    a = mdb.models[model].rootAssembly
    myRP = a.ReferencePoint(point=(x, y, z))
    r = a.referencePoints
    myRP_Position = r.findAt((x, y, z),)    
    refPoints1=(myRP_Position, )
    a.Set(referencePoints=refPoints1, name=setname)
    return myRP,myRP_Position

#------------------------------------------------------------------------------
# Reference point
myRP1,myRP_Position1 = Create_Reference_Point(0,myC_H/2,myLoad_D+myB_Depth/2,myString,'RP-1')
myRP2,myRP_Position2 = Create_Reference_Point(0,myC_H,0,myString,'RP-2')
myRP3,myRP_Position3 = Create_Reference_Point(0,0,0,myString,'RP-3')
#------------------------------------------------------------------------------
def Create_Edge_Set(model, part, points, set_name):
    a = mdb.models[model].rootAssembly
    s1 = a.instances[part].edges  # Access the edges for the part
    edgeSelection = s1.findAt(*points)  # Select edges based on the provided points
    a.Set(edges=edgeSelection, name=set_name)  # Correct keyword is 'edges'

# Assuming you are defining a single point, make sure it's nested properly as a tuple
Create_Edge_Set(myString, myPart_2, points=(((0.0, myC_H/2, myLoad_D+myC_Web_H/2 + myC_FlangeBotom_T),),((myB_FlangeTop_W/3, myC_H/2+myB_Web_H/2+myB_FlangeTop_T/2, myLoad_D+myC_Web_H/2 + myC_FlangeBotom_T),),((myB_FlangeTop_W/3, myC_H/2-myB_Web_H/2-myB_FlangeBotom_T/2, myLoad_D+myC_Web_H/2 + myC_FlangeBotom_T),),((-myB_FlangeTop_W/3, myC_H/2+myB_Web_H/2+myB_FlangeTop_T/2, myLoad_D+myC_Web_H/2 + myC_FlangeBotom_T),),((-myB_FlangeTop_W/3, myC_H/2-myB_Web_H/2-myB_FlangeBotom_T/2, myLoad_D+myC_Web_H/2 + myC_FlangeBotom_T),),), set_name='Beam_Set_RP-1')

#------------------------------------------------------------------------------
def Create_Edge_Set(model, part, points, set_name):
    a = mdb.models[model].rootAssembly
    s1 = a.instances[part].faces  # Access the edges for the part
    faceSelection = s1.findAt(*points)  # Select edges based on the provided points
    a.Set(faces=faceSelection, name=set_name)  # Correct keyword is 'edges'

Create_Edge_Set(myString, myPart_1, points=(((0.0, 0.0, 0.0),),((0.0, 0.0, myC_Web_H/2+myC_FlangeBotom_T/2),),((0.0, 0.0, -(myC_Web_H/2+myC_FlangeBotom_T/2)),),(((myC_FlangeBottom_W/2-1), 0.0, myC_Web_H/2+myC_FlangeBotom_T/2),),(((myC_FlangeBottom_W/2-1), 0.0, -(myC_Web_H/2+myC_FlangeBotom_T/2)),),((-(myC_FlangeBottom_W/2-1), 0.0, myC_Web_H/2+myC_FlangeBotom_T/2),),((-(myC_FlangeBottom_W/2-1), 0.0, -(myC_Web_H/2+myC_FlangeBotom_T/2)),),), set_name='Beam_Set_RP-3')
Create_Edge_Set(myString, myPart_1, points=(((0.0, myC_H, 0.0),),((0.0, myC_H, myC_Web_H/2+myC_FlangeBotom_T/2),),((0.0, myC_H, -(myC_Web_H/2+myC_FlangeBotom_T/2)),),(((myC_FlangeBottom_W/2-1), myC_H, myC_Web_H/2+myC_FlangeBotom_T/2),),(((myC_FlangeBottom_W/2-1), myC_H, -(myC_Web_H/2+myC_FlangeBotom_T/2)),),((-(myC_FlangeBottom_W/2-1), myC_H, myC_Web_H/2+myC_FlangeBotom_T/2),),((-(myC_FlangeBottom_W/2-1), myC_H, -(myC_Web_H/2+myC_FlangeBotom_T/2)),),), set_name='Beam_Set_RP-2')
#------------------------------------------------------------------------------
#RP to Bolt Rigid Body
#------------------------------------------------------------------------------
def Create_Interaction_Rigid_Column(model,rp_name, set_name,rigidbody_name):
    a = mdb.models[model].rootAssembly
    region1=a.sets[rp_name]
    region2=a.sets[set_name]
    mdb.models[model].RigidBody(name=rigidbody_name, 
        refPointRegion=region1, pinRegion=region2)

Create_Interaction_Rigid_Column(myString,"RP-1",'Beam_Set_RP-1',"Beam to RP-1",)
Create_Interaction_Rigid_Column(myString,"RP-2",'Beam_Set_RP-2',"Column to RP-2",)
Create_Interaction_Rigid_Column(myString,"RP-3",'Beam_Set_RP-3',"Column to RP-3")
#------------------------------------------------------------------------------
myBeamMesh_Size = 40.0
myBeamMeshEdge_Size = 8.0
myBolrMesh_Size = 5.0
myEPEdge_num = 2
myEPMesh_Size = 8.0
myColumnMesh_Size = 40.0
myColumnMeshEdge_Size = 8.0
myColumnEdge_num = 2
#------------------------------------------------------------------------------
#Mesh Beam
def Create_Mesh_Beam(model,part,mesh_size_b,mesh_size_b_edge):
    p = mdb.models[model].parts[part]
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#7fff ]', ), )
    p.setMeshControls(regions=pickedRegions, minTransition=ON)
    p.seedPart(size=mesh_size_b, deviationFactor=0.1, minSizeFactor=0.1)
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#c001c000 #3f ]', ), )
    p.seedEdgeBySize(edges=pickedEdges, size=mesh_size_b_edge, deviationFactor=0.1, 
        constraint=FINER)
    p.generateMesh()

Create_Mesh_Beam(myString,myPart_2,myBeamMesh_Size,myBeamMeshEdge_Size)
#------------------------------------------------------------------------------
def Create_Mesh_Bolt(model,part,bolt_size):
    p = mdb.models[model].parts[part]
    c = p.cells
    pickedRegions = c.getSequenceFromMask(mask=('[#fff ]', ), )
    p.setMeshControls(regions=pickedRegions, elemShape=HEX_DOMINATED, 
        technique=SWEEP, algorithm=MEDIAL_AXIS)
    p.seedPart(size=bolt_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

Create_Mesh_Bolt(myString, myPart_4,myBolrMesh_Size)
#------------------------------------------------------------------------------
def Create_Mesh_EP(model,part,ep_edge_num,ep_size):
    p = mdb.models[model].parts[part]
    e = p.edges
    c = p.cells
    pickedRegions = c.getSequenceFromMask(mask=('[#ffff ]', ), )
    p.setMeshControls(regions=pickedRegions, technique=SWEEP, 
        algorithm=MEDIAL_AXIS)
    pickedEdges = e.getSequenceFromMask(mask=(
        '[#1015545 #81002080 #104004a #c80a281d #dc500451 #38 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=ep_edge_num, constraint=FINER)
    p.seedPart(size=ep_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

Create_Mesh_EP(myString, myPart_3,myEPEdge_num,myEPMesh_Size)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def Create_Mesh_Column(model,part,edge_num,edge_size,mesh_size):
    p = mdb.models[model].parts[part]
    c = p.cells
    pickedRegions = c.getSequenceFromMask(mask=('[#ffffffff #3ff ]', ), )
    p.setMeshControls(regions=pickedRegions, technique=SWEEP, 
        algorithm=MEDIAL_AXIS)
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=(
        '[#8000095 #20800ab #80a5 #4400042a #4aaa0a42 #1040492 #41040811', 
        ' #92514410 #21148004 #1000000 #11020920 #132048 ]', ), )
    p.seedEdgeByNumber(edges=pickedEdges, number=edge_num, constraint=FINER)
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=(
        '[#feffffff #bfff7fff #77ffffff #eeffffff #be01ffff #e000147 #fff98000', 
        ' #ff00005d #10bf600b #c0fffcc0 #5c1009 #ebf491f ]', ), )
    p.seedEdgeBySize(edges=pickedEdges, size=edge_size, deviationFactor=0.1, 
        constraint=FINER)
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()


Create_Mesh_Column(myString, myPart_1,myColumnEdge_num,myColumnMeshEdge_Size,myColumnMesh_Size)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def Create_Amp(model,amp_name,t_1, amp_1,t_2, amp_2):
    mdb.models[model].TabularAmplitude(name=amp_name, 
        timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((t_1, amp_1), (t_2, amp_2)))

Create_Amp(myString,'Constant_Amp_Load',0,1,1,1)
Create_Amp(myString,'Ramp_Amp_Def',0,0,1,1)

#------------------------------------------------------------------------------

def Column_Top_Load(model,rp_name,load_name,load):
    a = mdb.models[model].rootAssembly
    region=a.sets[rp_name]
    mdb.models[model].ConcentratedForce(name=load_name, createStepName='Loading', region=region, cf2=load, amplitude='Constant_Amp_Load', distributionType=UNIFORM, field='', 
    localCsys=None)

Column_Top_Load(myString,'RP-3','Column_Top_Load',20000)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

def Create_Column_Bottom_Fixed(model,rp_name,bc_name,u_2):
    a = mdb.models[model].rootAssembly
    region=a.sets[rp_name]
    mdb.models[model].DisplacementBC(name=bc_name, 
        createStepName='Initial', region=region, u1=SET, u2=u_2, u3=SET, ur1=SET, 
        ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)

Create_Column_Bottom_Fixed(myString,'RP-2','Column_Top_Fixed',SET)
Create_Column_Bottom_Fixed(myString,'RP-3','Column_Bottom_Fixed',SET)
myBeamDisplacement = -300.0
#------------------------------------------------------------------------------
def Create_Beam_Def(model,rp_name,bc_name,def_y):
    a = mdb.models[model].rootAssembly
    region=a.sets[rp_name]
    mdb.models[model].DisplacementBC(name=bc_name, 
        createStepName='Loading', region=region, u1=0.0, u2=def_y, u3=UNSET, 
        ur1=UNSET, ur2=0.0, ur3=0.0, amplitude='Ramp_Amp_Def', fixed=OFF, 
        distributionType=UNIFORM, fieldName='', localCsys=None)

Create_Beam_Def(myString,'RP-1','Beam_Deflection',myBeamDisplacement)

#------------------------------------------------------------------------------

#----------------------------------------------------------------------------
def Create_Job(model,job_name):
    mdb.Job(name=job_name, model=model, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
    memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
    multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)

Create_Job(myString,myJobName)
#----------------------------------------------------------------------------
#------------------------------------------------------------------------------

mdb.saveAs(
    pathName='E:/M.Sc Thesis/Paper writting/Journal Paper-1/Parametric Study/Most Important/ID-26')
#: The model database has been saved to "F:\Civil Engineering Tutorials\Abaqus\Tutorial\Abaqus Column Model with python Scripts\Trial_1.cae".
#mdb.jobs[myString].submit(consistencyChecking=OFF)

#E:\M.Sc Thesis\Paper writting\Journal Paper-1\Parametric Study\Most Important
#------------------------------------------------------------------------------