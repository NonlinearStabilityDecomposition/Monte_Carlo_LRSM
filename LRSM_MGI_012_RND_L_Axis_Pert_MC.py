# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 23:40:26 2024

@author: User
"""
from abaqus import *
from abaqusConstants import *
import regionToolset
import __main__
import section
import regionToolset
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import odbAccess
from operator import add

import numpy as np
import matplotlib.pyplot as plt
from random import randint
import os
from random import choice
import random
import re

def Create_Node_Set_ByBoundingCylinder(model,part,x,y,z,length,radius,set_name):
    p = mdb.models[model].parts[part]
    n = p.nodes
    nodes = n.getByBoundingCylinder((x-0.01,y-0.01,z-0.01),(0.0,0.0,length-0.01),radius-0.01)
    p.Set(nodes=nodes, name=set_name)
    

def Create_Node_Set_ByBoundingCylinder2(model,part,x,y,z,length,radius,set_name):
    p = mdb.models[model].parts[part]
    e = p.elements
    elements = e.getByBoundingCylinder((x-0.01,y-0.01,z-0.01),(0.0,0.0,length-0.01),radius-0.01)
    p.Set(elements=elements, name=set_name)


#--------------------------------------------------------------------------

def CLT(myE11,myE22,myG12,myNu12,myLaminate,myLaminateThickness,myPlyNumber):
    #--------------------------------------------------------------------------
    #
    # This function calculates the ABD Composite Stiffness Matrix Components
    #
    #--------------------------------------------------------------------------

    i=0    
    myThickness = myPlyNumber*myLaminateThickness
    Z = []
    Z.append(-myThickness/2.0)
    
    for j in range(1,myPlyNumber+1,1):
        Z.append(Z[j-1]+myLaminateThickness)
    
    # Calcuate Reduced Stiffnesses
    
    myQ11 = []
    myQ12 = []
    myQ16 = []
    myQ22 = []
    myQ26 = []
    myQ66 = []
    
    myQ11_S = []
    myQ12_S = []
    myQ16_S = []
    myQ22_S = []
    myQ26_S = []
    myQ66_S = []
    
    for i in range(0,1,1):
        myQ11.append((myE11**2)/(myE11-myNu12**2*myE22))
        myQ12.append((myNu12*myE11*myE22)/(myE11-myNu12**2*myE22))
        myQ16.append(0)
        myQ22.append((myE11*myE22)/(myE11-myNu12**2*myE22))
        myQ26.append(0)
        myQ66.append(myG12)
    
    i=0    
    for j in range(0,myPlyNumber,1):
        myQ11_S.append(myQ11[i]*np.cos(myLaminate[j]*np.pi/180.0)**4 + 2*(myQ12[i]+2*myQ66[i])*(np.cos(myLaminate[j]*np.pi/180.0))**2*(np.sin(myLaminate[j]*np.pi/180.0))**2 + myQ22[i]*(np.sin(myLaminate[j]*np.pi/180.0))**4)
        myQ12_S.append(myQ12[i]*((np.cos(myLaminate[j]*np.pi/180.0)**4)+(np.sin(myLaminate[j]*np.pi/180.0)**4))+ (myQ11[i]+myQ22[i]-4*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)**2*np.sin(myLaminate[j]*np.pi/180.0)**2)
        myQ16_S.append((myQ11[i]-myQ12[i]-2*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)**3*np.sin(myLaminate[j]*np.pi/180.0) - (myQ22[i] - myQ12[i] - 2*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)*np.sin(myLaminate[j]*np.pi/180.0)**3)
        myQ22_S.append(myQ11[i]*np.sin(myLaminate[j]*np.pi/180.0)**4 + 2*(myQ12[i]+2*myQ66[i])*(np.cos(myLaminate[j]*np.pi/180.0))**2*(np.sin(myLaminate[j]*np.pi/180.0))**2 + myQ22[i]*(np.cos(myLaminate[j]*np.pi/180.0))**4)
        myQ26_S.append((myQ11[i]-myQ12[i]-2*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)*np.sin(myLaminate[j]*np.pi/180.0)**3 - (myQ22[i] - myQ12[i] - 2*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)**3*np.sin(myLaminate[j]*np.pi/180.0))
        myQ66_S.append((myQ11[i] + myQ22[i] - 2*myQ12[i] - 2*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)**2*np.sin(myLaminate[j]*np.pi/180.0)**2 + myQ66[i]*(np.cos(myLaminate[j]*np.pi/180.0)**4+np.sin(myLaminate[j]*np.pi/180.0)**4))
    
    
    
    # Calcualte A Matrix
    #-------------
            
    A11_v = []
    A12_v = []
    A16_v = []
    A22_v = []
    A26_v = []
    A66_v = []
    i = 1
    for j in range(0,myPlyNumber,1):
        A11_v.append(myQ11_S[j]*(Z[i]-Z[i-1]))
        A12_v.append(myQ12_S[j]*(Z[i]-Z[i-1]))
        A16_v.append(myQ16_S[j]*(Z[i]-Z[i-1]))
        A22_v.append(myQ22_S[j]*(Z[i]-Z[i-1]))
        A26_v.append(myQ26_S[j]*(Z[i]-Z[i-1]))
        A66_v.append(myQ66_S[j]*(Z[i]-Z[i-1]))
        i = i+1
    
    A11 = sum(A11_v)
    A12 = sum(A12_v)
    A16 = sum(A16_v)
    A22 = sum(A22_v)
    A26 = sum(A26_v)
    A66 = sum(A66_v)
    
    # Calcualte B Matrix
    #-------------
            
    B11_v = []
    B12_v = []
    B16_v = []
    B22_v = []
    B26_v = []
    B66_v = []
    i = 1
    for j in range(0,myPlyNumber,1):
        B11_v.append(0.5*myQ11_S[j]*(Z[i]**2-Z[i-1]**2))
        B12_v.append(0.5*myQ12_S[j]*(Z[i]**2-Z[i-1]**2))
        B16_v.append(0.5*myQ16_S[j]*(Z[i]**2-Z[i-1]**2))
        B22_v.append(0.5*myQ22_S[j]*(Z[i]**2-Z[i-1]**2))
        B26_v.append(0.5*myQ26_S[j]*(Z[i]**2-Z[i-1]**2))
        B66_v.append(0.5*myQ66_S[j]*(Z[i]**2-Z[i-1]**2))
        i = i+1
    
    
    
    
    B11 = sum(B11_v)
    B12 = sum(B12_v)
    B16 = sum(B16_v)
    B22 = sum(B22_v)
    B26 = sum(B26_v)
    B66 = sum(B66_v)
    
    # Calcualte D Matrix
    #-------------
            
    D11_v = []
    D12_v = []
    D16_v = []
    D22_v = []
    D26_v = []
    D66_v = []
    i = 1
    for j in range(0,myPlyNumber,1):
        D11_v.append((1/3.0)*myQ11_S[j]*(Z[i]**3-Z[i-1]**3))
        D12_v.append((1/3.0)*myQ12_S[j]*(Z[i]**3-Z[i-1]**3))
        D16_v.append((1/3.0)*myQ16_S[j]*(Z[i]**3-Z[i-1]**3))
        D22_v.append((1/3.0)*myQ22_S[j]*(Z[i]**3-Z[i-1]**3))
        D26_v.append((1/3.0)*myQ26_S[j]*(Z[i]**3-Z[i-1]**3))
        D66_v.append((1/3.0)*myQ66_S[j]*(Z[i]**3-Z[i-1]**3))
        i = i+1
    
    
    D11 = sum(D11_v)
    D12 = sum(D12_v)
    D16 = sum(D16_v)
    D22 = sum(D22_v)
    D26 = sum(D26_v)
    D66 = sum(D66_v)
    return A11,A12,A16,A22,A26,A66,B11,B12,B16,B22,B26,B66,D11,D12,D16,D22,D26,D66    

#------------------------------------------------------------------------------ 
    
def create_GeneralStiffness(model,part,Stiff_Name,A11, A12, A22, A13, A23, A33, B11, B21, B31, D11, B12, B22, B32, D12, D22, B13, B23, B33, D13, D23, D33):
    p = mdb.models[model].parts[part]
    mdb.models[model].GeneralStiffnessSection(name=Stiff_Name, referenceTemperature=None, stiffnessMatrix=(A11, A12, A22, A13, A23, A33, B11, B21, B31, D11, B12, B22, B32, D12, D22, B13, B23, B33, D13, D23, D33), applyThermalStress=0, poissonDefinition=DEFAULT, useDensity=OFF)
          
#------------------------------------------------------------------------------

def Create_Node_Set_ByBoundingSphere(model,part,x1,y1,z1,radius,set_name):
    p = mdb.models[model].parts[part]
    e = p.elements
    elements = e.getByBoundingSphere((x1,y1,z1),radius)
    p.Set(elements=elements, name=set_name)

#------------------------------------------------------------------------------  

def Create_Plot(data_1,data_2,max_data_1,max_data_2,name_figure,name):
    fig, ax = plt.subplots()
    ax.plot(data_1, data_2, 'g--', label=name)
    legend = ax.legend(loc='lower right', shadow=True)
    plt.ylabel('Force [N]')
    plt.xlabel('Displacement [mm]')
    plt.grid(True)
    plt.axis([0.0, max_data_1, 0.0, max_data_2])
    plt.savefig(name_figure+str('.png'))

#------------------------------------------------------------------------------  

def Open_ODB_and_Write_NodeSet_data_to_text(model,step_name,variable_name,set_name,Variable_component):
    # open ODB file - ABAQUS Result file
    odb = session.openOdb(str(model)+'.odb')
    
    # list for the VARIABLE you want to evaluate
    Variable_v = []
    
    # analysis step for your VARIABLE
    lastStep=odb.steps[step_name]
    
    #loop over all increments of the analysis step and save VARIABLE information from each increment
    for x in range(len(lastStep.frames)):
        lastFrame = lastStep.frames[x]
        Variable = lastFrame.fieldOutputs[variable_name]
        center = odb.rootAssembly.nodeSets[set_name]
        centerRForce = Variable.getSubset(region=center)
       
        # loop over the VARIABLE and save component (x,y,z - 0,1,2) to list
        for i in centerRForce.values:
            Variable_vr = [i.data[Variable_component]]
            Variable_v = Variable_v + Variable_vr
    
    # Max value of Variable_v
    Max_Variable = [np.max(Variable_v)] 
    Max_Variable_v = [Max_Variable]
            
    # write VARIABLE - component to text file
    
    np.savetxt(str(variable_name)+'_'+str(myString)+'.txt',Variable_v)
    return Variable_v



#session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('CF', NODAL, ((COMPONENT, 'CF1'), )), ), nodeSets=("PIPE-1.VER_RIGHT", ))
#------------------------------------------------------------------------------

def Write_Variable_to_text(variable,variable_name):
         
    # write VARIABLE - component to text file
    
    np.savetxt(str(variable_name)+'_'+str(myString)+'.txt',variable)    
    
#------------------------------------------------------------------------------      

def CreateJob(model, job_name, cpu):
    # mdb.models[model].fieldOutputRequests['F-Output-1'].setValues(variables=PRESELECT, timeInterval=0),
    # mdb.models[model].historyOutputRequests['H-Output-1'].setValues(variables=PRESELECT, timeInterval=0),
    mdb.Job(
        name=job_name,
        model=model,
        description='',
        type=ANALYSIS,
        atTime=None,
        waitMinutes=0,
        waitHours=0,
        queue=None,
        memory=90,
        memoryUnits=PERCENTAGE,

        getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE,
        echoPrint=OFF,      # Deaktiviert echoPrint
        modelPrint=OFF,     # Deaktiviert modelPrint
        contactPrint=OFF,   # Deaktiviert contactPrint
        historyPrint=OFF,   # Deaktiviert historyPrint
        userSubroutine='',
        scratch='',
        resultsFormat=ODB,
        multiprocessingMode=DEFAULT,
        numCpus=cpu,
        numDomains=cpu,
        numGPUs=0
    )
#------------------------------------------------------------------------------  

def SubmitJob(job_name):
    mdb.jobs[job_name].submit()
    mdb.jobs[job_name].waitForCompletion()
    
#------------------------------------------------------------------------------    

def Open_ODB_and_Write_NodeSet_data_to_text2(model,step_name,variable_name,set_name,Variable_component):
    # open ODB file - ABAQUS Result file
    odb = session.openOdb(str(model)+'.odb')
    
    # list for the VARIABLE you want to evaluate
    Variable_v = []
    
    # analysis step for your VARIABLE
    lastStep=odb.steps[step_name]
    
    #loop over all increments of the analysis step and save VARIABLE information from each increment
    for x in range(len(lastStep.frames)):
        lastFrame = lastStep.frames[x]
        Variable = lastFrame.fieldOutputs[variable_name]
        center = odb.rootAssembly.nodeSets[set_name]
        centerRForce = Variable.getSubset(region=center)
       
        # loop over the VARIABLE and save component (x,y,z - 0,1,2) to list
        for i in centerRForce.values:
            Variable_vr = [abs(i.data[Variable_component])]
            Variable_v = Variable_v + Variable_vr
    
    # Max value of Variable_v
    Max_Variable = [np.max(Variable_v)] 
    Max_Variable_v = [Max_Variable]
            
    # write VARIABLE - component to text file
    
    np.savetxt(str(variable_name)+'_'+str(myString)+'.txt',Variable_v)
    return Max_Variable

#------------------------------------------------------------------------------

def DeleteSTT(job_name):
    """
    Deletes the .stt file generated by ABAQUS for the specified job.
    
    Parameters:
        job_name (str): Name of the ABAQUS job.
    """
    stt_file = "{}.stt".format(job_name)  # Klassische String-Formatierung
    if os.path.exists(stt_file):
        try:
            os.remove(stt_file)
            print("Deleted: {}".format(stt_file))
        except Exception as e:
            print("Error deleting {}: {}".format(stt_file, e))
    else:
        print("No .stt file found for job: {}".format(job_name))
        
#------------------------------------------------------------------------------


        
#------------------------------------------------------------------------------
   
# variables

#------------------------------------------------------------------------------



myName = ['IW1_17_t1_LR3']


kkk  =250

Pert_START = [0]
Pert_END   = [kkk]

myAngle2 = []
for i in range(0,kkk,1):
    myAngle2.append(0+i*360.0/kkk)

myRadius = 33.0

myAngle = -2.0
PiDef = np.pi
Limit = 1

# material parameters

myE11 = 208000
myE22 = 208000
myNu12 = 0.3
myG12 = 80000
myG13 = 80000
myG23 = 80000


myLaminate = [45,-45,0,90,90,0,-45,45]

myPlyNumber = len(myLaminate)
myThickness = 0.1
Limit = 1


#pert = 0.0866*PiDef/180.0






# LR dependent data

#myLength = 13.15 # Z50
myLength = 18.59 # Z100
#myLength = 26.3 # Z200
#myLength = 37.19 # Z400

#myLength = 49.2 # Z700

#myLength = 58.81 # Z1000

#myLength = 72.03 # Z1500
#myLength = 83.17 # Z2000
#myLength = 101.87 # Z3000
#myLength = 131.51 # Z5000
#myLength = 185.99 # Z10000
myRadius = 33.0
myThickness = 0.1



###############################################################################
# Gegebene pert-Ranges
pert_ranges = {
    10000: (0.16, 0.32),
    5000: (0.15, 0.22),
    3000: (0.14, 0.20),
    2000: (0.13, 0.18),
    1500: (0.12, 0.17),
    1000: (0.11, 0.16),
    700: (0.09, 0.15),
    400: (0.08, 0.14),
    200: (0.08, 0.13),
    100: (0.07, 0.12),
     50: (0.06, 0.11),
}

def calculate_Z(L, R, t):
    """
    Berechnet Z basierend auf Länge (L), Radius (R) und Dicke (t).
    """
    return L**2 * (1 - 0.3**2)**0.5 / (R * t)

def generate_pert_for_Z(Z, num_samples):
    """
    Generiert Zufallswerte für pert basierend auf Z und den definierten Bereichen.
    """
    # Bereiche für Z durchlaufen
    for key in sorted(pert_ranges.keys(), reverse=True):
        if Z >= key:
            a, b = pert_ranges[key]
            # Gleichverteilte Zufallswerte erzeugen
            return [random.uniform(a, b) for _ in range(num_samples)]
    raise ValueError("Kein passender Bereich für Z={} gefunden.".format(Z))




# Z berechnen
Z = calculate_Z(myLength, myRadius, myThickness)

Z = 100
print("Berechneter Z-Wert: {:.2f}".format(Z))

# Zufallswerte für pert generieren
num_samples = kkk  # Anzahl der Zufallswerte
try:
    pert = generate_pert_for_Z(Z, num_samples)
    print("Generierte pert-Werte: {}".format(pert[:10]))  # Zeige die ersten 10 Werte
except ValueError as e:
    print(e)

###############################################################################


ucl = myLength*myThickness/(myRadius*(np.sqrt(3*(1-0.3**2))))

disp = ucl*4
disp = 1



myPart = 'PART-1'

cpu = 8
#mesh_S = 0.92
#Height_Axial_Stiffener = 3.0

LRSM_Factor = 1000





Load = 8000



# #Liste der verfügbaren Modelle
available_models = [
    # 'IW1_17_scaled',
    # 'IW1_18_scaled',
    # 'IW1_19_scaled',
    # 'IW1_20_scaled',
    # 'IW1_21_scaled',
    # 'IW1_22_scaled',
    # 'IW1_23_scaled',
    # 'IW1_24_scaled',
    # 'IW1_25_scaled',
    # 'IW1_26_scaled',
    # 'IW1_27_scaled',
    # 'IW1_28_scaled',
    # 'IW1_29_scaled',
    # 'IW1_30_scaled',
    # 'IW1_31_scaled',
    # 'IW1_32_scaled',
    # 'IW1_33_scaled',
    # 'IW1_34_scaled',
    # 'IW1_35_scaled',
    # 'IW1_36_scaled',
    # 'IW1_37_scaled',
    # 'IW1_38_scaled',
    # 'IW1_39_scaled',
    # 'IW1_40_scaled',
    # 'IW1_41_scaled',
    # 'IW1_42_scaled',
    # 'IW1_43_scaled',
    # 'IW1_44_scaled',
    # 'IW1_45_scaled',
    # 'IW1_46_scaled',
    # 'IW1_47_scaled',
    'IW1_48_scaled',
    'IW1_49_scaled',
    'IW1_50_scaled',
    'IW1_51_scaled',
    'IW1_52_scaled',
    'IW1_53_scaled',
    'IW1_54_scaled',
    'IW1_55_scaled',
    'IW1_56_scaled',
    'IW1_57_scaled',
    'IW1_58_scaled',
    'IW1_59_scaled',
    'IW1_60_scaled',
   
]


# available_models = [
#     'Z1_1_scaled',
#     'Z1_2_scaled',
#     'Z1_3_scaled',
#     'Z1_4_scaled',
#     'Z1_5_scaled',
#     'Z1_6_scaled',
#     'Z1_7_scaled',
#     'Z1_8_scaled',
#     'Z1_9_scaled',
#     'Z1_10_scaled',
#     'Z1_11_scaled',
#     'Z1_12_scaled',
#     'Z1_13_scaled',
#     'Z1_14_scaled',
#     #'Z1_15_scaled',
#     'Z1_16_scaled',
#     'Z1_17_scaled',
#     'Z1_18_scaled',
#     'Z1_19_scaled',
#     'Z1_20_scaled',
#     'Z1_21_scaled',
#     'Z1_22_scaled',
#     'Z1_23_scaled',
#     'Z1_24_scaled',
#     'Z1_25_scaled',
#     'Z1_26_scaled',
#     'Z1_27_scaled',
#     'Z1_28_scaled'
    

   
# ]

# Anzahl der Positionen adaptiv wählen basierend auf der Modelllänge und/oder MC-Größe
num_positions = max(kkk // 10, kkk)  # Mindestens 50 Positionen oder 10% von kkk

# Zufällige Verteilung der Höhen entlang der Schale (von einem Ende bis zum anderen)
z_pos = [random.uniform(0, myLength) for _ in range(num_positions)]

# Zufällige Verteilung der Winkel
myAngles = [random.uniform(0, 360) for _ in range(num_positions)]


# Normierte Höhen speichern und Winkel ebenfalls sichern
normed_z_pos = [z / myLength for z in z_pos]
with open("positions_and_angles.txt", "w") as file:
    file.write("Normierte Höhe (z/L), Winkel (Grad)\n")
    for z, angle in zip(normed_z_pos, myAngles):
        file.write("{:.4f}, {:.2f}\n".format(z, angle))

# Äußere Schleife für Parameterkombinationen
N = []
P = []
loop_counter = 1  # Zähler für die Loops

for jc in range(0, Limit, 1):
    somefloats_f = []
    somefloats_f2 = []
    somefloats_f3 = []
    pert_L_v = []
    for ic in range(Pert_START[0], Pert_END[0], 1):
        # Zufällige Modellauswahl für jede Simulation
        model_type = choice(available_models)
        print("Verwende Modell: {}".format(model_type))  # Ändere auf .format()
        model_number = re.search(r'\d+', model_type).group()  # Extrahiert die erste Zahl
        # Speichern des Modells in eine Textdatei
        with open("selected_models.txt", "a") as file:
            file.write("'{}'\n".format(model_type))
        myString = 'RSN_Loop_{}'.format(loop_counter)
        loop_counter += 1
        myModel = mdb.Model(name=myString, objectToCopy=mdb.models[model_type])
        LRSM_Radius = pert[ic]
        LRSM_Radius = LRSM_Radius * myRadius
        myAngle = myAngle2[ic]
        Create_Node_Set_ByBoundingSphere(myString,"PART-1",myRadius * np.cos(myAngles[ic] * np.pi / 180.0),myRadius * np.sin(myAngles[ic] * np.pi / 180.0),z_pos[ic], LRSM_Radius,"LRSM")
        a = mdb.models[myString].parts[myPart]
        a.SetByBoolean(name='Set-Shell-LRSM', operation=DIFFERENCE, sets=(a.sets['ESHELL'], a.sets['LRSM'], ))
        p = mdb.models[myString].parts[myPart]
        region = p.sets['Set-Shell-LRSM']
        mdb.models[myString].parts[myPart].sectionAssignments[0].setValues(region=region)
        A11, A12, A13, A22, A23, A33, B11, B12, B13, B22, B23, B33, D11, D12, D13, D22, D23, D33 = CLT(myE11, myE22, myG12, myNu12, myLaminate, myThickness / myPlyNumber, myPlyNumber)
        create_GeneralStiffness(myString, myPart, 'Section-Reference', A11, A12, A22, A13, A23, A33, B11, B12, B13, D11, B12, B22, B23, D12, D22, B13, B23, B33, D13, D23, D33)
        create_GeneralStiffness(myString, myPart, 'Section-Reduced_Stiffness_LRSM', A11 / LRSM_Factor, A12 / LRSM_Factor, A22 / LRSM_Factor, A13 / LRSM_Factor, A23 / LRSM_Factor, A33 / LRSM_Factor, 0, 0, 0, D11, 0, 0, 0, D12, D22, 0, 0, 0, D13, D23, D33)
        region = p.sets['LRSM']
        p.SectionAssignment(region=region, sectionName='Section-Reduced_Stiffness_LRSM', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
        # Aktualisierung der Randbedingungen und Job-Erstellung
        mdb.models[myString].boundaryConditions['Disp-BC-5'].setValues(u3=disp)
        CreateJob(myString, myString, cpu)
        SubmitJob(myString)
        DeleteSTT(myString)
        N.append(Open_ODB_and_Write_NodeSet_data_to_text2(myString, "Step-1", "RF", "RP-1", 2))
        P.append(LRSM_Radius / myRadius)
        Write_Variable_to_text(N, "Buckling Load")
        Write_Variable_to_text(P, "LRSM_Radius_to_Shell_Radius_ratio")


           

