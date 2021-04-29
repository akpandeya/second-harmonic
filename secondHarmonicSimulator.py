#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 31 12:17:50 2020

@author: apandeya
"""
import os
import numpy as np
from secondHarmonic import secondHarmonic
import matplotlib.pyplot as plt
from matplotlib import cm
##########################################
path = "toProcess"
##########################################

fileExtensionCode = 2  # 1: "dat", 2: "csv"
delimiterCode = 2  # 1: Tab, 2: Comma, 3: Semicolon
numberOfHeaderLines = 2

##########################################
"""
Selecting the correct axes for plotting
"""

fieldAxis = 0
currentAxis = 1
RAheAxis = 2	 
rAheErrorAxis = 3	 
rPheAxis = 4	 
rPheErrorAxis = 5	 
fieldLikeAxis = 6 
fieldLikeErrorAxis = 7	 
dAndThermalLikeAxis = 8	 
dAmpingLikeErrorAxis = 9	 
phi0Axis = 10	 
phi0ErrorAxis = 11	 
phi1 = 12	 
phi1Error = 13

##########################################
if fileExtensionCode == 1:
    fileExtension = ".dat"
    extensionPosition = -4
elif fileExtensionCode == 2:
    fileExtension = ".csv"
    extensionPosition = -4

# Select correct delimiter in the file
    
if delimiterCode == 1:
    delimit = "\t"
elif delimiterCode == 2:
    delimit = ","
elif delimiterCode == 3:
    delimit = ";"

##########################################
    
nameFileArray = tuple(path + '/' + file for file in os.listdir(path) if file.endswith(fileExtension))

##########################################
dampFig, dampAx = plt.subplots()
flFig, flAx = plt.subplots()
for inputFileName in nameFileArray:
    print ("##########################################")
    print (inputFileName + "\n")
    fileArray = np.genfromtxt(inputFileName,  skip_header=numberOfHeaderLines, delimiter = delimit)
    for index, row in enumerate(fileArray):
        colorNumber = np.linspace(0, 1, len(row)) 

        colors = [ cm.ocean(x) for x in colorNumber ]
        secondHarmonicObject = secondHarmonic()
        secondHarmonicObject.Bext = row[fieldAxis]
        secondHarmonicObject.Rahe = row[RAheAxis]
        secondHarmonicObject.Rphe = row[rPheAxis]
        secondHarmonicObject.phi0 = row[phi0Axis]
        secondHarmonicObject.adCoeff = row[dAndThermalLikeAxis]
        secondHarmonicObject.flCoeff = row[fieldLikeAxis]
        secondHarmonicObject.phi1 = row[phi1]
        
        secondHarmonicObject.dampingGen()
        secondHarmonicObject.flGen()
        
        
        
        
        dampAx.plot(secondHarmonicObject.angleArray*180/np.pi, secondHarmonicObject.dampingArray/1e-3, color = colors[index], label = f"{int(secondHarmonicObject.Bext)}")
        flAx.plot(secondHarmonicObject.angleArray*180/np.pi, secondHarmonicObject.flArray/1e-6, color = colors[index], label = f"{int(secondHarmonicObject.Bext)}")
        
        #dampAx.ticklabel_format(style = "sci", scilimits = (0,0), axis="y", useMathText = True)
        
        dampAx.set_xlabel(r"Angle (degree)")
        dampAx.set_ylabel(r"Damping Like Contribution (m$\Omega$)")
        dampAx.legend()
        
        flAx.set_xlabel(r"Angle (degree)")
        flAx.set_ylabel(r"Field Like Contribution ($\mu\Omega$)")
        flAx.legend()
        #flAx.annotate("Field Like Contribution", xy = (100,-20))
        
        #dampFig.savefig( inputFileName + "damping.pdf", bbox_inches="tight", dpi=600)
        #flFig.savefig( inputFileName + "fieldLike.pdf", bbox_inches="tight", dpi=600)

plt.show()
        
        