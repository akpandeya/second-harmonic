#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 21:40:59 2019

@author: apandeya
"""

import numpy as np
import matplotlib.pyplot as plt
from secondHarmonic import secondHarmonic


import tkinter.filedialog
import tkinter as tk




#File Format Inputs
aAxis = 3 # Angle: degree
rwAxis = 5 # Voltage: Volt
r2wAxis = 7 # Voltage: Volt
hAxis = 2 # Unit: Oe
currentAxis = 4 # Unit: A
fixedCurrent = 1.7e-2 # Unit A works if the current in the file is 0 A or useCurrent Flag is set to 1



xShift = 0 # Unit degree
xmin = -4 # Unit degree
xmax = 364 # Unit degree
xRange = xmax-xmin+1


"""
Does the fitting need to be done after calculatind odd and even part of the raw data? Select the right option
"""
oddEvenSecondHarmonic = 1
oddEvenFirstHarmonic = 0

useCurrent = 1

##########################################
#path = "toProcess"
root = tk.Tk()
root.withdraw()
#path = tkinter.filedialog.askdirectory(parent=root, title='Choose an XRR file')
nameFileArray =  tkinter.filedialog.askopenfilenames(title = "Select file",filetypes = (("CSV files","*.csv"),("all files","*.*")))
##########################################

seondHarmonicAnalysisArray = []
for inputFileName in nameFileArray:
#    with open(inputFileName, "r") as inputFile:
#        line = inputFile.readline().split()
#        print(len(line))
    #if inputFileName != path + '/FititngParam.csv':
    if inputFileName[-3:] == 'dat':
        delim = '\t'
    elif inputFileName[-3:] == 'csv':
        delim = ','
    fileArray = np.genfromtxt(inputFileName,  skip_header=1, delimiter = delim)
    bExt = round(fileArray[0, hAxis])
    angleArray = np.array(fileArray[:, aAxis])
    if (fileArray[0,currentAxis] == 0) or (useCurrent == 1):
        current = fixedCurrent
    else:
        current = fileArray[0,currentAxis]
    rwArray = np.array(fileArray[:, rwAxis])/current
    r2wArray = np.array(fileArray[:, r2wAxis])/current
    
    
    #fig, ax = plt.subplots()
    
    angleArrayInRange = np.array([])
    rwArrayInRange = np.array([])
    r2wArrayInRange = np.array([])
    xRange = angleArray.max() + 1
    for angle, rw, r2w in zip(angleArray, rwArray, r2wArray):
        if int(angle) > xmin-1 and int(angle) < xRange:
            angleArrayInRange = np.append(angleArrayInRange, angle)
            rwArrayInRange = np.append(rwArrayInRange, rw)
            r2wArrayInRange = np.append(r2wArrayInRange, r2w)
    
    # angleArrayInRange = np.sort((angleArrayInRange + xShift)%xRange)
    # rwArrayInRange = rwArrayInRange[np.argsort((angleArrayInRange + xShift)%xRange)]
    # r2wArrayInRange = r2wArrayInRange[np.argsort((angleArrayInRange + xShift)%xRange)]
    
    # angleArray = angleArray + xShift
    
    seondHarmonicAnalysisArray.append(secondHarmonic(angleArrayInRange, rwArrayInRange, r2wArrayInRange, bExt, inputFileName=inputFileName))
    #seondHarmonicAnalysisArray.append(secondHarmonic(angleArray, rwArray, r2wArray, bExt))

with open(nameFileArray[0]+'FititngParam.csv', "a") as outFile1:
    outFile1.write("Field,Current,Rahe, RaheError, Rphe, RpheError, FieldLike, FieldLikeError, DandThermalLike, DampingLikeError, Phi0, Phi0Error, Phi1, Phi1Error, y0, y0Error\n")
    outFile1.write("Oe,A,\Omega,\Omega,\Omega,\Omega,,,,,/degree,/degree,/degree,/degree,\Omega,\Omega \n")
    
    
for s1 in seondHarmonicAnalysisArray:
    if s1.Bext > 0:
        for s2 in seondHarmonicAnalysisArray:
            if abs(s2.Bext + s1.Bext) < 3:
                fig, ax = plt.subplots(6, 1, constrained_layout = True)
                fig.set_size_inches(5,20)
                oddR2wArray = (s1.rw2 - s2.rw2)/2
                evenR2wArray = (s1.rw2 + s2.rw2)/2
                oddRwArray = (s1.rw1 - s2.rw1)/2
                evenRwArray = (s1.rw1 + s2.rw1)/2
                inputFileName = s1.fileName[:-10] + '_' + str(s1.Bext)
                oddEvenHeader = 'Angle,Rw1Positive, Rw1Negative, Rw2Positive, Rw2Negative, RwOddPart, RwEvenPart Rw2OddPart, Rw2evenPart\n'
                oddEvenHeader = oddEvenHeader + '\degree,\Omega,\Omega,\Omega,\Omega,\Omega,\Omega,\Omega,\Omega \n'
                oddEvenHeader = oddEvenHeader + str(s1.Bext) + ' Oe,'+ str(s1.Bext) + ' Oe,' + str(s1.Bext) + ' Oe,' +str(s1.Bext) + ' Oe,' +str(s1.Bext) +  ' Oe,' +  str(s1.Bext) + ' Oe,' + str(s1.Bext) + ' Oe,' +  str(s1.Bext) + ' Oe,' + str(s1.Bext) + ' Oe,'
                oddEvenFileName = inputFileName + 'Oe_oddAndEvenPart.csv'
                
                oddEvenArray = np.array([s1.angles, s1.rw1,s2.rw1, s1.rw2, s2.rw2, oddRwArray, evenRwArray, oddR2wArray, evenR2wArray]).transpose()
                
                np.savetxt(oddEvenFileName, oddEvenArray , header = oddEvenHeader, comments = '', delimiter = ',')
                if oddEvenSecondHarmonic:
                    s1.rw2 = oddR2wArray
                    ax[0].set_ylabel(r"Odd part of $R^{2\omega}_{xy} (\Omega)$")
                    ax[1].set_ylabel(r"Error in fitting Odd part of $R^{2\omega}_{xy} (\Omega)$")
                    
                else:
                    ax[0].set_ylabel(r"$R^{2\omega}_{xy} (\Omega)$")
                    ax[1].set_ylabel(r"Error in fitting Odd part of $R^{2\omega}_{xy} (\Omega)$")
                    
                if oddEvenFirstHarmonic:
                    s1.rw1 = evenRwArray
                    ax[2].set_ylabel(r"Even part of $R^{\omega}_{xy} (\Omega)$")
                    ax[3].set_ylabel(r"Error in fitting Odd part of $R^{\omega}_{xy} (\Omega)$")
                    
                else:
                    ax[2].set_ylabel(r"$R^{\omega}_{xy} (\Omega)$")
                    ax[3].set_ylabel(r"Error in fitting $R^{\omega}_{xy} (\Omega)$")
                
                s1.fitFirstHarmonic()
                s1.fitSecondHarmonic()
                
                 
                
                fittingHeader = 'Angle,FirstHarmonic, FirstHarmonicError,SecondHarmonic, SecondHarmonicError\n'
                fittingHeader = fittingHeader + '\degree,Ohm, Ohm, Ohm, Ohm\n'
                fittingHeader = fittingHeader + str(s1.Bext) + ' Oe,' + str(s1.Bext) + ' Oe,' + str(s1.Bext) + ' Oe,' +str(s1.Bext) + ' Oe,' + str(s1.Bext) + ' Oe,'
                fittingFileName = inputFileName + 'Oe_fitting.csv'
                fittingArray = np.array([ s1.anglesArrayDeg, s1.simrwArray, s1.rwError, s1.simr2wArray, s1.r2wError]).transpose()
                
                np.savetxt(fittingFileName,fittingArray , header = fittingHeader, comments = '', delimiter = ',')
                
                
                
                #Plot Odd Part of Second Harmonic With Fitting
                
                
                ax[0].plot(s1.angles*180/np.pi, oddR2wArray,"or", label =   str(current) + ' A,' + str(s1.Bext) + " Oe")
                ax[0].plot(s1.anglesArrayDeg, s1.simr2wArray, label = r"$R_{fl}=$" + "%.2e $\pm$ %.2e \n" % (s1.flCoeff, s1.flCoeffErr) + 
                  r"$ R_{dl}=$" + "%.2e $\pm$ %.2e " % (s1.adCoeff, s1.adCoeffErr))
                #ax[0].plot(s1.anglesArrayDeg, s1.simr2wArrayPError)
                #ax[0].plot(s1.anglesArrayDeg, s1.simr2wArrayNError)
                
                ax[0].set_xlabel(r"$Angle (\degree)$")
                
                ax[1].plot(s1.angles*180/np.pi, oddR2wArray - s1.secondHarmonicfun(s1.angles, s1.adCoeff, s1.flCoeff, s1.y0, s1.phi1),"or", label =   "Error in fitting second harmonic")
                
                ax[2].plot(s1.angles*180/np.pi, s1.rw1,"or", label =   str(current) + ' A,' + str(s1.Bext) + " Oe")
                ax[2].plot(s1.anglesArrayDeg, s1.simrwArray, label =   r"$R_{AHE}=$" + "%.2e $\pm$ %.2e \n" %(s1.Rahe, s1.RaheErr )+ 
                  r"$ R_{PHE}=$" + "%.2e $\pm$ %.2e " % (s1.Rphe, s1.RpheErr))
                
                ax[2].set_xlabel(r"$Angle (\degree)$")
                
                ax[3].plot(s1.angles*180/np.pi, s1.rw1 - s1.firstHarmonicfun(s1.angles, s1.Rahe, s1.Rphe, s1.phi0),"or", label =  "Error in fitting first harmonic")
                
                ax[4].plot(s1.angles*180/np.pi, evenR2wArray,"or", label =   str(current) + ' A,' + str(s1.Bext) + " Oe")
                ax[4].set_ylabel(r"Even part of $R^{2\omega}_{xy} (\Omega)$")
                ax[4].set_xlabel(r"$Angle (\degree)$")
                
                ax[5].plot(s1.angles*180/np.pi, oddRwArray,"or", label =   str(current) + ' A,' + str(s1.Bext) + " Oe")
                ax[5].set_ylabel(r"Even part of $R^{\omega}_{xy} (\Omega)$")
                ax[5].set_xlabel(r"$Angle (\degree)$")
                
                
                
                
                for graph in ax:
                    graph.legend()
                    graph.tick_params(direction = "in", bottom = True, left = True, right = True, top=True)
                    graph.ticklabel_format(style = "sci", scilimits = (0,0), axis="y", useMathText = True)
                
                fig.savefig(inputFileName + ".pdf", dpi = 600, bbox_inches="tight")
                
                #plt.show()
                plt.close()
                
                with open(nameFileArray[0]+'FititngParam.csv', "a") as outFile1:
                    outFile1.write(str(s1.Bext) + ',' + str(current) + ',' + 
                                   str(s1.Rahe) + ','  + str(s1.RaheErr)+ ',' +
                                   str(s1.Rphe) + ',' + str(s1.RpheErr)+ ',' +
                                   str(s1.flCoeff) + ','+ str(s1.flCoeffErr)+ ',' +
                                   str(s1.adCoeff) + ','+ str(s1.adCoeffErr)+ ',' +
                                   str(s1.phi0*180/np.pi) + ',' + str(s1.phi0Err*180/np.pi)+ ',' +
                                   str(s1.phi1*180/np.pi) + ',' + str(s1.phi1Err*180/np.pi)+ ',' + 
                                   str(s1.y0) + ',' + str(s1.y0Err) + '\n')                    


