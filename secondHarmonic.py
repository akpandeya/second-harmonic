# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 15:07:54 2019

@author: apandeya
"""
import numpy as np
from scipy.optimize import curve_fit

class secondHarmonic():
    def __init__(self, angles = [0], rw1 = [0], rw2 = [0], Bext = 100, Rahe = 0.5, Rphe = 0.5, phi0 = 0, adCoeff = -0.0009, flCoeff = 0.0005, y0 = 0, phi1 = 0, current = 0.004, inputFileName="Input FiLe Name"):
        self.angles = np.array(angles)*np.pi/180
        self.rw1 = rw1
        self.rw2 = rw2
        self.Bext = Bext
        
        self.current = current
        
        self.Rahe = Rahe
        self.Rphe = Rphe
        self.phi0 = phi0
        
        self.adCoeff = adCoeff
        self.flCoeff = flCoeff
        self.y0 = y0
        self.phi1 = phi1
        
        self.firstParams = [self.Rahe, self.Rphe, self.phi0]
        self.bounds1 = ([-np.inf, -np.inf, -np.pi],[np.inf, np.inf, np.pi])
        self.secondParams = [self.adCoeff, self.flCoeff, self.y0, self.phi1]
        self.bounds2 = ([-np.inf, -np.inf, -np.inf, 0],[np.inf, np.inf, np.inf, np.pi/2])
        
        self.anglesArray = np.linspace(self.angles[0], self.angles[-1], 1200)
        self.anglesArrayDeg = self.anglesArray*180/np.pi
        
        self.simrwArray = np.array([])
        self.simr2wArray = np.array([])
        self.fileName = inputFileName

    def firstHarmonicfun(self,phi, Rahe, Rphe, phi0):
        return Rahe + Rphe * np.sin(2*(phi+phi0))
    
    def dampingfun(self, ad, phi):
        return ad*np.cos(phi)
    
    def generateAngle(self, start=0, stop=360, degree=True):
        if degree:
            start = start*np.pi/180
            stop = stop*np.pi/180
        if start > stop:
            start, stop = stop, start
        angleArray = np.arange(start, stop, (stop-start)/3600)
        
        return angleArray
        
    
    def dampingGen(self, start = 0, stop = 2*np.pi, degree=False):
        
        self.angleArray = self.generateAngle(start, stop, degree)
        self.dampingArray = np.array([self.dampingfun(self.adCoeff, phi) for phi in self.angleArray])
        
    def flGen(self, start = 0, stop = 2*np.pi, degree=False):
        self.angleArray = self.generateAngle(start, stop, degree)
        self.flArray = np.array([self.flfun( self.flCoeff, phi) for phi in self.angleArray])
        
        
    
    def flfun(self, fl, phi):
        return fl*(np.cos(3*phi) + np.cos(phi))
    
    def secondHarmonicfun(self, phi, ad, fl, y0, phi1):
        phi = phi+self.phi0 + phi1 # - np.pi/9  # Banabir: Phase of second harmonic in radians Remove the hash before pi
        return  y0 + ad*np.cos(phi) + fl*(np.cos(3*phi) + np.cos(phi))
    
    def paramPmSigma(self, param, sigma):
        
        
        for counter, parameterSigma in enumerate(zip(param, sigma)):
            if counter == 0:
                ParamSigma = [[]]
            else:
                ParamSigma.append([])
            pmArray = np.array([parameterSigma[0] + parameterSigma[1], parameterSigma[0] - parameterSigma[1]])
            ParamSigma[counter].append(pmArray)
        return ParamSigma
    
    def combination(self, paramPm):
        return np.array(np.meshgrid(paramPm[0],paramPm[1],paramPm[2],paramPm[3])).T.reshape(-1,4)
        
    def fitFirstHarmonic(self):
        self.fitFirstHarmonicParams, self.firstHarmCov = curve_fit(self.firstHarmonicfun, self.angles, self.rw1, self.firstParams, bounds = self.bounds1, ftol = 1e-10 )
        self.Rahe = self.fitFirstHarmonicParams[0]
        self.Rphe = self.fitFirstHarmonicParams[1]
        self.phi0 = self.fitFirstHarmonicParams[2]
        
        self.RaheErr = np.sqrt(self.firstHarmCov[0,0])
        self.RpheErr = np.sqrt(self.firstHarmCov[1,1])
        self.phi0Err = np.sqrt(self.firstHarmCov[2,2])
        
        self.paramPm = np.array(self.paramPmSigma(self.fitFirstHarmonicParams, self.firstHarmCov.diagonal()))
        combinationArray = combinationArray = np.array(np.meshgrid(self.paramPm[0],self.paramPm[1],self.paramPm[2])).T.reshape(-1,3)
        r1wValueArray = []
        
        for counter, parameter in enumerate(combinationArray):
            r1wValueArray.append([])
            for angle in self.anglesArray:
                r1wValueArray[counter].append(self.firstHarmonicfun(angle, parameter[0], parameter[1], parameter[2]))
        
        self.rwError = np.array(np.std(r1wValueArray, axis=0))*5
        
        self.simrwArray = np.array([])
        
        for angle in self.anglesArray:
            self.simrwArray = np.append(self.simrwArray, self.firstHarmonicfun(angle, self.Rahe, self.Rphe, self.phi0))
        
        self.simrwArrayPError =  self.simrwArray + self.rwError
        self.simrwArrayNError =  self.simrwArray - self.rwError
    
    def fitSecondHarmonic(self):
        self.fitSecondHarmonicParams, self.secondHarmCov = curve_fit(self.secondHarmonicfun, self.angles, self.rw2, self.secondParams, bounds = self.bounds2, ftol = 1e-10, max_nfev= 10000000)
        self.adCoeff = self.fitSecondHarmonicParams[0]
        self.flCoeff = self.fitSecondHarmonicParams[1]
        self.y0 = self.fitSecondHarmonicParams[2]
        self.phi1 = self.fitSecondHarmonicParams[3]
        
        self.adCoeffErr = np.sqrt(self.secondHarmCov[0,0])
        self.flCoeffErr = np.sqrt(self.secondHarmCov[1,1])
        self.y0Err = np.sqrt(self.secondHarmCov[2,2])
        self.phi1Err = np.sqrt(self.secondHarmCov[3,3])
        
        #print (self.secondHarmCov.diagonal())
        
        self.paramPm = np.array(self.paramPmSigma(self.fitSecondHarmonicParams, self.secondHarmCov.diagonal()))
        combinationArray = np.array(np.meshgrid(self.paramPm[0],self.paramPm[1],self.paramPm[2],self.paramPm[3])).T.reshape(-1,4)
        r2wValueArray = []
        
        for counter, parameter in enumerate(combinationArray):
            r2wValueArray.append([])
            for angle in self.anglesArray:
                r2wValueArray[counter].append(self.secondHarmonicfun(angle, parameter[0], parameter[1], parameter[2], parameter[3]))
        
        self.r2wError = np.array(np.std(r2wValueArray, axis=0))*5
        
        
        self.simr2wArray = np.array([])
        
        for angle in self.anglesArray:
            self.simr2wArray = np.append(self.simr2wArray, self.secondHarmonicfun(angle, self.adCoeff, self.flCoeff, self.y0, self.phi1))
            
        self.simr2wArrayPError =  self.simr2wArray + self.r2wError
        self.simr2wArrayNError =  self.simr2wArray - self.r2wError


        
    """    
    #Algorithm
    
    #1 Fit Rw vs Angle to extract Rahe and Rphe at a given Bext
    #2 Fit R2w vs Angle to extrat coefficient of cos phi and cos 3phi + cos phi
    #2a Coefficient of Cos phi dependence has ANE  + AHE contribution and the other one is FL + Oe contribution
    #3  Differentiate Rw vs angle in oop measurement
    #4 
    y0 + Rahe*cos((phi+phi0)*PI/180)*(Bad/Bext) + 
            Rphe*(cos(3*(phi+phi0)*PI/180) + cos((phi+phi0)*PI/180))*(Bfl + Boe)/(Bext+Ba*(cos((phi+phi0+phi1)*PI/180  )))
    """