#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 12:23:05 2023

@author: ida
"""

import numpy as np
#import matplotlib.pyplot as plt

#from prim.helper import slit_conv, StopExecution,save3d,load3d
from prim.prisma import PRISMA
from prim.uas import UnitAbsorption

class MaFt():
    """
    This class contains the matched filter retrieval. 

    Args:
        settings (object): An object containing various settings and configuration data.
        location (object): An object representing the location information.

    Attributes:
        location (object): An instance of the 'Location' class containing geographic information.
        settings (object): An object containing various settings and configuration data.
        wave_lbl (array-like): An array of wavelength values specified in the settings.
        nlay (int): The number of atmospheric layers defined in the settings.
        atm (Atmosphere): An instance of the 'Atmosphere' class to manage atmospheric data.
        surface (Surface): An instance of the 'Surface' class to manage surface data.
        optics (Optics): An instance of the 'Optics' class to manage optical properties.
        measurement (PRISMA): An instance of the 'PRISMA' class to handle PRISMA data.

    Methods:
        - transmission: Calculate transmission.
        - spectra: Calculate unit absorption spectrum
        """
    def __init__(self, settings, location):
        self.location = location
        self.settings = settings
        self.measurement = PRISMA(settings, location)
        meas = self.measurement.read_netcdf4(0, 0)
        self.wave_meas = meas["wave"]
        self.unit_absorption_spectrum = None
        
    def set_unit_absorption_spectrum(self, unit_absorption_spectrum):
        self.unit_absorption_spectrum = unit_absorption_spectrum
        
    def r_albedo(self, x, mean):
        r=(x.T @ mean)/(mean.T @ mean)
        return r
    
    def a_hat(self, x, mean, Sy, t, r):
        a = ((x - mean).T @ np.linalg.inv(Sy) @ t) / (r * t.T @ np.linalg.inv(Sy) @ t)
        return a
    
    def a_err(self, t, Sy):
        return 1/(np.sqrt(t.T @ np.linalg.pinv(Sy) @ t))
    
    def calculate_covariance_matrix(self, y, m, n, k, l):
        len_column = 1000
        mean_arr = np.zeros([n,len(self.wave_meas)])
        C = np.zeros([n,len(self.wave_meas),len(self.wave_meas)])
        len_background = len_column-l-3
        background = np.zeros([n,len(self.wave_meas),len_background])
        #if destripe:
        #    y=load3d("../data_final/detrended_data/"+location+"spectrum_low_pass_filter_across_along_0.5_MaFt.txt",997)
        for i in range(m,m+n):
            a = np.arange(k,k+l)
            b = np.array([0,1,999])
            index = np.concatenate((a,b), axis=0)
            background[i-m,:,:] = np.delete(y[i-m,:,:],index,1)
            mean_arr[i-m,:] = np.mean(background[i-m,:,:],axis=1) #mean per column len_column
            for j in range(0,len_background):
                #print(i-m)
                # covariance matrix sum of the outer product of each individual 
                # spectrum in a column and the mean of the column
                C[i-m,:,:] = C[i-m,:,:] + np.outer((background[i-m,:,j]-mean_arr[i-m,:]),#j-k
                                               (background[i-m,:,j]-mean_arr[i-m,:]))
            C[i-m,:,:] = 1/(len_background) * C[i-m,:,:] #normalization
        return mean_arr,C,background
    
    def calculate_normalized_covariance_matrix(self, m, n, k, l):
        y = np.zeros([n, len(self.wave_meas), 1000]) #shape: sample, nwave,line 0:1000
        len_column = 1000
        for i in range(m, m+n):
            meas = self.measurement.read_netcdf4(i, 0)
            y[i-m,:,:] = meas["spec_per_column"][:,:] 
            for j in range(0, len_column):
                    y[i-m,:,j] = y[i-m,:,j]/np.max(y[i-m,:,j])
        mean_arr, C, background = self.calculate_covariance_matrix(y,m,n,k,l)
        return C
    
    def background_per_column(self, m, n, k, l):
        y = np.zeros([n,len(self.wave_meas),1000]) #shape: sample, nwave,line 0:1000
        target_spectrum = np.zeros([n, len(self.wave_meas)])
        for i in range(m,m+n):
            meas = self.measurement.read_netcdf4(i,0)
            y[i-m,:,:] = meas["spec_per_column"][:,:] 
        mean_arr,C,background = self.calculate_covariance_matrix(y,m,n,k,l)
        for i in range(m,m+n):
            target_spectrum[i-m,:] = mean_arr[i-m,:]*self.unit_absorption_spectrum
        return mean_arr, C, target_spectrum, background
            
    def run(self, m, n, k, l):
        try:
            self.unit_absorption_spectrum
        except:
            print("No unit absorption spectrum as input, calculation intern")
            self.unit_absorption_spectrum = UnitAbsorption(self.settings,self.location).spectrum(m,k)
        mean, Sy, target_spectrum, background = self.background_per_column(m,n,k,l)
        a_lst = np.zeros((n,l))
        delta_a = np.zeros((n,l))
        for i in range(m,m+n):
            print(i)
            for j in range(k,k+l):
                    #if destripe:
                    #    x=y[i-m,:,j-2]#meas["spec"]/np.max(meas["spec"])
                meas = self.measurement.read_netcdf4(i,j)
                x = meas["spec"]
                r = self.r_albedo(x,mean[i-m,:])
                a_lst[i-m,j-k] = self.a_hat(x,mean[i-m,:],Sy[i-m,:,:],target_spectrum[i-m,:],r)
                delta_a[i-m,j-k] = self.a_err(target_spectrum[i-m,:],Sy[i-m,:,:])
        return a_lst, delta_a
            
        
    