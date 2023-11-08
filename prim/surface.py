# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 20:47:05 2023

@author: ida
"""
import os

import numpy as np

from prim.helper import StopExecution

class Surface:
    """
    The surface_prop class collects methods to
    calculate surface properties
    
    CONTAINS
    
    method __init__(self,wave)
    method get_albedo_flat(alb_prior)
    method get_albedo_CUSTOM(filename,sfname)
    method get_albedo_ECOSTRESS(filename)
    """
    ###########################################################
    def __init__(self, settings):
        """
        init class
        
        Parameters 
    	----------
        settings (object): An object containing various settings and configuration data. Here: array of wavelengths [wavelength] [nm]
        
        Attributes
        ----------
            wave (array): The wavelength grid of the wavlength for the albedo
            alb (array): A zeros array with lenght of wave. 
            
        """
        self.wave = settings.wave_alb
        self.alb = np.zeros_like(self.wave)
        
    ###########################################################
    def get_albedo_flat(self, alb_prior):
        """
        Generate spectrally flat albedo array
           
        Parameters 
    	----------
        alb_prior: albedo value to be used throughout spectral range
        Returns
    	------- 
        alb: constant albedo array [wavelength]
        """    
        #check whether input is in range
        while True:
            if 0. <= alb_prior <= 1.:
               break
            else:
               print("ERROR! surface_prop.get_albedo_flat: albedo prior needs to be in [0...1].")
               raise StopExecution
    
        # Read data from file
        self.alb = np.array([alb_prior for wi in self.wave])
        
    
    ###########################################################
    def get_albedo_poly(self, c_lst):
        """
         Parameters
        ----------
        c_lst: list of constants

        Returns
        -------
        albedo polynomial dependent on wavelenght

        """   
        albedo = np.zeros(len(self.wave))  
        for i in range(0,len(c_lst)):
            albedo=albedo + c_lst[i]*(self.wave)**(i) 
        self.alb=albedo
        return albedo
        
    ###########################################################
    def get_albedo_CUSTOM(self,filename,sfname):
        """
        Read albedo from custom database. This is generic typical 
        data. For a comprehensive albedo database, have a look 
        at: https://speclib.jpl.nasa.gov/
           
        Parameters 
    	----------
        filename: str
        	file with albedo database
        sfname: str
        	name of surface type
        	[sand,soil,snow,vegetation,water]
        Returns
    	------- 
        alb: albedo array interpolated to wavelength [wavelength]
        """    
        #check whether input is in range
        sftypes=['sand','soil','snow','vegetation','water']
        while True:
            if os.path.exists(filename) and sfname in sftypes:
               break
            else:
               print("ERROR! surface_prop.get_albedo_CUSTOM: input out of range.")
               raise StopExecution
    
        # Read data from file
        raw = np.genfromtxt(filename,skip_header=15) # read file into numpy array, skip 15 header lines
        
        # Index of surface type
        isf = sftypes.index(sfname)
        
        # Interpolate albedo to wavelength array
        wave_org = raw[:,0]
        alb_org = raw[:,isf+1]
        self.alb = np.interp(self.wave,wave_org,alb_org)

    ###########################################################
    def get_albedo_ECOSTRESS(self,filename):
        """
        Read albedo from ECOSTRESS database
        at: https://speclib.jpl.nasa.gov/
        
        Parameters 
    	----------
        filename: file with albedo database
        Returns
    	-------  
        alb: albedo array interpolated to wavelength [wavelength]
        """ 
        # Check whether input is in range
        while True:
            if os.path.exists(filename):
               break
            else:
               print("ERROR! surface_prop.get_albedo_ECOSTRESS: input out of range.")
               raise StopExecution
            
        # Read data from file
        raw = np.genfromtxt(filename,skip_header=21,unpack=True) 
    
        wv_in = np.array([a*1E3 for a in raw[0,:]]) # wavelength [nm]
        alb_in = np.array([a/1E2 for a in raw[1,:]]) # albedo [0...1]
        # Check if wavelength in ascending order. If not, flip arrays.
        if wv_in[0] > wv_in[-1]:
           # Interpolate albedo to wavelength array
           self.alb = np.interp(self.wave,np.flip(wv_in),np.flip(alb_in))
        else:
           self.alb = np.interp(self.wave,wv_in,alb_in)
