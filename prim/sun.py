# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 19:40:51 2023

@author: idajandl@gmx.de
"""
import os

import netCDF4 as nc
from numpy import interp, genfromtxt

from prim.helper import StopExecution

class Sun:
    def __init__(self, settings, wave_meas):
        file = settings.sunpath
        if file[0] =="TSIS1HSRS":
            sun = self.read_sun_spectrum_TSIS1HSRS(file[1])
        elif file[0] == "S5P":
            sun = self.read_sun_spectrum_S5P(file[1])
        self.lbl, self.meas = self.interpolate_sun(sun, settings.wave_lbl, wave_meas)
    
    def read_sun_spectrum_S5P(self, filepath):
        """
        Read sun spectrum Sentinel-5 Precursor (S5P) format
        
        Parameters 
        ----------
        filename: str
        	filepath to solar spectrum
        	
        Returns
        ------- 	
        dictionary with wavelength [nm], irradiance [mW nm-1 m-2], irradiance [ph s-1 cm-2 nm-1]
        """
        # check whether input is in range
        while True:
            if os.path.exists(filepath):
               break
            else:
               print("ERROR! read_spectrum_S5P: filename does not exist.")
               raise StopExecution
        
        # Read data from file
        dataset_sun = genfromtxt(filepath,skip_header=42,unpack=True)
    
        # Write sun spectrum in dictionary
        sun = {};
        sun['wl'] = dataset_sun[0,:]
        sun['Wm2nm'] = dataset_sun[1,:]*10**(-3)
        sun['phscm2nm'] = dataset_sun[2,:]
    
        return sun
    
    def read_sun_spectrum_TSIS1HSRS(self, filepath):
        """
        Read sun spectrum TSIS-1 HSRS, downloaded from 
        https://lasp.colorado.edu/lisird/data/tsis1_hsrs, 
        Coddington et al., GRL, 2021, https://doi.org/10.1029/2020GL091709
        NETCDF format: 'pip install netCDF4'
        
        Parameters 
        ----------
        filename: str
        	filepath to solar spectrum
        
        Returns
        ------- 
        dictionary with wavelength [nm], irradiance [W m-2 nm-1]
        """    
        # check whether input is in range
        while True:
            if os.path.exists(filepath):
               break
            else:
               print("ERROR! read_spectrum_TSIS1HSRS: filename does not exist.")
               raise StopExecution
        
        # Open netcdf file
        dataset_sun = nc.Dataset(filepath)
        #print(ds.variables)
        
        # Write sun spectrum in dictionary
        sun = {};
        sun['Wm2nm'] = dataset_sun['SSI'][:] # Solar spectral irradiance [W m-2 nm-1]
        sun['wl'] = dataset_sun['Vacuum Wavelength'][:] # Vacuum Wavelength [nm]
        
        dataset_sun.close
        
        return sun
    
    def interpolate_sun(self, sun, wave_lbl, wave_meas):
        """
        Interpolates the suns measurement to the line-by-line wavlength grid and to the PRISMA measurement wavlength grid. 
        
        Parameters 
        ----------
        sun: array
        	 solar spectrum
        wave_lbl: array 
                line-by-line wavlength grid
        wave_meas: array
                PRISMA measurements wavlength grid
        
        Returns
        ------- 
        Interpolated solar spectrum for line-by-line wavlength grid and PRISMA measurement wavelength grid in units of [Wm2um].
        """    
        spectra = sun['Wm2nm']*10**3 #Wm2um
        wave_sun= sun['wl'][:]
        sun_lbl = interp(wave_lbl,wave_sun,spectra)
        sun_meas = interp(wave_meas,wave_sun,spectra) # interpolated sun spectra with satellite data to calculate  
                                                     	  # later the first guess of the albedo 
        return sun_lbl,sun_meas