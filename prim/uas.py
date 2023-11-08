#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 12:49:46 2023

@author: ida
"""
import os

import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl


from prim.helper import slit_conv, StopExecution
from prim.optics import Optics
from prim.sun import Sun
from prim.atmosphere import Atmosphere
from prim.surface import Surface
from prim.prisma import PRISMA

class UnitAbsorption:
    """
    Represents a unit responsible for calculating absorption properties for atmospheric species.

    This class handles the calculation of absorption properties for various atmospheric species
    using data provided by the PRISMA satellite and atmospheric information. It utilizes the
    'Optics' class to manage optical properties, and the 'Atmosphere' and 'Surface' classes to
    gather atmospheric and surface properties, respectively.

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
        

    Example:
        settings_data = Settings()  # Replace with your actual settings object.
        location_data = Location("Sample_Location", 38.558, 54.2019, "aster_file.nc", ('cams', "cams_file.nc"), [585, 70, 280, 45])
        unit_absorption_calculator = UnitAbsorption(settings_data, location_data)
    """
    def __init__(self, settings, location):
        self.location = location
        self.settings = settings
        self.wave_lbl = settings.wave_lbl
        self.nlay = self.settings.atm_dict['nlay']

        self.atm = Atmosphere(settings, location)
        self.surface = Surface(settings)

        self.optics = Optics(settings)
        print(self.optics.species)
        self.measurement = PRISMA(settings, location)

        meas = self.measurement.read_netcdf4(0, 0)
        self.wave_meas = meas["wave"]
        self.fwhm = meas["fwhm"]
        self.sun_lbl = Sun(self.settings, meas['wave']).lbl
        self.pklfile = '../../tmp/optics_prop' + self.location.name + '.pkl'  # TODO:
        if os.path.exists(self.pklfile):
            self.optics.prop = pkl.load(open(self.pklfile, 'rb'))
        else:
            if not os.path.exists("tmp"):
                os.makedirs("tmp")
            self.optics.cal_molec_xsec(self.atm)
            pkl.dump(self.optics.prop, open(self.pklfile, 'wb'))
        self.ref_CO2, self.ref_CH4, self.ref_H2O, self.ref_N2O, self.ref_NH3 = [
            self.atm.atmo['CO2'], self.atm.atmo['CH4'], self.atm.atmo['H2O'], self.atm.atmo['N2O'], self.atm.atmo['NH3']]    
    def transmission(self, mu0, muv):
        """
        Calculate transmission solution given
        geometry (mu0,muv) using matrix algebra and albedo = 100%. 
        
        Parameters 
        ----------
        optics: optic_prop object
        
        mu0: cosine of the solar zenith angle [-]
        muv: cosine of the viewing zenith angle [-]
        
        Returns
        ------- 
        rad_trans: single scattering relative radiance [wavelength] [1/sr]
        """    
        while True:
            if 0. <= mu0 <= 1. and -1. <= muv <= 1.:
               break
            else:
               print("ERROR! transmission: input out of range.")
               raise StopExecution

        
        tautot    = self.optics.prop['tautot'][:]
        mueff     = abs(1./mu0) + abs(1./muv)
        exptot    = np.exp(-tautot*mueff)
        rad_trans = self.sun_lbl*1*exptot 
        return rad_trans
    
    def spectra(self, i, j, ch4=True, co2=False, pertubation=None, plot=False):
        """
        Calculate unit absorption spectrum given
         indecees of the PRISMA data image (for muv and mu0),
        perturbation in ppm for CO2 and ppb for CH4, and plots the result, if necessary.
        
        Parameters 
        ----------
        i (int): An index related to across-track axis in the PRISMA data.
        j (int): Another index related to along-track axis in the PRISMA data.
        ch4 (boolean): if true: CH4 is the molecule of interest
        co2 (boolean): if true: CO2 is the molecule of interest 
        perturbation (int): default: None, means perturbation is 5000 ppb for methane and 50ppm for co2. 
                        if set manually, perturbation for methane is given in ppb and perturbation is given in ppm. 
        plot (boolean): default: False, means no plotting of the unit absorption spectrum. If plot=True,
                        than unit absorption spectrum is plotted. 
    
        
        Returns
        ------- 
        S_conv : unit absorption spectrum [m^2/ppb] [m^2/ppm]
        """    
        self.optics.set_opt_depth_species(self.atm.atmo)
        meas = self.measurement.read_netcdf4(i,j)
        slit = slit_conv(self.fwhm, self.wave_lbl, self.wave_meas)
        mu0 = meas["mu0"]
        muv = meas["muv"]
        if pertubation == None: 
            if ch4:
                pertubation = 1000 #ppb
            
            if co2: 
                pertubation = 50 #ppm
        
        I_0 = self.transmission(mu0, muv)
        I_0_conv = slit.dot(I_0)
        
        mueff     = abs(1./mu0) + abs(1./muv)
        if ch4 and co2: 
            print("WARNING! Methane and Carbon dioxide are both true. Only one is allowed. Set CH4=False or CO2=False")
        if not ch4 and not co2:
            print("WARNING! Methane and Carbon dioxide are both false. Only one is allowed. Set CH4=True or CO2=True")
        if ch4:
            I_pert = I_0*np.exp(-pertubation*1E-9*mueff*self.optics.prop["CH4"]['xsec'][:,self.nlay-1]*self.atm.atmo['AIR'][self.nlay-1]*1.E-4)
        if co2:
            I_pert = I_0*np.exp(-pertubation*1E-6*mueff*self.optics.prop["CO2"]['xsec'][:,self.nlay-1]*self.atm.atmo['AIR'][self.nlay-1]*1.E-4)
        I_pert_conv = slit.dot(I_pert)
        
        
        unit_absorption_conv = np.zeros(len(I_0_conv))
        for i in range(0, len(I_0_conv)):
            if (I_pert_conv[i] != 0 and I_0_conv[i] != 0):
                unit_absorption_conv[i] = 1/(pertubation*self.atm.atmo['AIR'][self.nlay-1]*1E-4) * np.log(I_pert_conv[i]/I_0_conv[i])
        if plot:
            plt.plot(self.wave_meas, unit_absorption_conv, label="convoluted fwhm=9.8")
        return unit_absorption_conv