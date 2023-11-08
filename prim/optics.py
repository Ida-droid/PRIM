# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 20:47:03 2023

@author: idajandl@gmx.de
"""

import numpy as np

import prim.hapi as hp
from prim.molecular_data import MolecularData

PSTD = 1013.25            # Standard pressure [hPa]

class Optics:
    """
    The optic_ssc_prop class collects methods to
    calculate optical scattering and absorption
    properties of the single scattering atmosphere
    
    CONTAINS
    
    method __init__(self, wave, zlay)
    method cal_molec_xsec(self, molec_data, atm_data)
    method set_opt_depth_species(self, atm_data,species)
    """
    def __init__(self, settings):
        """
        init class
        Arguments
        ---------
        settings: 
            
        Parameters 
    	----------
        prop: dictionary of contributing phenomena
        prop[wave]: array of wavelengths [wavelength] [nm]
        prop[zlay]: array of vertical height layers, midpoints [nlay] [m]
        """
        self.prop   = {}  
        self.wave   = settings.wave_lbl  
        self.zlay   = settings.atm_dict['zlay']  
        self.species = settings.species
        self.molec  = MolecularData(settings)
        
    def cal_molec_xsec(self, atm_data):
        """
        Calculates molecular absorption cross sections
        
        Parameters 
    	----------
        atm_data: atmosphere_data object
        
        Returns
    	------- 
        prop['molec_XX']: dictionary with optical properties with XXXX HITRAN identifier code
        prop['molec_XX']['xsec']: absorption optical thickness [wavelength, nlay] [-]
        """
        nlay  = self.zlay.size
        nwave = self.wave.size
        nu_samp = 0.005 # Wavenumber sampling [1/cm] of cross sections
        # Loop over all molecules, id = molecule name
        for id in self.molec.xsdb.keys():            
            name = id
            isotopologues = self.molec.xsdb[id]['isotopologues'] #list of all isotopologues from the molecule
 
            # Write absorption optical depth [nwave,nlay] in dictionary / per molecule
            self.prop[name]={}
            self.prop[name]['xsec']     = np.zeros((nwave,nlay))
            self.prop[name]['isotopologues']  = isotopologues
            # Dictionary which looks like {'CH4': {'xsec': array([[0,0,0....]...]), 'isotopologues': [32,33], 'H2O':......}
            # Check whether absorber type is in the atmospheric data structure
            if name not in atm_data.atmo.keys():
               print("WARNING! optic_prop.cal_molec: absorber type not in atmospheric data.",id,isotopologues)
            else:
               # Loop over all atmospheric layers
               for ki,pi,Ti in zip(range(len(atm_data.atmo['zlay'])),atm_data.atmo['play'],atm_data.atmo['tlay']):
                   
                   if name=="H2O":
                       print("yes")
                       vmr_h2o=atm_data.atmo["H2O"][ki]/atm_data.atmo["AIR"][ki]
                       pi_corr=pi*(1+4*vmr_h2o)
                       nu,xs = hp.absorptionCoefficient_Voigt(SourceTables=self.molec.xsdb[id]['name'], Environment={'p':pi_corr/PSTD, 'T':Ti},
                                                              WavenumberStep=nu_samp )#,Diluent={"air":(1-vmr_h2o),"self":vmr_h2o})#,GammaL='gamma_air'
                       	
                       self.prop[name]['xsec'][:,ki] = np.interp(self.wave,np.flip(1E7/nu),np.flip(xs))
                   else:# Calculate absorption cross section for layer
                       #    vmr_ch4=atm_data.atmo["CH4"][ki]/atm_data.atmo["AIR"][ki]
                       nu,xs = hp.absorptionCoefficient_Voigt(SourceTables=self.molec.xsdb[id]['name'], Environment={'p':pi/PSTD, 'T':Ti},
                                                              WavenumberStep=nu_samp)#,Diluent={"air":(1-vmr_ch4),"self":vmr_ch4})
                   	   # Interpolate on wavelength grid provided on input
                       self.prop[name]['xsec'][:,ki] = np.interp(self.wave,np.flip(1E7/nu),np.flip(xs))
                   

    def set_opt_depth_species(self, atm_data):
        """
        calaculates absorption optical depth from the various species and combines those as specified
        
        Parameters 
    	----------
        atmospheric input and species to be combined to total optical depth
        Returns
        	-------
            prop[taua]: total absorption optical thickness array [wavelength, nlay] [-]
            """
            #check whether input is in range
        #print("bin drinnen,wirklich, CH4 ")
        species = self.species
        nlay  = self.zlay.size
        nwave = self.wave.size

        conv = 1.E-4   #cross sections are given in cm^2, atmospheric densities in m^2 
        #correction for water
        """
        for k in zip(range(len(atm_data.atmo["zlay"]))):
            self.prop["H2O"]['xsec'][:,ki]=self.prop["H2O"]['xsec'][:,ki]+
        
        do k = 1, natm!*** tau + N*dtau/dN
        cross_section(1:nwave,k,i) = cross_section(1:nwave,k,i) + &
	          			      vmr_h2o(k)*(cross_section_h2o_per(1:nwave,k,i)-cross_section(1:nwave,k,i))/dvmr
        """
        for name in species:
                self.prop[name]['tau_per_molec']= np.zeros((nwave,nlay))
                for ki in range(nlay):
                    self.prop[name]['tau_per_molec'][:,ki] = self.prop[name]['xsec'][:,ki] * atm_data[name][ki]*conv #Tau = sigma*n
            
        # tau gesamt für jedes nlay
        
        self.prop['tau_per_nlay'] = np.zeros((nwave,nlay))    
        
        for name in species:
            self.prop['tau_per_nlay'] = self.prop['tau_per_nlay'] + self.prop[name]['tau_per_molec']
          
        # Total vertical optical thickness per layer (Delta tau_k) [nwave,nlay]
        tauk = self.prop['tau_per_nlay'] 
        
        # total optical thickness per spectral bin [nwave] 
        self.prop['tautot'] = np.zeros([nwave]) 
        self.prop['tautot'][:] = np.sum(tauk, axis=1)