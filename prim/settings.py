#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 17:48:17 2023

@author: ida
"""
from numpy import arange

class Settings:
    """
    The Settings class collects setting properties for the retrievals. 
    
    Args:
        iso_id (lst): list of isotopologue ids, format [(name1, id1),(name2, id2) ...](see hp.gethelp(hp.ISO_ID))
        wave_edges (tuple): Tuple of the wave_edges in format [wave_initial, wave_final]
        sunpath(Tuple): A tuple containing two elements. The first element is a string
                             specifying the data source (e.g., 'TSIS1HSRS'), and the second element
                             is the file path to the solar data.
        
    Optional: 
        dzlay (int): integer how thin the atmospheric layer should be in meter. Default is 1E3 meter. 
        nlay (int): Number of layers in the atmosphere. Default is 20. 
        dwave_lbl (int): Spectral resolution of the line by line grid. Given in nm. Default is 0.01 nm. 
        
        
    Attributes:
        atm_dict (dict): A dictionary containing atmospheric information. 
        iso_id (lst): lst of isotopologue ids
        species (lst): Lst of species that are used. 
        wave_lbl (array-like): Line by line wavelength grid 
        nlay (int): The number of atmospheric layers defined in the settings.
        wave_edges(Tuple): Tuple of wave edges [wave_initial, wave_final]
        wave_alb(lst): Lst of normalized wavelength grid used for the albedo calculation in the physics-based retrievals.
        sunpath(Tuple): A tuple containing two elements. The first element is a string
                             specifying the data source (e.g., 'TSIS1HSRS'), and the second element
                             is the file path to the solar data.
    
    CONTAINS
    
    method __init__(self,iso_id,wave_edges,sunpath,dzlay=1E3,nlay=20,nu=0,dwave_lbl=0.01)
    
    Example: for ch4
        
        wave_edges = [2110, 2450]
        iso_ids = [('CH4', [32, 33, 34, 35]), ('H2O', [1, 2, 3, 4, 5, 6, 129])]
        spectra = ('TSIS1HSRS','../../data/hybrid_reference_spectrum_c2021_03_04_with_unc.nc')
        setting = prim.Settings(iso_ids, wave_edges, spectra) 
    
    
    """
    def __init__(self, iso_id, wave_edges, sunpath, dzlay=1E3, nlay=20, dwave_lbl=0.01):
    
        zlay = (arange(nlay-1, -1, -1) + 0.5) * dzlay  #altitude of layer midpoint
        zlev = arange(nlay, -1, -1) * dzlay        #altitude of layer interfaces = levels
        self.atm_dict = {}
        self.atm_dict['nlay']  = nlay
        self.atm_dict['nlev']  = nlay+1
        self.atm_dict['zlay']  = zlay
        self.atm_dict['zlev']  = zlev
     
        self.iso_id = iso_id
        self.species = [iso[0] for iso in iso_id]
        self.wave_edges = wave_edges
        wave_extend = 20 #nm
        self.wave_lbl = arange(wave_edges[0] - wave_extend, wave_edges[1] + wave_extend,dwave_lbl) #nm
        dwv = self.wave_lbl[-1] - self.wave_lbl[0]
        self.wave_alb = (self.wave_lbl - self.wave_lbl[0])/dwv
        self.sunpath = sunpath


        
