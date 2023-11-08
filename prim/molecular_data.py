# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 20:44:31 2023

@author: ida
"""
import os

import prim.hapi as hp
from prim.helper import StopExecution

class MolecularData:
    """
    The molecular_data class collects method for calculating
    the absorption cross sections of molecular absorbers
    
    CONTAINS
    
    method __init__(self,wave)
    
    method get_data_HITRAN(self,xsdbpath, hp_ids)
    """
    ###########################################################
    def __init__(self,settings):
        """
        init class
        
        arguments: 
        wave: array of wavelengths [wavelength] [nm]
        xsdb: dictionary with cross section data
        """
        self.xsdb = {}
        self.wave = settings.wave_lbl
        
        self.get_data_HITRAN(settings.iso_id)
    ###########################################################
    def get_data_HITRAN(self,hp_ids):
        """
        Download line parameters from HITRAN web ressource via
        the hapi tools, needs hapy.py in the same directory
        
        Parameters 
    	----------
        xsdbpath: path to location where to store the absorption data
        hp_ids: list of isotopologue ids, format [(name1, id1),(name2, id2) ...](see hp.gethelp(hp.ISO_ID))
        
        Returns
    	-------
        
        xsdb[id][path]: dictionary with paths to HITRAN parameter files
        """
        #check whether input is in range
        while True:
            if len(hp_ids)>0:
               break
            else:
               print("ERROR! molecular_data.get_data_HITRAN: provide at least one species.")
               raise StopExecution
               
        
        xsdbpath = '../../data'
        wv_start = self.wave[0]
        wv_stop = self.wave[-1]
        
        #hp.gethelp(hp.ISO_ID)
        for id in hp_ids:
            key = '%s'%id[0] #string: molecule name (example: "CH4")
            self.xsdb[key]={}
            self.xsdb[key]['isotopologues']=id[1] #list of isotopologue of the desired molecule (example for CH4: [32,33,34] )
            self.xsdb[key]['name']='ID%s_WV%d-%d'%(id[0],wv_start,wv_stop)# write 1 file per molekule (included all isopotopologues that you want)
            # (file)name looks like "IDCH4_WV1300-1400"
            # Check if data files are already inplace, if not: download
            if (not os.path.exists(os.path.join(xsdbpath,self.xsdb[key]['name']+'.data')) and not os.path.exists(os.path.join(xsdbpath,self.xsdb[key]['name']+'.header'))):
                hp.fetch_by_ids(self.xsdb[key]['name'],id[1],1.E7/wv_stop,1.E7/wv_start) # wavelength input is [nm], hapi requires wavenumbers [1/cm]
                    
