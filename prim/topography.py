# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 19:35:08 2023

@author: ida
"""
import netCDF4 as nc
import numpy as np

g0   = 9.80665            # Standard gravity at ground [m/s2]
class Topography:
    """
    Calculate surface pressure for a given longitude and latitude using meteorological data from CAMS
    and topography information from ASTER.

    This class calculates the surface pressure for a given location by utilizing meteorological
    data from CAMS (Copernicus Atmosphere Monitoring Service) and topography information from
    ASTER (Advanced Spaceborne Thermal Emission and Reflection Radiometer).

    Args:
        location (Location): An instance of the 'Location' class containing the geographic
            information necessary for calculations.

    Attributes:
        topo (dict): A dictionary containing topography information, including longitude [-180, 180], latitude,
            and longitude within the range [0, 360].
        location (Location): An instance of the 'Location' class representing the location data.
        
    Methods:
        - get_meterological_data_cams: Retrieve meteorological data from CAMS.
        - get_height: Obtain topographic height data.
        - calculate_surf_press: Calculate surface pressure using the acquired data.

    Example:
        location_data = Location("Sample_Location", 38.558, 54.2019, "aster_file.nc", ('cams', "cams_file.nc"), [585, 70, 280, 45])
        topography_calculator = Topography(location_data)
        topography_calculator.calculate_surf_press()
    """
    def __init__(self, location):
        self.topo = {}
        self.topo['lon'] = location.lon % 360
        self.topo['lat'] = location.lat
        self.topo['lon_360'] = location.lon
        self.location = location
        file = location.filename_atm

        self.get_meterological_data_cams(file[1])
        self.get_height()
        self.calculate_surf_press()

    def get_meterological_data_cams(self,filepath):
        """
        read out meterological information from CAMS 
         
        
        Parameters 
    	----------
        filename: str 
        	path to CAMS file
                    
        Returns
    	-------
        topo[height_1]: float 
        	height at lowest atmospheric layer at given latitude and longitude 
            #TODO: do flag if mean height of one picture or only height at one pixel is necessary
        	mean value is better if one chooses one (mean) surface pressure for whole image (because of buildings and small deviations)
        topo[press_1]: float 
            pressure at lowest atmospheric layer at given latitude and longitude
        topo[temp_1]: float
            temperature at lowest atmospheric layer at given latitude and longitude
        """
        ds = nc.Dataset(filepath)
        #print(ds.variables)
    
        # Select month index, latitude/longitude index (next neighbour)
        itime = 0
        ilat = np.argmin(abs(ds['latitude'][:] - self.topo['lat']))
        ilon = np.argmin(abs(ds['longitude'][:] - self.topo['lon_360']))
        
        # ECMWF: Geopotential [m2 s-2] converted to height [m], approximate use of g0 
        zalt_in = np.array([d/g0 for d in ds['z'][itime,:,ilat,ilon]])
        # ECMWF: Pressure [hPa]
        press_in = ds['level'][:]
        # ECMWF: Temperature [K]
        temp_in = ds['t'][itime,:,ilat,ilon]
        
        ds.close
        self.topo['height_1'] = zalt_in[-1]
        self.topo['press_1'] = press_in[-1]
        self.topo['temp_1'] = temp_in[-1]

    def get_heigth(self): 
        """
        read out height information from ASTER 
         
        
        Parameters 
    	----------
        None
                    
        Returns
    	-------
        topo[height_2]: float 
        	height at latitude and longitude #TODO: do flag if mean height of one picture or only height at one pixel is necessary
        	mean value is better if one chooses one (mean) surface pressure for whole image (because of buildings and small deviations)
        
        """
        #file="../data/ASTER/ASTGTMV003_N51E019_dem.nc"
        fp = nc.Dataset(self.location.filename_aster)
        
       
        ilat = np.argmin(abs(fp['lat'][:] - self.topo['lat']))
        ilon = np.argmin(abs(fp['lon'][:] - self.topo['lon_360']))
        
        #print(fp['Band1'][ilat,ilon])
    
        #self.topo['height_2']=np.mean(fp['Band1'][ilat-20:ilat+20,ilon-20:ilon+20]) #calculate mean height over an area
        self.topo['height_2'] = fp['Band1'][ilat,ilon] #calculate height at this lat and height
        #print(height)
        print("height",self.topo['height_2'])
        
    def def_height(self,height):
        """
        if no ASTER file is available, then one can define a height 
        """
        self.topo['height_2'] = height
    
    def calculate_surf_press(self):
        """
        Calculate surface pressure at given height 
        https://rechneronline.de/barometer/hoehe.php
        Luftdruck auf Zielhöhe = Luftdruck auf aktueller Höhe * (1-Temperaturgradient*Höhenunterschied/Temperatur auf aktueller Höhe in Kelvin)^(0,03416/Temperaturgradient)
        pressure at target altitude = pressure at actual height * (1- temperature gradient * height difference / temperature at actual height in K ) ^ (0,03416 / temperatur gradient )
        
        Parameters 
   	    ----------
        None
                   
        Returns
    	------- 
        topo[surf_press]: float 
        	surface pressure at given latitude and longitude 
        """
        temp_gradient = 0.0065 #K/m
        height_diff = self.topo['height_2']-self.topo['height_1']
        self.topo['surf_press'] = self.topo['press_1']* (1- temp_gradient*height_diff/self.topo['temp_1'])**(0.03416/temp_gradient)
        