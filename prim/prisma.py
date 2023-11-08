# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 18:44:58 2023

@author: ida
"""
import os 
import netCDF4 as nc
import numpy as np
import h5py

from prim.helper import StopExecution

def read_dataset(H5_object, path):
    data = H5_object[path]
    return np.array(data)

class PRISMA:
    def __init__(self,settings,location):
        self.wave_edges = settings.wave_edges
        self.lat = location.lat
        self.lon = location.lon
        self.name = location.name
        
    def create_netcdf4(filepath_l1, filepath_l2, name):
        """
        Read spectrum PRISMA
        
        Parameters
        ----------
        filenamel1 : str 
        	filepath to level 1 h5py file
        
        filenamel2 : str 
        	filepath to level 2 h5py file 
        		            
        name : str 
        	name of the wished netcdf4 file which is then stored at ../data/
        
        
        
        Returns
        -------
        a netcdf file stored at ../data
                   
        
        
        """
    
        if not (os.path.exists(filepath_l1) or os.path.exists(filepath_l2)):
            print("ERROR! read_spectrum_S5P: filename does not exist.")
            raise StopExecution
    
        # Read Level 1 data
        with h5py.File(filepath_l1, mode="r") as f1:
            data = read_dataset(f1, "/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/SWIR_Cube")
            # As described in the PRISMA product specification (p.107) 
            # Radiance = DN/ ScaleFactor - Offset
            # default value of ScaleFactor=100 "to be used in order to transform 
            # Swir Cube from Digital Number 16 bit to radiance W/m2srum" (p.107)
            # default value of Offset=0
            PRISMA_spec = (
                data[:, :, :] / f1.attrs["ScaleFactor_Swir"] - f1.attrs["Offset_Swir"]
            )  # ignore first two values (always zero)
    
            # list of wavelenght used for SWIR spectra is given in descending order. 
            # Flip the list to have them in ascending order.
            # For the same reason flip fwhm
            spec = np.copy(PRISMA_spec[:, :, :]) # shape sample, band, line [1000*][173][nHypAlongPixel]
            for i in range(0, 999):
                for j in range(0, 999):
                    spec[i, :, j] = np.flip(PRISMA_spec[i, :, j])
    
            wave_spec = np.flip(
                f1.attrs["List_Cw_Swir"]
            )  # ignore first two values (always zero)
    
            fwhm_SWIR = np.flip(f1.attrs["List_Fwhm_Swir"])
            cloud_SWIR = read_dataset(
                f1, "/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/Cloud_Mask"
            )[:, :]
            land_SWIR = read_dataset(
                f1, "/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/LandCover_Mask"
            )[:, :]
            #time_SWIR=Integration_Time#read_dataset(f1,"/HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Time")[:]
    
        # Read Level 2B data
        with h5py.File(filepath_l2, mode="r") as f2:
    
            latitude = f2["/HDFEOS/SWATHS/PRS_L2B_HCO/Geolocation Fields/Latitude"][:]
            latitude = latitude[:, :]
    
            longitude = f2["/HDFEOS/SWATHS/PRS_L2B_HCO/Geolocation Fields/Longitude"][:]
            longitude = longitude[:, :]
            solar_zenith_angle = read_dataset(
                f2, "/HDFEOS/SWATHS/PRS_L2B_HCO/Geometric Fields/Solar_Zenith_Angle"
            )[:, :]  # Dimension lines samples
            viewing_zenith_angle = read_dataset(
                f2, "/HDFEOS/SWATHS/PRS_L2B_HCO/Geometric Fields/Observing_Angle"
            )[:, :]  # Dimension lines samples
         
        os.chdir("/home/ijandl/PRISMA/") # TODO: the path should be adapted to your path
        # Save L1 and L2B data in one netcdf4 file called "name"
        fn = "data/" + name + ".nc"
        ds = nc.Dataset(fn, "w", format="NETCDF4")
    
        xpixel = ds.createDimension("lines", 1000)
        bands = ds.createDimension("bands", 173)
        ypixel = ds.createDimension("samples", 1000) #i,nAcrossPixel  SAMPLE
        spectra = ds.createVariable("spectra", "f4", ("lines", "bands", "samples"))
        spectra.units = "W/(str*m^2*um)"
        spectra[:, :, :] = spec
    
        wavelenght = ds.createVariable("wavelength", "f4", ("bands"))
        wavelenght.units = "nm"
        wavelenght[:] = wave_spec
    
        fwhm = ds.createVariable("fwhm", "f4", ("bands"))
        fwhm.units = "nm"
        fwhm[:] = fwhm_SWIR
    
        lat = ds.createVariable("lat", "f4", ("lines", "samples"))
        lat.units = "None"
        lat[:, :] = latitude
    
        sza = ds.createVariable("sza", "f4", ("lines", "samples"))
        sza.units = "°"
        sza[:, :] = solar_zenith_angle
    
        vza = ds.createVariable("vza", "f4", ("lines", "samples"))
        vza.units = "°"
        vza[:, :] = viewing_zenith_angle
    
        lon = ds.createVariable("lon", "f4", ("lines", "samples"))
        lon.units = "None"
        lon[:, :] = longitude
    
        cloud = ds.createVariable("cloud", "f4", ("lines", "samples"))
        cloud.units = "None"
        cloud[:, :] = cloud_SWIR
    
        land = ds.createVariable("land", "f4", ("lines", "samples"))
        land.units = "None"
        land[:, :] = land_SWIR
        
        ##time = ds.createVariable("time", "f4", ("samples"))
        #time.units = "MJD2000Decimaldays"
        #time[:] = time_SWIR
    
        ds.close()
        return

    def read_netcdf4(self, i, j,hide=True):
    
        """
        function that opens the stored .nc file and selects the spectra 
        for i (nAcrossPixel) and j (nAlongPixel) 
        in a desired wavelenght range [wave_start,wave_end]
        
        Parameters 
        ----------
        name: str 
        	name of .nc data file stored under ../data/
        wave_start: int 
        	desired start of wavelenght grid
        	
        wave_end: int 
        	desired end of wavlenght grid
        	
        i: int 
        	desired nAcrossPixel  SAMPLE
        j: int 
        	desired nAlongPixel LINE
        	
        lat_station: float 
        	latitude where to find the plume
        	
        lon_station: float
        	longitude where to find the plume
        	
        hide: boolean, 
        	if it is set False, than all the pixels are printed where to find the plume. Default value: False
        
        Returns
        ------- 
        meas: dictionary, 
        	captures all relevant Level 1 and Level 2 data from PRISMA at given Pixel 
        meas["wave"]: lst
        	wavelenght range of measured spectra
        	
        meas["spec"]: lst
        	spectra of PRISMA in the wavelenght range
        	
        meas["fwhm"]: float
        	fwhm of the Spectra 
        
        meas["sza"]: float
        	solar zenith angle
        	
        meas["vza"]: float 
        	viewing zenith angle 
        	
        meas["cloud"]: int 
        	information of cloudy pixel (0 if there is no cloud, 1 for cloud, 255 for error) 
        	
        meas["irr_2342nm"]: 2-d array 
        	spectra of infrared image in 2340nm range for quick look at the image to observe CH4 absorption
        	
        meas["irr_2061nm"]: 2-d array 
        	spectra of infrared image in 2000nm range for quick look at the image to observe CO2 absorption
        
        meas["lat"]: float 
        	latutiude at [i,j]
        	
        meas["lon"]: float 
        	longitude at [i,j]
        	
        	
        	
        """
        wave_start, wave_end = self.wave_edges
        lat_station = self.lat
        lon_station = self.lon
        fp = "../../data/locations/" + self.name + ".nc"
        
        n = nc.Dataset(fp)
        spectra = n["spectra"]
        wave_spec = n["wavelength"]
        wave_spec = np.asarray(wave_spec)
        start = np.argmin(abs(wave_spec - wave_start))
        end = np.argmin(abs(wave_spec - wave_end))
        
        meas={}
        meas["spec"] = spectra[i, start:end+1, j]
        meas["spec_per_column"] = spectra[i, start:end+1, 0:1000]
        meas["wave"] = wave_spec[start:end+1]
        
        meas["irr_2342nm"] = spectra[:, 150, :]
        meas["irr_2061nm"] =spectra[:, 114, :]
        meas["fwhm"] = n["fwhm"][start:end+1]
        meas["sza"] = n["sza"][i, j]
        meas["mu0"] = np.cos(np.deg2rad(meas["sza"]))   
        meas["vza"] = n["vza"][i, j]
        meas["muv"] = np.cos(np.deg2rad(meas["vza"]))
        meas["cloud"] = n["cloud"][i, j]
        longitude=n['lon'][:,:]
        latitude=n['lat'][:,:]
        #meas["time"]=n['time'][:]
        meas["lon"]=n['lon'][i,j] #longitude at that point
        meas["lat"]=n['lat'][i,j] #latitude at that point
        if hide==False:
                for i in range(0,999):
                    for j in range(0,999):
                        if (abs(longitude[i][j]-lon_station)<0.0008 and abs(latitude[i][j]-lat_station)<0.0008):
                            print(i,j)
                            m=i
                            n=j
                            print('latitude_l2',latitude[m][n])
                            print('longitude_l2',longitude[m][n])
        
        return meas
