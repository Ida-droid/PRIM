# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 19:54:13 2023

@author: idajandl@gmx.de
"""
import os

import numpy as np
from numpy import exp
import netCDF4 as nc

from prim.helper import StopExecution
from prim.topography import Topography

# Some global constants
hplanck = 6.62607015E-34  # Planck's constant [J/s]
kboltzmann = 1.380649E-23 # Boltzmann's constant [J/K]
clight = 2.99792458E8     # Speed of light [m/s]
e0   = 1.6021766E-19      # Electron charge [C]
g0   = 9.80665            # Standard gravity at ground [m/s2]
NA   = 6.02214076E23      # Avogadro number [#/mol]
Rgas = 8.314462           # Universal gas constant [J/mol/K]
XO2  = 0.2095             # Molecular oxygen mole fraction [-]
XCH4 = 1.8E-6             # std CH4 volume mixing ratio 1800 ppb
XCO2 = 410E-6
MDRYAIR = 28.9647E-3      # Molar mass dry air [kg/mol]
MH2O = 18.01528E-3        # Molar mass water [kg/mol]
MCO2 = 44.01E-3           # Molar mass CO2 [kg/mol]
MCH4 = 16.04E-3           # Molar mass CH4 [kg/mol]
LH2O = 2.5E6              # Latent heat of condensation for water [J/kg]
PSTD = 1013.25            # Standard pressure [hPa]
   
    
class Atmosphere:
    """
    The atmosphere_data class collects data of
    the thermodynamic state and the composition of
    the atmosphere.
    
    CONTAINS
    
    method __init__(self,zlay,zlev) |
    method get_data_AFGL(self,filename) |
    method get_data_ECMWF_ads_egg4(self,filename, month, longitude, latitude)
    """
    ###########################################################
    def __init__(self,settings,location):
        """
        init class
        
        Parameters 
    	----------
        zlay: array of vertical height layers, mid-points [nlay] [m]
        zlev: array of vertical height levels, boundaries [nlev=nlay+1] [m]
        psurf: scalar of surface pressure [hPa]
        
        Returns
    	-------
        atmo: dictionary with atmospheric data
        atmo[zlev]: array of vertical height levels, boundaries [nlev=nlay+1] [m]
        atmo[zlay]: array of vertical height layers [nlay] [m]
        """
        self.atmo   = settings.atm_dict
        topography = Topography(location)
        self.lon= location.lon % 360  
        self.lat = location.lat
        self.atmo['psurf']  = topography.topo['surf_press']
        file= location.filename_atm
        if file[0].lower()=="afgl":
            self.get_data_AFGL(file[1])
        elif file[0]=="cams":
            self.get_data_cams(file[1])
    
    def clausius_clapeyron(self,T):
        """
        Calculates saturation water vapor pressure according to
        Clausius-Clapeyron equation, exponential approximation
        
        Parameters
        ----------
        Temperature : float [K]
        
        Returns
        -------
        saturation water vapor pressure [hPa]
        """
        # check whether input is in range
        while True:
            if T > 0.:
               break
            else:
               print("ERROR! clausius_clapeyron: input out of range.")
               raise StopExecution
                
        T0 = 273.16 # ref. temperature[K], triple point
        p0 = 611.66 # ref. pressure [Pa], triple point
        
        # Saturation water vapor pressure according to
        # Clausius-Clapeyron equation, exponential approximation
        pS = p0 * exp(-LH2O*MH2O/Rgas*(1./T - 1./T0))
        
        return pS/100. # [hPa]
    
    ###########################################################
    def get_data_AFGL(self,filename):
        """
        Read atmospheric data from AFGL database
        file, interpolate on output height grid
        
        Parameters 
    	----------
        filename: str 
        	file with AFGL data
        	
        Returns
    	-------  
        atmo[tlev]: temperature level profile [nlev] [K]
        atmo[tlay]: temperature layer profile [nlev] [K]
        atmo[plev:  pressure level profile [nlev] [hPa]
        atmo[play:  pressure layer profile [nlay] [hPa]
        atmo[AIR]:  air partial column profile [nlay] [#/m^2]
        atmo[O3]:   o3 partial column profile [nlay] [#/m^2]
        atmo[H2O]:  h2o prtial column profile [nlay] [#/m^2]
        atmo[CO2]:  co2 partial column profile [nlay] [#/m^2]
        atmo[NO2]:  no2 partial column profile [nlay] [#/m^2]
        atmo[O2]:   o2 partial column profile [nlay] [#/m^2]
        atmo[CH4]:  CH4 partial column profile [nlay] [#/m^2]
        """
        
        #some constants
        Rgas = 8.3144598    #universal gas constant [J/(mol*K)]
        grav = 9.80665      #gravitational acceleration [m/s2]
        Mair = 0.0289644    #molar mass of Earth's air [kg/mol]
        Avog = 6.022E23     #Avogadro constant [part./mol]
        #check whether input is in range
        while True:
            if os.path.exists(filename):
               break
            else:
               print("ERROR! atmosphere_data.get_data_AFGL: file does not exist.")
               raise StopExecution
        
        # Read AFGL file
        atm_in = np.genfromtxt(filename,skip_header=2) 

        zalt_in  = atm_in[:,0]*1.E3         # height [km] -> [m]
    
        press_in = atm_in[:,1]              # pressure [hPa]
        temp_in  = atm_in[:,2]              # temperature [K]
        air_in   = atm_in[:,3]              # air number density [#/cm^3] 
        o3_in    = atm_in[:,4]/air_in       # o3 number density -> mole fraction [-]
        o2_in    = atm_in[:,5]/air_in       # o2 number density -> mole fraction [-]
        h2o_in   = atm_in[:,6]/air_in       # h2o number density -> mole fraction [-]
        co2_in   = atm_in[:,7]/air_in       # co2 number density -> mole fraction [-]
        no2_in   = atm_in[:,8]/air_in       # no2 number density -> mole fraction [-]
        nlev_in  = zalt_in.size             # number of input levels

        sp = NA/(MDRYAIR*g0)*1.E-2 # [#/m^2 * 1/hPa] air column above P is P*NA/Mair/g from p = m*g/area

        #truncate or extrapolate the AFGL profile depending on psurf
        #print(press_in[press_in.size-1], self.atmo['psurf'])
        if press_in[press_in.size-1] < self.atmo['psurf']:
            #extrapolation required
            dz       = np.log(self.atmo['psurf']/press_in[press_in.size-1])*Rgas*temp_in[temp_in.size-1]/(grav*Mair)
            press_in = np.append(press_in, self.atmo['psurf'])
            zalt_in  = np.append(zalt_in-dz, 0.)
            temp_in  = np.append(temp_in, temp_in[temp_in.size-1]) 
            air_in   = np.append(air_in,press_in[press_in.size-1]*Avog*1.E-4/(Rgas *temp_in[temp_in.size-1]))
            o3_in    = np.append(o3_in, o3_in[o3_in.size-1]) 
            o2_in    = np.append(o2_in, o2_in[o2_in.size-1]) 
            h2o_in   = np.append(h2o_in, h2o_in[h2o_in.size-1]) 
            co2_in   = np.append(co2_in, co2_in[co2_in.size-1]) 
            no2_in   = np.append(no2_in, no2_in[no2_in.size-1]) 
            nlev_in  = nlev_in-1
        elif press_in[press_in.size-1] > self.atmo['psurf']:
        
            #print('yes')
            #interpolation required
            intv = np.searchsorted(press_in,self.atmo['psurf'])  # self.atmo['psurf'] is in the interval [press_in[intv], press_in[intv-1]]
            press_in = np.append(press_in[0:intv],self.atmo['psurf'])
            temp_in  = temp_in[0:intv+1]
            air_in   = np.append(air_in[0:intv],press_in[press_in.size-1]*Avog*1.E-4/(Rgas *temp_in[temp_in.size-1]))
            o3_in    = o3_in[0:intv+1]
            o2_in    = o2_in[0:intv+1]
            h2o_in   = h2o_in[0:intv+1]
            co2_in   = co2_in[0:intv+1]
            no2_in   = no2_in[0:intv+1]
            zalt_in  = zalt_in[0:intv+1]
            dz = np.log(press_in[press_in.size-1]/press_in[press_in.size-2])*Rgas*temp_in[temp_in.size-1]/(grav*Mair)
            zalt_in  = np.append(zalt_in[0:intv]-zalt_in[intv-1]+dz,0)
            
        # Interpolate temperature [K] on output layers
        # Flip arrays because our heights are descending
        # (from top to bottom), while np.interp expects ascending order
        
        self.atmo['tlev'] = np.flip(np.interp(np.flip(self.atmo['zlev']),np.flip(zalt_in),np.flip(temp_in)))
        self.atmo['tlay'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(temp_in)))
        
        # Calculate pressure [hPa] on output levels and layers
        self.atmo['plev'] = np.flip(np.interp(np.flip(self.atmo['zlev']),np.flip(zalt_in),np.flip(press_in)))
        self.atmo['play'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(press_in)))
        
        # Calculate the vertical column of air above pressure level 
        # and use this to calculate the partial vertical air columns per layer [#/m^2].
        # Partial columns have the advantage that multiplication with cross sections
        # yields optical depth.
        
        nlev = len(self.atmo['zlev'])
        sp = NA/(MDRYAIR*g0)*1.E-2 # [#/m^2 * 1/hPa] air column above P is P*NA/Mair/g from p = m*g/area
        vc_air = sp*self.atmo['plev'] # air column [#/m^2] above pressure level
        print(len(self.atmo['plev']),nlev)
        self.atmo['AIR'] = (vc_air[1:nlev]-vc_air[0:nlev-1])/ (1.+h2o_in[0:nlev-1]/1.60855)#[#/m^2]
        #self.atmo['AIR'][0] = vc_air[0]#[#/m^2] uppermost layer extends to infinity in terms of number of molecules
        self.atmo['AIR'][0] =self.atmo['AIR'][0] + vc_air[0]/ (1.+h2o_in[0]/1.60855)
        
        # Interpolate mole fractions on output height grid
        # and then calculate partial columns per layer [#/m^2]
        # ozone
        self.atmo['O3']  = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(o3_in)))*self.atmo['AIR']
        # water vapor
        self.atmo['H2O'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(h2o_in)))*self.atmo['AIR']
        # co2
        self.atmo['CO2'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(co2_in)))*self.atmo['AIR']
        # no2
        self.atmo['NO2'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(no2_in)))*self.atmo['AIR']
        # ch4 use a constant mixing ratio
        self.atmo['CH4'] = np.repeat(XCH4,20)*self.atmo['AIR']
        

    ###########################################################
    def get_data_cams(self,filename):
        """
        Read atmospheric data provided by ECMWF ADS EGG4 run:
        https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-ghg-reanalysis-egg4-monthly
         
        Download script:
         
        import cdsapi
        c = cdsapi.Client()
        c.retrieve(
             'cams-global-ghg-reanalysis-egg4-monthly',
             {
                 'variable': [
                     'carbon_dioxide', 'geopotential', 'methane',
                     'relative_humidity', 'temperature',
                 ],
                 'pressure_level': [
                     '1', '2', '3',
                     '5', '7', '10',
                     '20', '30', '50',
                     '70', '100', '150',
                     '200', '250', '300',
                     '400', '500', '600',
                     '700', '800', '850',
                     '900', '925', '950',
                     '1000',
                 ],
                 'year': '2016',
                 'month': [
                     '01', '02', '03',
                     '04', '05', '06',
                     '07', '08', '09',
                     '10', '11', '12',
                 ],
                 'product_type': 'monthly_mean',
                 'format': 'netcdf',
             },
             'download.nc'
        
        Parameters 
    	----------
        filename: str 
        	filepath to ECMWF netcdf file
        month: int  
        	[01 ... 12]
        longitude: int
        	[0 ... 360 degree]
        latitude: int 
        	 [-90 ... 90 degree]
        Returns
    	------- 
        temp: temperature profile [nlay] [K]
        plev: pressure level profile [nlev] [hPa]
        atmo[tlev]: temperature level profile [nlev] [K]
        atmo[tlay]: temperature layer profile [nlev] [K]
        atmo[plev]: pressure level profile [nlev] [hPa]
        atmo[play]: pressure layer profile [nlay] [hPa]
        atmo[AIR]: air partial column profile [nlay] [#/m^2]
        atmo[H2O]: h2o partial column profile [nlay] [#/m^2]
        atmo[CO2]: co2 partial column profile [nlay] [#/m^2]
        atmo[CH4]: no2 partial column profile [nlay] [#/m^2]
        atmo[O2]: o2 partial column profile [nlay] [#/m^2]
        """    
        # check whether input is in range
        while True:
            if os.path.exists(filename)  and -90. <= self.lat <= 90. and 0. <= self.lon <= 360.:
               break
            elif not os.path.exists(filename):
               print("ERROR! read_ecmwf_ads_egg4: filename does not exist.")
               raise StopExecution
            else:
               print("ERROR! read_ecmwf_ads_egg4: input out of range.")
               raise StopExecution
        Rgas = 8.3144598    #universal gas constant [J/(mol*K)]
        grav = 9.80665      #gravitational acceleration [m/s2]
        Mair = 0.0289644    #molar mass of Earth's air [kg/mol]
        Avog = 6.022E23     #Avogadro constant [part./mol]    
        # Open netcdf file
        ds = nc.Dataset(filename)
        #print(ds.variables)

        # Select month index, latitude/longitude index (next neighbour)
        itime = 0
        ilat = np.argmin(abs(ds['latitude'][:] - self.lat))
        ilon = np.argmin(abs(ds['longitude'][:] - self.lon))
        
        # ECMWF: Geopotential [m2 s-2] converted to height [m], approximate use of g0 
        zalt_in = np.array([d/g0 for d in ds['z'][itime,:,ilat,ilon]])
        nlev_in = zalt_in.size
        # ECMWF: Pressure [hPa]
        press_in = ds['level'][:]
        # ECMWF: Temperature [K]
        temp_in = ds['t'][itime,:,ilat,ilon]
        mean_t=np.mean(temp_in)

        # ECMWF: Humidity [%] converted to water vapor mole fraction [mol/mol] via Clausius-Clapeyron equation
        pS = [self.clausius_clapeyron(Ti) for Ti in temp_in]
        # ECMWF: Mole fraction is partial pressure over dry total pressure, partial pressure is rel. hum. * sat. vapor pressure
        h2o_in = np.array([d/100.*pSi/(pi - d/100.*pSi) for d,pSi,pi in zip(ds['r'][itime,:,ilat,ilon],pS,press_in)])
        # ECMWF: Carbon dioxide mass mixing ratio [kg kg-1] converted to mole fraction [mol/mol]
        co2_in = np.array([d/MCO2*MDRYAIR for d in ds['co2'][itime,:,ilat,ilon]])
        # ECMWF: Methane mass mixing ratio [kg kg-1] converted to mole fraction [mol/mol]
        ch4_in = np.array([d/MCH4*MDRYAIR for d in ds['ch4'][itime,:,ilat,ilon]])
        #print("press_in",press_in,"temp_in",temp_in, "h20_in",h2o_in, "co2_in",co2_in, "ch4_in",ch4_in,"zalt_in",zalt_in)
        ds.close
        
        
        ########################################################################################################################################
        
        #air_in=sp
        #truncate or extrapolate the AFGL profile depending on psurf
        
        if press_in[press_in.size-1] < self.atmo['psurf']:
            #extrapolation required
            dz       = np.log(self.atmo['psurf']/press_in[press_in.size-1])*Rgas*temp_in[temp_in.size-1]/(grav*Mair)
            press_in = np.append(press_in, self.atmo['psurf'])
            zalt_in  = np.append(zalt_in+dz, 0.) #TODO
            temp_in  = np.append(temp_in, temp_in[temp_in.size-1]) 
            #air_in   = np.append(air_in,press_in[press_in.size-1]*Avog*1.E-4/(Rgas *temp_in[temp_in.size-1]))
            #o3_in    = np.append(o3_in, o3_in[o3_in.size-1]) 
            #o2_in    = np.append(o2_in, o2_in[o2_in.size-1]) 
            h2o_in   = np.append(h2o_in, h2o_in[h2o_in.size-1]) 
            co2_in   = np.append(co2_in, co2_in[co2_in.size-1]) 
            ch4_in   = np.append(ch4_in, ch4_in[ch4_in.size-1])
            #no2_in   = np.append(no2_in, no2_in[no2_in.size-1]) 
            nlev_in  = nlev_in-1
        elif press_in[press_in.size-1] > self.atmo['psurf']:
            #interpolation required
            print("yes")
            intv = np.searchsorted(press_in,self.atmo['psurf'])  # self.atmo['psurf'] is in the interval [press_in[intv], press_in[intv-1]]
            press_in = np.append(press_in[0:intv],self.atmo['psurf'])
            temp_in  = temp_in[0:intv+1]
            #air_in   = np.append(air_in[0:intv],press_in[press_in.size-1]*Avog*1.E-4/(Rgas *temp_in[temp_in.size-1]))
            #o3_in    = o3_in[0:intv+1]
            #o2_in    = o2_in[0:intv+1]
            h2o_in   = h2o_in[0:intv+1]
            co2_in   = co2_in[0:intv+1]
            ch4_in   = ch4_in[0:intv+1]
            #no2_in   = no2_in[0:intv+1]
            zalt_in  = zalt_in[0:intv+1]
            dz = np.log(press_in[press_in.size-1]/press_in[press_in.size-2])*Rgas*temp_in[temp_in.size-1]/(grav*Mair)
            zalt_in  = np.append(zalt_in[0:intv]-zalt_in[intv-1]+dz,0)
        
        # Interpolate temperature [K] on output layers
        # Flip arrays because our heights are descending
        # (from top to bottom), while np.interp expects ascending order
        self.atmo['tlay'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(temp_in)))
        self.atmo['tlev'] = np.flip(np.interp(np.flip(self.atmo['zlev']),np.flip(zalt_in),np.flip(temp_in)))
        
        # Calculate pressure [hPa]
        self.atmo['play'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(press_in)))
        self.atmo['plev'] = np.flip(np.interp(np.flip(self.atmo['zlev']),np.flip(zalt_in),np.flip(press_in)))
        
        # Calculate the vertical column of air above pressure level 
        # and use this to calculate the partial vertical air columns per layer [#/m^2].
        # Partial columns have the advantage that multiplication with cross sections
        # yields optical depth.
        
        h2o_in_corr=np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(h2o_in)))
        nlev = len(self.atmo['zlev'])
        sp = NA/(MDRYAIR*g0)*1.E+2 # [#/m^2 * 1/hPa] air column above P is P*NA/Mair/g from p = m*g/area
        vc_air = sp*self.atmo['plev'] # air column [#/m^2] above pressure level
        
        self.atmo['AIR'] = (vc_air[1:nlev]-vc_air[0:nlev-1])/ (1.+h2o_in_corr[0:nlev-1]/1.60855)#[#/m^2]
        
        #self.atmo['AIR'][0] = vc_air[0]#[#/m^2] uppermost layer extends to infinity in terms of number of molecules
        self.atmo['AIR'][0] =self.atmo['AIR'][0] + vc_air[0]/ (1.+h2o_in_corr[0]/1.60855)
              
        # Interpolate mole fractions on output height grid
        # and then calculate partial columns per layer [#/m^2]
        # water vapor
        self.atmo['H2O'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(h2o_in)))*self.atmo['AIR']
        # co2
        self.atmo['CO2'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(co2_in)))*self.atmo['AIR']
        # no2
        self.atmo['CH4'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(ch4_in)))*self.atmo['AIR']
        #n2o
        self.atmo['N2O'] = np.flip(np.interp(np.flip(self.atmo['zlay']),np.flip(zalt_in),np.flip(co2_in)))*self.atmo['AIR']/414.5*0.33
        
        #self.atmo['CH4'] = XCH4*self.atmo['AIR']
        # o2 use a constant mixing ratio
        self.atmo['O2']  = XO2*self.atmo['AIR']
        