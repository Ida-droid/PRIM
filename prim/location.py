#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 17:36:48 2023

@author: ida
"""

class Location:
    def __init__(self, name, lat, lon, filepath_aster, filepath_atm, lst_index):
        """
        Represents a geographic location (e.g. point source location) with associated data.

        Args:
            name (str): The name or identifier for the location.
            lat (float): The latitude of the location.
            lon (float): The longitude of the location.
            filepath_aster (str): The file path to ASTER data associated with the location.
            filepath_atm (tuple): A tuple containing two elements. The first element is a string
                                 specifying the data source (e.g., 'cams'), and the second element
                                 is the file path to atmospheric data associated with the location.
            lst_index (list): A list containing four integer values that represent the PRISMA datacubes indecees
                             associated with the location [y_0, length(y), x_0, length(x)]

        Attributes:
            name (str): The name or identifier for the location.
            lat (float): The latitude of the location.
            lon (float): The longitude of the location.
            filepath_aster (str): The file path to ASTER data associated with the location.
            filepath_atm (tuple): A tuple containing two elements. The first element is a string
                                 specifying the data source, and the second element is the file path
                                 to atmospheric data associated with the location.
            lst_index (list): A list containing four integer values representing PRISMA datacubes indecees, where the plume is assumed.

        Example:
            location = Location("Turkmenistan_2020_07_21",
                                38.558, 54.2019,
                                "../../data/ASTER/ASTGTMV003_N38E054_dem.nc",
                                ('cams', "../../data/atmosphere/downloadTurkmenistanpressureleveldaily.nc"),
                                [585, 70, 280, 45]
                                )
        """
        self.name = name
        self.lat = lat
        self.lon = lon
        self.filepath_aster = filepath_aster
        self.filepath_atm = filepath_atm
        self.lst_index = lst_index

        
        
class ExampleLocations:
    def __init__(self):
        self.Turkmenistan_2020_07_21 = Location("Turkmenistan_2020_07_21",
                                                38.558,54.2019,
                                                "../../data/ASTER/ASTGTMV003_N38E054_dem.nc",
                                                ('cams',"../../data/atmosphere/downloadTurkmenistanpressureleveldaily.nc"),
                                                [585,70,280,45]
                                                )
        self.Turkmenistan_2020_10_10 = Location("Turkmenistan_2020_10_10",
                                                38.49393,54.19764,
                                                "../../data/ASTER/ASTGTMV003_N38E054_dem.nc",
                                                ('cams',"../../data/atmosphere/download_Turkmenistan_2020_10_10.nc"),
                                                [540,70,550,50]
                                                )
        self.Turkmenistan_2020_03_27 = Location("Turkmenistan_2020_03_27",
                                                39.4968,53.63771,
                                                "../../data/ASTER/ASTGTMV003_N39E053_dem.nc",
                                                ('cams',"../../data/atmosphere/download_Turkmenistan_2020_03_27.nc"),
                                                [430,70,550,100]
                                                )
        self.China_Guanter = Location("China_Guanter",
                                                36.233333,112.946389,
                                                "../../data/ASTER/ASTGTMV003_N36E112_dem.nc",
                                                ('cams',"../../data/atmosphere/download_China_Guanter.nc"),
                                                [300,120,400,125]
                                                )
        self.InterMountainPowerplant = Location("InterMountainPowerplant",
                                                39.364358,-112.535342,
                                                "../../data/ASTER/ASTGTMV003_N39W113_dem.nc",
                                                ('cams',"../../data/atmosphere/download_Intermountain_pressurelevel_daily.nc"),
                                                [440,80,350,130]
                                                )
        self.Belchatow = Location("Belchatow",
                                                51.268526,19.331206,
                                                "../../data/ASTER/ASTGTMV003_N51E019_dem.nc",
                                                ('cams',"../../data/atmosphere/download_Belchatow.nc"),
                                                [450,80,400,130]
                                                )
        self.Four_Corners_2020_07_31 = Location("Four_Corners_2020_07_31",
                                                36.685618,-108.4808,
                                                "../../data/ASTER/ASTGTMV003_N36W109_dem.nc",
                                                ('cams',"../../data/atmosphere/download_Four_Corners_2020_07_31.nc"),
                                                [420,200,580,270]
                                                )
        self.San_Juan_Power_Plant = Location("San_Juan_Power_Plant",
                                                36.803985,-108.440156,
                                                "../../data/ASTER/ASTGTMV003_N36W109_dem.nc",
                                                ('cams',"../../data/atmosphere/download_Four_Corners_2020_07_31.nc"),
                                                [470,80,220,130]
                                                )
        
        
        
