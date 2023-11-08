# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 20:44:31 2023

@author: ida
"""

from prim.sun import Sun
from prim.prisma import PRISMA
import numpy as np




def initial_guess_albedo_0(settings, location, i, j):
    """
    Calculate an initial guess for the zeroth order albedo (c_0) based on PRISMA data.

    This function computes an initial guess for the albedo (c_0) by utilizing PRISMA
    measurements and sun-related data.

    Args:
        settings (object): An object containing various settings and configuration data.
        location (object): An object representing the location information.
        i (int): An index related to across-track axis in the PRISMA data.
        j (int): Another index related to along-track axis in the PRISMA data.

    Returns:
        float: An initial guess for the albedo (c_0) based on the provided data.

    Example:
        # Define the 'settings' and 'location' objects.
        initial_albedo_guess = initial_guess_albedo_0(settings_object, location_object, 1, 2)
        # This will calculate an initial albedo guess based on the specified settings 
        # and location and for the across-track 1 and along-track 2.
    """
    # Create a PRISMA measurement object.
    measurement = PRISMA(settings, location)
    # Read PRISMA data at specific indices.
    meas = measurement.read_netcdf4(i, j)
    # Calculate sun measurements using the measurement wavelength data.
    sun_meas = Sun(settings, meas["wave"]).meas
    # Extract values from PRISMA measurements.
    mu0 = meas["mu0"]
    y = meas["spec"]
    # Calculate the initial guess for albedo (A).
    c_0 = (mu0 * np.max(sun_meas) / (np.pi * np.max(y))) ** (-1)
    return c_0
    
    