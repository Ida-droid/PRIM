# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 19:33:21 2023

@author: ida
"""
import numpy as np

def slit_conv(fwhm, wave, wave_meas): 
    """
    Perform slit convolution on a set of wavelength measurements.
    
    Slit convolution is a mathematical operation used in spectroscopy to simulate
    the effect of a finite-width slit on spectral data.
    
    Args:
        fwhm (array-like): An array of Full Width at Half Maximum (FWHM) values,
            representing the width of the slit for each wavelength measurement.
        wave (array-like): An array of reference wavelengths (e.g. line by line wavelength grid).
        wave_meas (array-like): An array of wavelengths for which the slit
            convolution will be performed.
    
    Returns:
        numpy.ndarray: A 2D numpy array representing the slit-convolved data.
            Each row corresponds to a wavelength measurement, and each column
            corresponds to a wavelength from the reference data.
    
    Note:
        - Slit convolution is performed for each wavelength measurement using
          the provided FWHM value for that specific wavelength.
        - The output array represents the effect of the slit on the spectral data.
    
    Example:
        fwhm_values = [2.0, 1.5, 2.5]
        reference_wavelengths = [500, 600, 700, 800]
        measured_wavelengths = [550, 650]
        slit_convolved_data = slit_conv(fwhm_values, reference_wavelengths, measured_wavelengths)
        sun_reference=slit_convolved_data.dot(sun_meas)
    """
   
    dmeas = wave_meas.size
    dlbl  = wave.size
    slit  = np.zeros(shape=(dmeas,dlbl))
    #const = fwhm**2/(4*np.log(2))
    for l,wmeas in enumerate(wave_meas):
        const = fwhm[l]**2/(4*np.log(2))
        wdiff = wave - wmeas
        slit[l,:] = np.exp(-wdiff**2/const)
        slit[l,:] = slit[l,:]/np.sum(slit[l,:])

    #rad_conv = slit.dot(rad_lbl)
    return slit

class StopExecution(Exception):
    """
    Class to stop execution in jupyter notebook and ipython
    Call via 'raise StopExecution'
    """    
    def _render_traceback_(self):
        pass

def save3d(array3d, txt):
    """
    Save a 3D NumPy array to a text file after reshaping it.
    
    This function reshapes the 3D array into a 2D array and then saves it to a
    text file with the given name. It is often used to save 3D data in a more
    human-readable format for later analysis.
    
    Args:
        array3d (numpy.ndarray): A 3D NumPy array to be saved.
        txt (str): A string representing the desired name of the text file.
            The file will be saved in the "../data/" directory.
    
    Returns:
        int: The third dimension size (number of layers) of the original 3D array.
    
    Example:
        data_3d = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
        num_layers = save3d(data_3d, "my_data")
        # This will save the reshaped data in a file called "my_data.txt" and return 2.
    """
    array_reshaped = array3d.reshape(array3d.shape[0], -1)
    np.savetxt("../data/" + txt+".txt", array_reshaped) 
    return  array3d.shape[2] 

def load3d(filepath, array3d_shape):
    """
    Load a 3D NumPy array from a text file and reshape it to the specified shape.

    This function loads a 3D array from a text file and reshapes it to the specified shape
    provided as 'array3d_shape'. It is often used to load previously saved 3D data from a
    text file and restore its original shape for further analysis.

    Args:
        filepath (str): The path to the text file containing the 3D data.
        array3d_shape (int): The third dimension size (number of layers) of the original 3D array.

    Returns:
        numpy.ndarray: A 3D NumPy array with the shape specified in 'array3d_shape',
            representing the loaded data.

    Example:
        # Assuming 'data.txt' contains previously saved 3D data, and the original shape
        # of the data was (2, 2, 2):
        restored_data = load3d("data.txt",  2)
        # This will load the data from the file and reshape it to the original shape.
    """
    loaded_arr = np.loadtxt(filepath) 
    res_shape = array3d_shape # entspricht residuum_array.shape[2] 
    load_original_arr = loaded_arr.reshape(loaded_arr.shape[0], loaded_arr.shape[1] // res_shape, res_shape) 
    return load_original_arr