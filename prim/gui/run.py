# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:20:29 2023

@author: kaiar
"""
import os
import sys
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import QSize
from PyQt5.QtWidgets import QApplication

from prim.gui.main import MainWindow

def run():
    myappid = "bh_gui"  # arbitrary string
    try:
        from ctypes import windll
        windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
        print('Windows detected')
    except:
        pass
    #set fontsize:
    # Get the default font
    # defaultFont = QtWidgets.QApplication.font()
    
    # # Define a scaling factor for the font size (adjust as needed)
    # scalingFactor = 1.2
    
    # # Calculate the scaled font size
    # scaledFontSize = defaultFont.pointSizeF() * scalingFactor
    
    # # Create a font with the scaled font size
    # scaledFont = defaultFont
    # scaledFont.setPointSizeF(scaledFontSize)
    
    # # Set the scaled font as the application font
    # QtWidgets.QApplication.setFont(scaledFont)

    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec_()