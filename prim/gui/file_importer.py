# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:30:47 2023

@author: kaiar
"""
from PyQt5 import QtCore,QtWidgets

from prim.gui.widgets import DragAndDropLabelwithButton

class FileImporter(QtWidgets.QWidget):
    def __init__(self):
        super(FileImporter,self).__init__()
        self.measurement = DragAndDropLabelwithButton("measurement file")
        self.sunpath = DragAndDropLabelwithButton("sun reference file")
        self.atm = DragAndDropLabelwithButton("atmosphere file")
        self.aster = DragAndDropLabelwithButton("ASTER file")
        
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.measurement)
        layout.addWidget(self.sunpath)
        layout.addWidget(self.atm)
        layout.addWidget(self.aster)
        
        self.setLayout(layout)
    