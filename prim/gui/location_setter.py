# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:36:37 2023

@author: kaiar
"""

from PyQt5 import QtCore,QtWidgets

# class FloatBox(QtWidgets.QDoubleSpinBox):
#     def __init__(self):
#         super(Number,self).__init__()     
#         self.

# class NumberBox(QtWidgets.QSpinBox):
#     def __init__(self):
#         super(Number,self).__init__()     
#         self.

class LocationSetter(QtWidgets.QWidget):
    def __init__(self):
        super(LocationSetter,self).__init__() 
        self.latitude = QtWidgets.QDoubleSpinBox()
        self.longitude = QtWidgets.QDoubleSpinBox()
        
        self.i_from = QtWidgets.QSpinBox()
        self.j_from = QtWidgets.QSpinBox()         
        
        self.i_range = QtWidgets.QSpinBox()
        self.j_range = QtWidgets.QSpinBox() 
        
        left = QtWidgets.QFormLayout()
        left.addRow("Latitude",self.latitude)
        left.addRow("Longitude",self.longitude)
        
        middle = QtWidgets.QFormLayout()
        middle.addRow("Across from",self.i_from)
        middle.addRow("Along from",self.j_from) 
        
        right = QtWidgets.QFormLayout()
        right.addRow("range",self.i_range) 
        right.addRow("range",self.j_range)         
        
        layout = QtWidgets.QHBoxLayout()
        layout.addLayout(left)
        layout.addLayout(middle)
        layout.addLayout(right)
        
        self.setLayout(layout)
        