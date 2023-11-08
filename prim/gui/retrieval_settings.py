# -*- coding: utf-8 -*-
"""
Created on Sat Jun 24 00:08:07 2023

@author: kaiar
"""

from PyQt5 import QtCore,QtWidgets

class WaveEdges(QtWidgets.QWidget):
    def __init__(self):
        super(WaveEdges,self).__init__() 
        self.wave_from = QtWidgets.QSpinBox()
        self.wave_to = QtWidgets.QSpinBox()
        
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(QtWidgets.QLabel("wave from"))
        layout.addWidget(self.wave_from)
        layout.addWidget(QtWidgets.QLabel("to"))
        layout.addWidget(self.wave_to)
        layout.addStretch()
        self.setLayout(layout)

class IsoID(QtWidgets.QWidget):
    def __init__(self):
        super(IsoID,self).__init__() 
        self.names = []
        self.ids = []
        
        for i in range(4):
            self.names.append(QtWidgets.QLineEdit())
            self.ids.append(QtWidgets.QLineEdit())
            
        layout_names = QtWidgets.QHBoxLayout()
        layout_names.addWidget(QtWidgets.QLabel("Names"),1)
        for name in self.names:
            layout_names.addWidget(name,2)
            
        layout_ids = QtWidgets.QHBoxLayout()
        layout_ids.addWidget(QtWidgets.QLabel("IDs"),1)
        for i in self.ids:
            layout_ids.addWidget(i,2)
            
        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(layout_names)
        layout.addLayout(layout_ids)
        
        self.setLayout(layout)
        
class X0(QtWidgets.QWidget):
    def __init__(self):
        super(X0,self).__init__()
        molec=['CH4','H2O']
        self.molecs = []
        self.molec_start = []
        self.albedo = []
        
        for i in range(2):
            self.molecs.append(QtWidgets.QLabel(molec[i]))
            self.molec_start.append(QtWidgets.QSpinBox())
            
        for i in range(4):  
            self.albedo.append(QtWidgets.QSpinBox())
            
        layout_names = QtWidgets.QHBoxLayout()
        for molec in self.molecs:
            layout_names.addWidget(molec)
        
        layout_molec = QtWidgets.QHBoxLayout()
        for molec in self.molec_start:
            layout_molec.addWidget(molec)
        
        layout_albedo_name = QtWidgets.QHBoxLayout()
        for i in range(4):
            layout_albedo_name.addWidget(QtWidgets.QLabel("Alb"+str(i)))
        
        layout_albedo = QtWidgets.QHBoxLayout()
        for alb in self.albedo:
            layout_albedo.addWidget(alb)
            
        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(layout_names)
        layout.addLayout(layout_molec)
        layout.addLayout(layout_albedo_name)
        layout.addLayout(layout_albedo)
        
        self.setLayout(layout)
        

class RetrievalSetter(QtWidgets.QWidget):
    def __init__(self):
        super(RetrievalSetter,self).__init__() 
        self.wave_edges = WaveEdges()
        self.iso_id = IsoID()
        self.x0 = X0()
        
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.wave_edges)
        layout.addWidget(self.iso_id)
        layout.addWidget(self.x0)
        
        self.setLayout(layout)