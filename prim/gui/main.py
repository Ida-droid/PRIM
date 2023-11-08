# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:08:40 2023

@author: kaiar
"""

from PyQt5 import QtCore,QtWidgets

from prim.gui.file_importer import FileImporter
from prim.gui.location_setter import LocationSetter
from prim.gui.retrieval_settings import RetrievalSetter
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.files = FileImporter()
        self.location = LocationSetter()
        self.settings = RetrievalSetter()
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(QtWidgets.QLabel('Import files'))
        layout.addWidget(self.files)
        layout.addWidget(QtWidgets.QLabel('Location settings'))
        layout.addWidget(self.location)   
        layout.addWidget(QtWidgets.QLabel('Retrieval settings'))
        layout.addWidget(self.settings)  
        layout.addStretch()
        widget = QtWidgets.QWidget()
        widget.setStyleSheet(" font-size: 12pt; ")
        widget.setLayout(layout)
        self.setCentralWidget(widget)
        self.show()
        
        