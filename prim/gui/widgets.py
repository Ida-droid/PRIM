#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 22:11:38 2023

@author: kai
"""
import os

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (
    QVBoxLayout,
    QLabel,
    QFrame,
    QPushButton,
    QFileDialog
)


class DragAndDropLabelwithButton(QFrame):
    fileDropped = pyqtSignal(str)

    def __init__(self, text, loadtext=None):
        super(DragAndDropLabelwithButton, self).__init__()
        self.setObjectName("DragAndDropLabelwithButton")
        self.file = text
        self.loadtext = loadtext
        self.setAcceptDrops(True)
        self.setAutoFillBackground(True)
        self.file_path = None
        self.setStyleSheet(
            "QFrame#DragAndDropLabelwithButton {border: 2px dashed #aaa; border-radius: 10px}"
        )
        self.label = QLabel("Drop " + self.file + " here\nor")
        self.label.setAlignment(Qt.AlignCenter)
        button = QPushButton("Select File")
        button.clicked.connect(self.chooseFile)

        layout = QVBoxLayout()
        layout.addStretch()
        layout.addWidget(self.label)
        layout.addWidget(button)
        layout.addStretch()
        self.setLayout(layout)

    def chooseFile(self):
        options = QFileDialog.Options()
        self.file_path, _ = QFileDialog.getOpenFileName(
            self,
            "QFileDialog.getSaveFileName()",
            "",
            "All Files (*);;Text Files (*.txt)",
            options=options,
        )
        if self.file_path:
            self.fileDropped.emit(self.file_path)
            self.setStyleSheet(
                "QFrame#DragAndDropLabelwithButton {border: 2px solid green; border-radius: 10px}"
            )
            filename = os.path.basename(self.file_path)
            if self.loadtext:
                self.label.setText(self.loadtext)
            else:
                self.label.setText(
                    filename
                    + " loaded\nDrag and drop anothor "
                    + self.file
                    + "\nto change the file\nor"
                )

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            self.setStyleSheet(
                "QFrame#DragAndDropLabelwithButton {border: 2px dashed #00f; border-radius: 10px}"
            )
            event.accept()
        else:
            event.ignore()

    def dragLeaveEvent(self, event):
        if not self.file_path:
            self.setStyleSheet(
                "QFrame#DragAndDropLabelwithButton {border: 2px dashed #aaa; border-radius: 10px}"
            )
        else:
            self.setStyleSheet(
                "QFrame#DragAndDropLabelwithButton {border: 2px solid green; border-radius: 10px}"
            )

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            self.file_path = str(url.toLocalFile())
            self.fileDropped.emit(self.file_path)
        self.setStyleSheet(
            "QFrame#DragAndDropLabelwithButton {border: 2px solid green; border-radius: 10px}"
        )
        filename = os.path.basename(self.file_path)
        if self.loadtext:
            self.label.setText(self.loadtext)
        else:
            self.label.setText(
                filename
                + " loaded\nDrag and drop anothor "
                + self.file
                + "\nto change the file\nor"
            )
        event.accept()
