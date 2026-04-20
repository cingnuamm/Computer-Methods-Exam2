# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Problem1.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(400, 300)
        self.verticalLayout = QtWidgets.QVBoxLayout(Form)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gb_Input = QtWidgets.QGroupBox(Form)
        self.gb_Input.setObjectName("gb_Input")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.gb_Input)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.layout_GridInput = QtWidgets.QGridLayout()
        self.layout_GridInput.setObjectName("layout_GridInput")
        self.verticalLayout_2.addLayout(self.layout_GridInput)
        self.verticalLayout.addWidget(self.gb_Input)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.gb_Input.setTitle(_translate("Form", "Input"))

