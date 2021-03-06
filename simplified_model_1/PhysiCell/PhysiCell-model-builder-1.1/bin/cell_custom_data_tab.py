"""
Authors:
Randy Heiland (heiland@iu.edu)
Adam Morrow, Grant Waldrow, Drew Willis, Kim Crevecoeur
Dr. Paul Macklin (macklinp@iu.edu)

--- Versions ---
0.1 - initial version
"""

import sys
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QDoubleValidator

class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)

# Overloading the QLineEdit widget to let us map it to its variable name. Ugh.
class MyQLineEdit(QLineEdit):
    # prev_vname = None  # previous variable name (before changed)
    vname = None  # variable name
    idx = None  # index

class CellCustomData(QWidget):
    def __init__(self):
        super().__init__()

        # self.current_param = None
        self.xml_root = None
        self.celldef_tab = None
        self.count = 0
        self.max_rows = 99  # initially
        self.custom_vars_to_delete = []

        #-------------------------------------------
        self.label_width = 150
        self.units_width = 90

        self.scroll_area = QScrollArea()
        # splitter.addWidget(self.scroll)
        # self.cell_def_horiz_layout.addWidget(self.scroll)

        self.custom_data_params = QWidget()
        self.main_layout = QVBoxLayout()
        # self.main_layout.addStretch(0)

        #------------------
        controls_hbox = QHBoxLayout()
        # self.new_button = QPushButton("New")
        self.new_button = QPushButton("Append 5 more rows")
        controls_hbox.addWidget(self.new_button)
        self.new_button.clicked.connect(self.append_more_cb)

        # self.copy_button = QPushButton("Copy")
        # controls_hbox.addWidget(self.copy_button)

        self.clear_button = QPushButton("Clear selected rows")
        controls_hbox.addWidget(self.clear_button)
        self.clear_button.clicked.connect(self.clear_rows_cb)

        #------------------
		# <random_seed type="int" units="dimensionless">0</random_seed> 
		# <cargo_signal_D type="double" units="micron/min^2">1e3</cargo_signal_D>

        # Fixed names for columns:
        hbox = QHBoxLayout()
        # self.select = QCheckBox("")
        col1 = QLabel("Name (required)")
        col1.setAlignment(QtCore.Qt.AlignCenter)
        hbox.addWidget(col1)
        # col2 = QLabel("Type")
        # col2.setAlignment(QtCore.Qt.AlignCenter)
        # hbox.addWidget(col2)
        col3 = QLabel("Default Value (floating point)")
        col3.setAlignment(QtCore.Qt.AlignCenter)
        hbox.addWidget(col3)
        col4 = QLabel("Units")
        col4.setFixedWidth(self.units_width)
        col4.setAlignment(QtCore.Qt.AlignCenter)
        hbox.addWidget(col4)
        # label.setFixedWidth(180)
        self.main_layout.addLayout(hbox)

        #------------------

        # Create lists for the various input boxes
        self.select = []
        self.name = []
        # self.type = []
        self.value = []
        self.units = []
        self.description = []

        # self.type_dropdown = QComboBox()
        # self.type_dropdown.setFixedWidth(300)
        # # self.type_dropdown.currentIndexChanged.connect(self.cycle_changed_cb)
        # self.type_dropdown.addItem("int")
        # self.type_dropdown.addItem("double")
        # self.type_dropdown.addItem("bool")
        # self.type_dropdown.addItem("text")

        for idx in range(self.max_rows):  # TODO: mismatch with # rows in Cell Types|Custom Data may cause problems.
            # self.main_layout.addLayout(NewUserParam(self))
            hbox = QHBoxLayout()
            w = QCheckBox("")
            self.select.append(w)
            hbox.addWidget(w)

            # Need to distinguish between editing an existing name, read in from .xml, and creating a new one.
            w_varname = MyQLineEdit()
            self.name.append(w_varname)
            w_varname.vname = w_varname  # ??
            w_varname.idx = idx
            w_varname.textChanged[str].connect(self.custom_data_name_changed)  # being explicit about passing a string 
            # self.diffusion_coef.enter.connect(self.save_xml)
            hbox.addWidget(w_varname)
            # if idx == 0:
            #     w.setText("random_seed")

            # w = QComboBox()
            # # xml2jupyter: {"double":"FloatText", "int":"IntText", "bool":"Checkbox", "string":"Text", "divider":""}
            # w.addItem("double")
            # w.addItem("int")
            # w.addItem("bool")
            # w.addItem("string")
            # if idx == 0:
            #     w.setCurrentIndex(1)
            # self.type.append(w)
            # hbox.addWidget(w)

            w = MyQLineEdit()
            w.setValidator(QtGui.QDoubleValidator())
            w.setText("0.0")
            w.vname = w_varname
            w.idx = idx
            w.textChanged[str].connect(self.custom_data_value_changed)  # being explicit about passing a string 
            self.value.append(w)
            # if idx == 0:
            #     w.setText("0")
            hbox.addWidget(w)

            w = QLineEdit()
            w.setFixedWidth(self.units_width)
            self.units.append(w)
            hbox.addWidget(w)

            # units = QLabel("micron^2/min")
            # units.setFixedWidth(units_width)
            # hbox.addWidget(units)
            self.main_layout.addLayout(hbox)

            #-----
            hbox = QHBoxLayout()
            w = QLabel("      Description:")
            hbox.addWidget(w)

            w = QLineEdit()
            self.description.append(w)
            hbox.addWidget(w)
            # w.setStyleSheet("background-color: lightgray")
            # w.setStyleSheet("background-color: #e4e4e4")
            self.main_layout.addLayout(hbox)

            #-----
            # self.vbox.addLayout(hbox)
            self.count = self.count + 1
            # print(self.count)


        #==================================================================
        self.custom_data_params.setLayout(self.main_layout)

        self.scroll_area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll_area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll_area.setWidgetResizable(True)

        self.scroll_area.setWidget(self.custom_data_params)

        self.layout = QVBoxLayout(self)

        self.layout.addLayout(controls_hbox)
        warning = QLabel("                                        Note: changing a default value here will also change it in each Cell Type")
        warning.setStyleSheet('color: red;')  # ???background-color: cyan???
        self.layout.addWidget(warning)
        self.layout.addWidget(self.scroll_area)

    # --- custom data (rwh: OMG, this took a lot of time to solve!)
    def custom_data_name_changed(self, text):
        print("--------- (master) custom_data tab: custom_data_name_changed() --------")
        # print("self.sender() = ", self.sender())
        vname = self.sender().vname.text()
        idx = self.sender().idx
        prev_vname = self.celldef_tab.custom_data_name[idx].text()
        print(" prev = ",prev_vname)
        # print("(master) prev_vname = ", self.sender().prev_vname)
        print("(master) vname = ", vname)
        print("(master) idx = ", idx)
        print("(master) custom_data_name_changed(): text = ", text)
        print()

        # Update the value in all Cell Types|Custom Data
        num_cell_defs = self.celldef_tab.tree.invisibleRootItem().childCount()
        print("  num_cell_defs =",num_cell_defs )
        for k in self.celldef_tab.param_d.keys():   # for all cell types
            # mydict[k_new] = mydict.pop(k_old)
            if len(prev_vname) > 0:
                self.celldef_tab.param_d[k]['custom_data'][vname] = self.celldef_tab.param_d[k]['custom_data'].pop(prev_vname)
            else:  # adding a new one
                self.celldef_tab.param_d[k]['custom_data'][vname] = "0.0"
                self.celldef_tab.custom_data_value[idx].setText("0.0")
                self.celldef_tab.custom_data_count += 1
        #     self.celldef_tab.param_d[k]['custom_data'][vname] = text
        #     print(" ===>>> ",k, " : ", self.celldef_tab.param_d[k])

            self.celldef_tab.custom_data_name[idx].setText(vname)
        #     print()

    # --- custom data (rwh: OMG, this took a lot of time to solve!)
    def custom_data_value_changed(self, text):
        print("--------- (master) custom_data tab: custom_data_value_changed() --------")
        # print("self.sender() = ", self.sender())
        vname = self.sender().vname.text()
        idx = self.sender().idx
        print("(master) vname = ", vname)
        print("(master) idx = ", idx)
        print("(master) custom_data_value_changed(): text = ", text)
        print()

        # Update the value in all Cell Types|Custom Data
        num_cell_defs = self.celldef_tab.tree.invisibleRootItem().childCount()
        print("  num_cell_defs =",num_cell_defs )
        for k in self.celldef_tab.param_d.keys():   # for all cell types
            self.celldef_tab.param_d[k]['custom_data'][vname] = text
            # print(" ===>>> ",k, " : ", self.celldef_tab.param_d[k])

            self.celldef_tab.custom_data_value[idx].setText(text)
            # print()


        # populate: self.param_d[cell_def_name]['custom_data'] =  {'cvar1': '42.0', 'cvar2': '0.42', 'cvar3': '0.042'}
        # self.param_d[self.current_cell_def]['custom_data']['cvar1'] = text
        # self.param_d[self.current_cell_def]['custom_data'][vname] = text
        # print(self.param_d[self.current_cell_def]['custom_data'])


    # @QtCore.Slot()
    def clear_rows_cb(self):
        print("----- clearing all selected rows")
        self.custom_vars_to_delete.clear()
        for idx in range(self.count):
            if self.select[idx].isChecked():
                self.custom_vars_to_delete.append(self.name[idx].text())
                self.name[idx].clear()
                # self.value[idx].clear()
                self.value[idx].setText("0.0")
                self.units[idx].clear()
                self.description[idx].clear()
                self.select[idx].setChecked(False)
        print("custom_vars_to_delete = ",self.custom_vars_to_delete)

        # rwh/todo: also update all cell types custom data

    # @QtCore.Slot()
    def append_more_cb(self):
        # print("---- append_more_cb()")
        # self.create_user_param()
        # self.scrollLayout.addRow(NewUserParam(self))
        for idx in range(5):
            # self.main_layout.addLayout(NewUserParam(self))
            hbox = QHBoxLayout()
            w = QCheckBox("")
            self.select.append(w)
            hbox.addWidget(w)

            w = QLineEdit()
            self.name.append(w)
            hbox.addWidget(w)

            w = QLineEdit()
            self.value.append(w)
            w.setValidator(QtGui.QDoubleValidator())
            hbox.addWidget(w)

            w = QLineEdit()
            w.setFixedWidth(self.units_width)
            self.units.append(w)
            hbox.addWidget(w)

            self.main_layout.addLayout(hbox)

            hbox = QHBoxLayout()
            w = QLabel("Desc:")
            hbox.addWidget(w)

            w = QLineEdit()
            self.description.append(w)
            hbox.addWidget(w)
            w.setStyleSheet("background-color: lightgray")


            # units = QLabel("micron^2/min")
            # units.setFixedWidth(units_width)
            # hbox.addWidget(units)
            self.main_layout.addLayout(hbox)

            self.count = self.count + 1
            print(self.count)
    #     # self.text.setText(random.choice(self.hello))
    #     pass

    def clear_gui(self, celldef_tab):
        # pass
        for idx in range(self.count):
            self.name[idx].setText("")
            self.value[idx].setText("0.0")
            self.units[idx].setText("")
            self.description[idx].setText("")
        # self.count = 0

        # self.name.clear()
        # self.value.clear()
        # self.units.clear()
        # self.description.clear()

        celldef_tab.clear_custom_data_tab()


    def fill_gui(self, celldef_tab):   # == populate()
        # pass
        uep_custom_data = self.xml_root.find(".//cell_definitions//cell_definition[1]//custom_data")
        print('fill_gui(): uep_custom_data=',uep_custom_data)

        idx = 0
        # rwh/TODO: if we have more vars than we initially created rows for, we'll need
        # to call 'append_more_cb' for the excess.

        # Should we also update the Cell Types | Custom Data tab entries?

        # for idx in range(self.custom_data_count):
        #     self.custom_data_name[idx].setReadOnly(False)
        #     self.custom_data_name[idx].setText("")
        #     self.custom_data_name[idx].setReadOnly(True)

        #     self.custom_data_value[idx].setText("")

        #     self.custom_data_units[idx].setReadOnly(False)
        #     self.custom_data_units[idx].setText("")
        #     self.custom_data_units[idx].setReadOnly(True)

        for var in uep_custom_data:
            # print(idx, ") ",var)
            self.name[idx].setText(var.tag)
            # print("tag=",var.tag)

            self.value[idx].setText(var.text)

            if 'units' in var.keys():
                self.units[idx].setText(var.attrib['units'])

            if 'description' in var.keys():
                self.description[idx].setText(var.attrib['description'])
            idx += 1
            # print("custom_data:  fill_gui(): idx=",idx,", count=",self.count)
            if idx > self.count:
                self.append_more_cb()

    # doesn't make sense for this tab; custom_data is saved in XML for each cell type
    # def fill_xml(self):
        # pass
    