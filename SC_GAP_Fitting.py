import sys
import json
import numpy as np

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtCore
from PyQt5 import QtGui

import matplotlib
from scipy import optimize
from scipy import interpolate as interp
from scipy.integrate import quad
from sympy import *
from Fitting_Method import *

matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import re

class PlotFrame_SC_GAP:

    def __init__(self):
        self.filename_present = None
        self.raw_data_dict = {}
        self.raw_bias_dict = {}
        self.process_bias_present = None  # 用来存放处理之后的偏压值(即增加能量点或能量范围延展)
        self.process_data_present = None  # 用来存放当前的处理后数据的平均AVG
        self.if_rm_bkgd = False  # 当前是否已经去掉背底倾斜
        self.fitting_methods = Fitting_Methods()
        self.listwidgetitemflags = QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsEnabled | \
                                   QtCore.Qt.ItemIsDragEnabled | QtCore.Qt.ItemIsUserCheckable
        self.current_fitting_type = ""
        self.slider_widget = {}
        self.lineEdit_widget = {}
        self.lower_bound = {}
        self.upper_bound = {}
        self.para_dict = {}
        self.para_fixed = {}
        self.axis = []
        self.setup_mainFrame()
        self.retranslateUi()
        self.setup_mainwindow()

    def setup_mainFrame(self):
        self.mainWindow = QMainWindow()
        self.mainWindow.setMinimumSize(1180, 850)
        self.mainWindow.setWindowTitle('SC_GAP_Fitting')
        self.main_frame = QWidget(self.mainWindow)
        font = QtGui.QFont()
        font.setFamily("Consola")
        font.setPointSize(9)
        if "groupBox_Load_Files":
            self.groupBox_Load_Files = QGroupBox(self.main_frame)
            self.groupBox_Load_Files.setFixedSize(250, 700)
            self.groupBox_Load_Files.setObjectName("groupBox_Single_Pixel")  # 不是界面上显示的
            self.Load_File_btn = QPushButton(self.main_frame)
            self.Load_File_btn.clicked.connect(self._handleOn_LoadFile)
            self.Load_File_btn.setFixedSize(115, 40)
            self.data_process_btn = QPushButton("process_all",self.main_frame)
            self.data_process_btn.clicked.connect(self._handleOn_Data_Process_All)
            self.data_process_btn.setFixedSize(115, 40)
            self.import_all_btn = QPushButton("import_all",self.main_frame)
            self.import_all_btn.clicked.connect(self._handleOn_ImportAllFiguresFromDataList)
            self.import_all_btn.setFixedSize(115, 40)
            self.delete_all_btn = QPushButton("delete_all",self.main_frame)
            self.delete_all_btn.clicked.connect(self._handleOn_DeleteAllFromDataList)
            self.delete_all_btn.setFixedSize(115, 40)
            HBox1 = QHBoxLayout()
            HBox2 = QHBoxLayout()
            HBox1.addWidget(self.Load_File_btn)
            HBox1.addWidget(self.delete_all_btn)
            HBox2.addWidget(self.data_process_btn)
            HBox2.addWidget(self.import_all_btn)

            self.data_QList = QListWidget(self.main_frame)
            self.data_QList.setFixedSize(230, 550)
            self.data_QList.itemDoubleClicked.connect(self._handleOnImportFigureFromDataList)
            self.data_QList.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
            data_QList_menu = QMenu(self.data_QList)
            self.data_QList_menu_import_figure_action = QAction(data_QList_menu)
            self.data_QList_menu_delete_data_action = QAction(data_QList_menu)
            self.data_QList_menu_import_figure_action.triggered.connect(self._handleOnImportFigureFromDataList)
            self.data_QList_menu_delete_data_action.triggered.connect(self._handleOnDeleteDataFromDataList)
            data_QList_menu.addAction(self.data_QList_menu_import_figure_action)
            data_QList_menu.addAction(self.data_QList_menu_delete_data_action)
            self.data_QList_menu_import_figure_action.setText("import figure")
            self.data_QList_menu_delete_data_action.setText("delete item")
            def data_QList_menu_show():
                if self.data_QList.currentItem() is None:
                    return
                data_QList_menu.exec_(QtGui.QCursor.pos())

            self.data_QList.customContextMenuRequested.connect(data_QList_menu_show)

            groupBox_Load_Files_Layout = QVBoxLayout(self.groupBox_Load_Files)
            groupBox_Load_Files_Layout.addLayout(HBox1, QtCore.Qt.AlignCenter)
            groupBox_Load_Files_Layout.addLayout(HBox2, QtCore.Qt.AlignCenter)
            groupBox_Load_Files_Layout.addWidget(self.data_QList)

        if "groupBox_Single_Pixel":
            self.groupBox_Single_Pixel = QGroupBox(self.main_frame)
            self.groupBox_Single_Pixel.setFixedSize(1150, 1150)
            self.groupBox_Single_Pixel.setObjectName("groupBox_Single_Pixel")

            self.button_remove_all = QPushButton("REMOVE ALL", self.groupBox_Single_Pixel)
            self.button_remove_all.setFixedSize(120, 32)
            self.button_remove_all.clicked.connect(self._handleOn_Remove_All)

            self.button_insert = QPushButton("INSERT", self.groupBox_Single_Pixel)
            self.button_insert.setFixedSize(120, 32)
            self.button_insert.clicked.connect(self._handleOn_Insert)

            self.fig = plt.figure(figsize=(16, 10), dpi=100)
            self.axes = plt.axes([0.10, 0.25, 0.85, 0.70])
            self.axes_error = plt.axes([0.10, 0.05, 0.85, 0.15])
            self.canvas = FigureCanvas(self.fig)
            self.canvas.setParent(self.groupBox_Single_Pixel)
            self.mpl_toolbar_origin = NavigationToolbar(self.canvas, self.groupBox_Single_Pixel)

            HBox = QHBoxLayout()
            HBox.addWidget(self.mpl_toolbar_origin)
            HBox.addWidget(self.button_insert)
            HBox.addWidget(self.button_remove_all)
            groupBox_Single_Pixel_Layout = QVBoxLayout(self.groupBox_Single_Pixel)
            groupBox_Single_Pixel_Layout.addLayout(HBox)
            groupBox_Single_Pixel_Layout.addWidget(self.canvas)

        if "groupBox_Fitting_method":
            self.groupBox_Fitting_method = QGroupBox(self.main_frame)
            self.groupBox_Fitting_method.setFixedSize(600, 180)

            self.checkBox_one_band = QCheckBox(self.groupBox_Fitting_method)
            self.checkBox_one_band.setGeometry(QtCore.QRect(10, 20, 130, 52))
            self.checkBox_one_band.setFont(font)
            self.comboBox_one_band = QComboBox(self.groupBox_Fitting_method)
            self.comboBox_one_band.setGeometry(QtCore.QRect(10, 70, 100, 38))
            self.comboBox_one_band.setFont(font)
            self.comboBox_one_band.addItem("s")
            self.comboBox_one_band.addItem("2")
            self.comboBox_one_band.addItem("4")
            self.comboBox_one_band.addItem("6")
            self.comboBox_one_band.addItem("s+2")
            self.comboBox_one_band.addItem("s+4")
            self.comboBox_one_band.addItem("s+6")
            self.comboBox_one_band.addItem("2+4")
            self.comboBox_one_band.addItem("2+6")
            self.comboBox_one_band.addItem("4+6")

            self.checkBox_two_band = QCheckBox(self.groupBox_Fitting_method)
            self.checkBox_two_band.setGeometry(QtCore.QRect(175, 20, 138, 52))
            self.checkBox_two_band.setFont(font)
            self.comboBox_two_band = QComboBox(self.groupBox_Fitting_method)
            self.comboBox_two_band.setGeometry(QtCore.QRect(175, 70, 100, 38))
            self.comboBox_two_band.addItem("2s")
            self.comboBox_two_band.addItem("s+2-two-band")
            self.comboBox_two_band.addItem("s+4-two-band")
            self.comboBox_two_band.addItem("s+6-two-band")
            self.comboBox_two_band.addItem("2+4-two-band")
            self.comboBox_two_band.addItem("2+6-two-band")
            self.comboBox_two_band.addItem("4+6-two-band")

            self.pushButton_fitting = QPushButton(self.groupBox_Fitting_method)
            self.pushButton_fitting.setGeometry(QtCore.QRect(95, 135, 105, 38))
            self.pushButton_fitting.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
            self.pushButton_fitting.clicked.connect(self._handleOn_Fitting_by_Methods)
            self.label_parameter_linewidth = QLabel("linewidth",self.groupBox_Fitting_method)
            self.label_parameter_linewidth.setGeometry(QtCore.QRect(220, 135, 100, 38))
            self.label_parameter_linewidth.setFont(font)
            self.lineEdit_parameter_linewidth = QLineEdit(self.groupBox_Fitting_method)
            self.lineEdit_parameter_linewidth.setGeometry(QtCore.QRect(325, 135, 50, 30))

        if "groupBox_Fitting_Parameters":
            self.groupBox_Fitting_Parameters = QGroupBox(self.main_frame)
            self.groupBox_Fitting_Parameters.setFixedSize(700, 600)
            self.groupBox_Fitting_Parameters.setObjectName("groupBox_Fitting_Parameters")  # setObjectName有什么作用??
            self.btn_fitting = QPushButton("fit", self.groupBox_Fitting_Parameters)
            self.btn_fitting.setFixedSize(100, 30)
            self.btn_fitting.clicked.connect(self._handleOn_Fitting_Fixed_Para)
            font.setPointSize(11)
            def setup_para_widget(para_name:str):
                para_label = QLabel(para_name, self.groupBox_Fitting_Parameters)
                para_label.setFixedSize(70, 30)
                para_label.setFont(font)
                para_linedit = QLineEdit(self.groupBox_Fitting_Parameters)
                para_linedit.setFixedSize(100, 30)
                para_slider = QSlider(self.groupBox_Fitting_Parameters)
                para_slider.setFixedSize(175, 35)
                para_slider.setOrientation(QtCore.Qt.Horizontal)
                para_slider.setMinimum(0)  # 设置最小值
                para_slider.setMaximum(1000)  # 设置最大值
                para_slider.setSingleStep(2)  # 设置步长值
                para_slider.setValue(0)  # 设置当前值
                para_min = QLineEdit(self.groupBox_Fitting_Parameters)
                para_min.setFixedSize(50, 35)
                para_max = QLineEdit(self.groupBox_Fitting_Parameters)
                para_max.setFixedSize(50, 35)
                self.slider_widget[para_name] = para_slider
                self.lineEdit_widget[para_name] = para_linedit
                self.lower_bound[para_name] = para_min
                self.upper_bound[para_name] = para_max
                HBOX_1 = QHBoxLayout()
                HBOX_1.addWidget(para_label)
                HBOX_1.addWidget(para_linedit)
                HBOX_2 = QHBoxLayout()
                HBOX_2.addWidget(para_min)
                HBOX_2.addWidget(para_slider)
                HBOX_2.addWidget(para_max)
                VBOX = QVBoxLayout()
                VBOX.addLayout(HBOX_1)
                VBOX.addLayout(HBOX_2)
                return VBOX

            if "para":
                VBOX_O = setup_para_widget("x")
                self.slider_widget["x"].valueChanged.connect(self._handleOn_delta0_Change)
                self.slider_widget["x"].valueChanged.connect(self._handleOn_x_Change)
                VBOX_A = setup_para_widget("delta0")
                self.slider_widget["delta0"].valueChanged.connect(self._handleOn_delta0_Change)
                VBOX_B = setup_para_widget("delta1")
                self.slider_widget["delta1"].valueChanged.connect(self._handleOn_delta1_Change)
                VBOX_C = setup_para_widget("eff_T1")
                self.slider_widget["eff_T1"].valueChanged.connect(self._handleOn_effective_Temperature1_Change)
                VBOX_D = setup_para_widget("gamma1")
                self.slider_widget["gamma1"].valueChanged.connect(self._handleOn_gamma1_Change)
                VBOX_I = setup_para_widget("phi")
                self.slider_widget["phi"].valueChanged.connect(self._handleOn_phi_Change)

                VBOX_E = setup_para_widget("delta2")
                self.slider_widget["delta2"].valueChanged.connect(self._handleOn_delta2_Change)
                VBOX_F = setup_para_widget("eff_T2")
                self.slider_widget["eff_T2"].valueChanged.connect(self._handleOn_effective_Temperature2_Change)
                VBOX_G = setup_para_widget("gamma2")
                self.slider_widget["gamma2"].valueChanged.connect(self._handleOn_gamma2_Change)
                VBOX_H = setup_para_widget("offset")
                self.slider_widget["offset"].valueChanged.connect(self._handleOn_offset_change)

            groupBox_Fitting_Parameters_Layout = QGridLayout(self.groupBox_Fitting_Parameters)
            groupBox_Fitting_Parameters_Layout1 = QVBoxLayout()
            groupBox_Fitting_Parameters_Layout1.addLayout(VBOX_A)
            groupBox_Fitting_Parameters_Layout1.addLayout(VBOX_B)
            groupBox_Fitting_Parameters_Layout1.addLayout(VBOX_C)
            groupBox_Fitting_Parameters_Layout1.addLayout(VBOX_D)
            groupBox_Fitting_Parameters_Layout1.addLayout(VBOX_I)
            groupBox_Fitting_Parameters_Layout2 = QVBoxLayout()
            groupBox_Fitting_Parameters_Layout2.addLayout(VBOX_E)
            groupBox_Fitting_Parameters_Layout2.addLayout(VBOX_F)
            groupBox_Fitting_Parameters_Layout2.addLayout(VBOX_G)
            self.groupBox_band_one = QGroupBox("band_one",self.groupBox_Fitting_Parameters)
            self.groupBox_band_one.setLayout(groupBox_Fitting_Parameters_Layout1)
            self.groupBox_band_two = QGroupBox("band_two",self.groupBox_Fitting_Parameters)
            self.groupBox_band_two.setLayout(groupBox_Fitting_Parameters_Layout2)
            groupBox_Fitting_Parameters_Layout.addWidget(self.btn_fitting,0,0,1,1)
            groupBox_Fitting_Parameters_Layout.addLayout(VBOX_O,1, 0, 1, 1, QtCore.Qt.AlignTop)
            groupBox_Fitting_Parameters_Layout.addLayout(VBOX_H,1, 1, 1, 2, QtCore.Qt.AlignTop)
            groupBox_Fitting_Parameters_Layout.addWidget(self.groupBox_band_one,2, 0, 1, 1, QtCore.Qt.AlignTop)
            groupBox_Fitting_Parameters_Layout.addWidget(self.groupBox_band_two,2, 1, 1, 1, QtCore.Qt.AlignTop)
            self.groupBox_Fitting_Parameters.setLayout(groupBox_Fitting_Parameters_Layout)

        if "groupBox_Data_PreProcess":
            self.groupBox_Data_PreProcess = QGroupBox(self.main_frame)
            self.groupBox_Data_PreProcess.setFixedSize(250, 300)
            self.groupBox_Data_PreProcess.setObjectName("groupBox_Data_PreProcess")  # setObjectName有什么作用??

            font.setPointSize(9)
            self.label_parameter_amplifier = QLabel(self.groupBox_Data_PreProcess)
            self.label_parameter_amplifier.setFixedSize(100, 32)
            self.label_parameter_amplifier.setFont(font)
            self.text_parameter_amplifier = QLineEdit(self.groupBox_Data_PreProcess)
            self.text_parameter_amplifier.setFixedSize(65, 32)
            self.text_parameter_amplifier.setFont(font)
            self.text_parameter_amplifier.editingFinished.connect(self._handleOn_Amplification)

            self.label_parameter_exp_Energy = QLabel(self.groupBox_Single_Pixel)
            self.label_parameter_exp_Energy.setFixedSize(250, 32)
            self.label_parameter_exp_Energy.setFont(font)

            self.label_parameter_delta_Energy = QLabel(self.groupBox_Data_PreProcess)
            self.label_parameter_delta_Energy.setFixedSize(100, 32)
            self.label_parameter_delta_Energy.setFont(font)
            self.lineEdit_parameter_new_delta_Energy = QLineEdit(self.groupBox_Data_PreProcess)
            self.lineEdit_parameter_new_delta_Energy.setFixedSize(65, 32)
            # self.lineEdit_parameter_new_delta_Energy.editingFinished.connect(self._handleOn_DataChange)

            self.label_parameter_background = QLabel(self.groupBox_Data_PreProcess)
            self.label_parameter_background.setFixedSize(80, 32)
            self.label_parameter_background.setFont(font)
            self.btn_bkgd_rm = QPushButton(self.groupBox_Data_PreProcess)
            self.btn_bkgd_rm.setFixedSize(120, 32)
            self.btn_bkgd_rm.setFont(font)
            self.btn_bkgd_rm.clicked.connect(self._handleOn_Remove_BackGround)
            self.label_parameter_symmetry = QLabel("symmetry",self.groupBox_Data_PreProcess)
            self.label_parameter_symmetry.setFixedSize(80, 32)
            self.label_parameter_symmetry.setFont(font)
            self.btn_symmetry = QPushButton("symm_process",self.groupBox_Data_PreProcess)
            self.btn_symmetry.setFixedSize(120, 32)
            self.btn_symmetry.setFont(font)
            # self.btn_symmetry.clicked.connect(self._handleOn_Symmetry_process)

            HBOX_1 = QHBoxLayout()
            HBOX_1.addWidget(self.label_parameter_amplifier)
            HBOX_1.addWidget(self.text_parameter_amplifier)
            HBOX_3 = QHBoxLayout()
            HBOX_3.addWidget(self.label_parameter_delta_Energy)
            HBOX_3.addWidget(self.lineEdit_parameter_new_delta_Energy)
            HBOX_4 = QHBoxLayout()
            HBOX_4.addWidget(self.label_parameter_background)
            HBOX_4.addWidget(self.btn_bkgd_rm)
            HBOX_5 = QHBoxLayout()
            HBOX_5.addWidget(self.label_parameter_symmetry)
            HBOX_5.addWidget(self.btn_symmetry)
            groupBox_Data_PreProcess_Layout = QVBoxLayout(self.groupBox_Data_PreProcess)
            groupBox_Data_PreProcess_Layout.addLayout(HBOX_1)
            # groupBox_Data_PreProcess_Layout.addLayout(HBOX_2)
            groupBox_Data_PreProcess_Layout.addWidget(self.label_parameter_exp_Energy)
            groupBox_Data_PreProcess_Layout.addLayout(HBOX_3)
            groupBox_Data_PreProcess_Layout.addLayout(HBOX_4)
            groupBox_Data_PreProcess_Layout.addLayout(HBOX_5)
            self.groupBox_Fitting_Parameters.setLayout(groupBox_Data_PreProcess_Layout)

        if "groupBox_fit_lines_list":
            self.groupBox_lines_list = QGroupBox(self.main_frame)
            self.groupBox_lines_list.setFixedSize(250, 355)
            self.save_all_btn = QPushButton("save_all",self.groupBox_lines_list)
            self.save_all_btn.clicked.connect(self._handleOnSave_All_LineListItem)
            self.save_all_btn.setFixedSize(115, 40)
            self.lines_QList = QListWidget(self.groupBox_lines_list)
            self.lines_QList.setFixedSize(230, 250)
            self.lines_QList.itemDoubleClicked.connect(self._handleOnImportFittingLinesFromLinesList)
            self.lines_QList.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
            lines_QList_menu = QMenu(self.data_QList)
            self.lines_QList_menu_import_figure_action = QAction(lines_QList_menu)
            self.lines_QList_menu_delete_data_action = QAction(lines_QList_menu)
            self.lines_QList_menu_hide_figure_action = QAction(lines_QList_menu)
            self.lines_QList_menu_save_data_action = QAction(lines_QList_menu)
            self.lines_QList_menu_import_figure_action.triggered.connect(self._handleOnImportFittingLinesFromLinesList)
            self.lines_QList_menu_delete_data_action.triggered.connect(self._handleOnDeleteDataFromLineList)
            self.lines_QList_menu_save_data_action.triggered.connect(self._handleOnSaveLineList)
            lines_QList_menu.addAction(self.lines_QList_menu_import_figure_action)
            lines_QList_menu.addAction(self.lines_QList_menu_delete_data_action)
            lines_QList_menu.addAction(self.lines_QList_menu_save_data_action)
            self.lines_QList_menu_import_figure_action.setText("import figure")
            self.lines_QList_menu_delete_data_action.setText("delete item")
            self.lines_QList_menu_save_data_action.setText("save item")
            def lines_QList_menu_show():
                if self.lines_QList.currentItem() is None:
                    return
                lines_QList_menu.exec_(QtGui.QCursor.pos())
            self.lines_QList.customContextMenuRequested.connect(lines_QList_menu_show)

            groupBox_lines_QList_Layout = QVBoxLayout(self.groupBox_lines_list)
            groupBox_lines_QList_Layout.addWidget(self.save_all_btn)
            groupBox_lines_QList_Layout.addWidget(self.lines_QList)

        if "mainlayout":
            VBOX_1 = QVBoxLayout()
            VBOX_1.addWidget(self.groupBox_Data_PreProcess)
            VBOX_1.addWidget(self.groupBox_Load_Files)
            VBOX_2 = QVBoxLayout()
            VBOX_2.addWidget(self.groupBox_Fitting_method)
            VBOX_2.addWidget(self.groupBox_lines_list)
            VBOX_2.addWidget(self.groupBox_Fitting_Parameters)
            main_frame_layout = QGridLayout(self.main_frame)
            main_frame_layout.setAlignment(QtCore.Qt.AlignTop)
            main_frame_layout.addLayout(VBOX_1, 0, 0, 2, 1, QtCore.Qt.AlignTop)
            main_frame_layout.addWidget(self.groupBox_Single_Pixel, 0, 1, 2, 2)
            main_frame_layout.addLayout(VBOX_2, 0, 3, 2, 1, QtCore.Qt.AlignTop)
            # main_frame_layout.addWidget(self.scroll_area)
            self.main_frame.setLayout(main_frame_layout)

            self.scroll_area = QScrollArea(self.mainWindow)
            self.scroll_area.setWidget(self.main_frame)
            self.mainWindow.setCentralWidget(self.scroll_area)

    def setup_mainwindow(self):
        # mainWindow的menuBar的设置
        def create_action(text, slot=None, shortcut=None,
                          icon=None, tip=None, checkable=False,
                          signal="triggered()"):  # 创建一个动作,
            action = QAction(text, self.mainWindow)
            if icon is not None:
                action.setIcon(QIcon(":/%s.png" % icon))
            if shortcut is not None:
                action.setShortcut(shortcut)
            if tip is not None:
                action.setToolTip(tip)
                action.setStatusTip(tip)
            if slot is not None:
                action.triggered.connect(slot)
            if checkable:
                action.setCheckable(True)
            return action

        def add_actions(target, actions):
            for action in actions:
                if action is None:
                    target.addSeparator()  # 两个动作之间画一条横线
                else:
                    target.addAction(action)

        def save_plot():  # file_menu的功能之一,用来保存当前的图像,并在statusBar中显示保存的路径
            file_choices = "PNG (*.png)|*.png"

            path, filetype = QFileDialog.getSaveFileName(self.mainWindow,
                                                         'Save file', '',
                                                         file_choices)
            if path:
                self.canvas.print_figure(path, dpi=100)
                self.mainWindow.statusBar().showMessage('Saved to %s' % path, 8000)  # 8000是指message显示的时间长度为8000毫秒

        file_menu = self.mainWindow.menuBar().addMenu("File")  # 在menuBar上添加File菜单
        save_file_action = create_action("Save plot", slot=save_plot, shortcut="Ctrl+S", tip="Save the plot")
        quit_action = create_action("Quit", slot=self.mainWindow.close, shortcut="Ctrl+Q", tip="Close the application")
        add_actions(file_menu, (save_file_action, None, quit_action))

        def on_about():  # help_menu中的功能,用来介绍该界面已经实现的功能
            msg = """ A SC_gap Fitting GUI of using PyQt with matplotlib:

                * Use the matplotlib navigation bar
                * Add values to the text box and press Enter (or click "Draw")
                * Show or hide the grid
                * Drag the slider to modify the width of the bars
                * Save the plot to a file using the File menu
                * Click on a bar to receive an informative message
               """
            QMessageBox.about(self.mainWindow, "About the demo", msg.strip())  # #

        help_menu = self.mainWindow.menuBar().addMenu("Help")
        about_action = create_action("About", slot=on_about, shortcut='F1', tip='About the demo')
        add_actions(help_menu, (about_action,))

        # mainWindow的status Bar的设定
        self.mainWindow.statusBar().showMessage("")

    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        if "groupBox_Load_Files":
            self.groupBox_Load_Files.setTitle(
                _translate("PlotFrame_groupBox_Load_Files_title", "Load_Files"))
            self.Load_File_btn.setText(_translate("PlotFrame_groupBox_Load_Files_text", "LOAD"))

        if "groupBox_Lines_List":
            self.groupBox_lines_list.setTitle(
                _translate("PlotFrame_groupBox_lines_list_title", "Lines_List"))

        if "groupBox_Single_Pixel":
            self.groupBox_Single_Pixel.setTitle(
                _translate("PlotFrame_groupBox_Single_Pixel_title", "Single_Pixel"))

        if "groupBox_Fitting_method":
            self.groupBox_Fitting_Parameters.setTitle(
                _translate("PlotFrame_groupBox_Fitting_method_title", "Fitting_Parameters"))
            self.pushButton_fitting.setText(_translate("PlotFrame_SC_GAP_pushButton_fitting_text", "拟合"))
            self.comboBox_one_band.setToolTip(
                _translate("PlotFrame_SC_GAP_comboBox_one_band_tip", "选择一种各项同性的拟合方式"))
            self.checkBox_one_band.setToolTip(_translate("PlotFrame_SC_GAP_checkBox_one_band_tip", "单带拟合"))
            self.checkBox_one_band.setText(_translate("PlotFrame_SC_GAP_checkBox_one_band_text", "ONE-BAND"))
            self.checkBox_two_band.setToolTip(_translate("PlotFrame_SC_GAP_checkBox_two_band_tip", "双带拟合"))
            self.checkBox_two_band.setText(_translate("PlotFrame_SC_GAP_checkBox_two_band_text", "TWO-BAND"))
            self.comboBox_two_band.setToolTip(
                _translate("PlotFrame_SC_GAP_comboBox_two_band_tip", "选择一种各项异性的拟合方式"))

        if "groupBox_Fitting_Parameters":
            self.groupBox_Fitting_method.setTitle(
                _translate("PlotFrame_SC_GAP_groupBox_Fitting_method_title", "Fitting_Methods"))

        if "groupBox_Data_PreProcess":
            self.groupBox_Data_PreProcess.setTitle(
                _translate("PlotFrame_SC_GAP_groupBox_Data_PreProcess_title", "Data_PreProcess"))
            self.label_parameter_amplifier.setToolTip(
                _translate("PlotFrame_SC_GAP_label_parameter_amplifier_tip", "信号放大器放大倍数"))
            self.label_parameter_amplifier.setText(
                _translate("PlotFrame_SC_GAP_label_parameter_amplifier_text", "AMPLIFIER:"))
            self.label_parameter_delta_Energy.setText(
                _translate("PlotFrame_SC_GAP_label_parameter_delta_Energy_text", "DELTA_E:"))
            self.label_parameter_exp_Energy.setText(
                _translate("PlotFrame_SC_GAP_label_parameter_exp_Energy_text", "ENERGY:"))
            self.label_parameter_exp_Energy.setText(
                _translate("PlotFrame_SC_GAP_label_parameter_exp_Energy_text", "ENERGY:"))
            self.label_parameter_background.setText(
                _translate("PlotFrame_SC_GAP_label_parameter_background_text", "Background:"))
            self.btn_bkgd_rm.setText(
                _translate("PlotFrame_SC_GAP_self.btn_bkgd_rm_text", "rm bkgd"))

    def parameter_change(self, para_name:str):
        para_dict = self.getParameters_bound_FromInput()
        step = (para_dict[para_name][1] - para_dict[para_name][0]) / \
               (self.slider_widget[para_name].maximum() - self.slider_widget[para_name].minimum())
        para = self.slider_widget[para_name].value() * step + para_dict[para_name][0]
        self.lineEdit_widget[para_name].setText(str(para))
        return para

    def _handleOn_x_Change(self):
        try:
            x = self.parameter_change("x")
            self.fitting_family({"x":x})
        except Exception as e:
            self.informMsg(str(e))
            return

    def _handleOn_delta0_Change(self):
        try:
            delta0 = self.parameter_change("delta0")
            self.fitting_family({"delta0":delta0})
        except Exception as e:
            self.informMsg(str(e))
            return

    def _handleOn_phi_Change(self):
        try:
            phi = self.parameter_change("phi")
            self.fitting_family({"phi":phi})
        except Exception as e:
            self.informMsg(str(e))
            return

    def _handleOn_delta1_Change(self):
        try:
            delta1 = self.parameter_change("delta1")
            self.fitting_family({"delta1":delta1})
        except Exception as e:
            self.informMsg(str(e))
            return

    def _handleOn_gamma1_Change(self):
        try:
            gamma1 = self.parameter_change("gamma1")
            self.fitting_family({"gamma1":gamma1})
        except Exception as e:
            self.informMsg(str(e))
            return

    def _handleOn_effective_Temperature1_Change(self):
        try:
            eff_T1 = self.parameter_change("eff_T1")
            self.fitting_family({"eff_T1":eff_T1})
        except Exception as e:
            self.informMsg(str(e))
            return

    def _handleOn_delta2_Change(self):
        try:
            delta2 = self.parameter_change("delta2")
            self.fitting_family({"delta2":delta2})
        except Exception as e:
            self.informMsg(str(e))
            return

    def _handleOn_gamma2_Change(self):
        try:
            gamma2 = self.parameter_change("gamma2")
            self.fitting_family({"gamma2":gamma2})
        except Exception as e:
            self.informMsg(str(e))
            return

    def _handleOn_effective_Temperature2_Change(self):
        try:
            eff_T2 = self.parameter_change("eff_T2")
            self.fitting_family({"eff_T2":eff_T2})
        except Exception as e:
            self.informMsg(str(e))
            return

    def _handleOn_offset_change(self):
        try:
            offset = self.parameter_change("offset")
            bias_new = self.process_bias_present + offset
            fit_type = self.current_fitting_type  # 只有进行拟合之后self.current_fitting_type才会有具体的值
            item = self._getItemFromName(self.lines_QList, fit_type)
            fit_data = item.data(1)
            item.data(4).remove()
            self.canvas_plot_from_data(bias_new, fit_data, item)
        except Exception as e:
            self.informMsg(str(e) + ":offset失败")
            return False

    # 得到一个QWidgetList的item
    def _getItemFromName(self, parent:QListWidget, name) -> QListWidgetItem:
        row = parent.count()
        for i in range(row):
            if parent.item(i).text() == name:
                return parent.item(i)
        item = QListWidgetItem(parent)
        item.setText(name)
        return item

    def _handleOn_LoadFile(self) -> bool:
        if "get fileName":
            fileNames_List, fileType = QFileDialog.getOpenFileNames(self.main_frame, r'Load json',
                                                             r'D:\Users\yg\PycharmProjects\spectra_data',
                                                             r'dat Files(*.dat)')
            print(fileNames_List)
            for fileName in fileNames_List:
                if fileName == "":
                    self.informMsg("未选择文件")
                    return False

                fileName_list = fileName.split('/')
                filename_new = fileName_list[-1]
                if filename_new in self.raw_bias_dict.keys():
                    reply = self.questionMsg("List中已经存在相同名称,是否进行覆盖？")
                    if reply == False:
                        return False

                if "get bias and DIDV_AVG":
                    with open(fileName, 'r') as of:
                        N = -1  # 用N记录当前读取的行号(空行并不算在内)
                        of.seek(0, 2)  # 指针移动到文件末尾
                        size = of.tell()  # 此时指针的位置即为文件末尾的位置
                        of.seek(0, 0)  # 把指针重新移回文件开头

                        # 开始读取dat文件中的数据
                        read_line = of.readline()
                        read_line = read_line.strip()
                        while read_line != '[DATA]':
                            if of.tell() >= size:
                                self.informMsg("END OF FILE")
                                break
                            read_line = of.readline()
                            read_line = read_line.strip()
                            if read_line.strip() == "":  # 读到的是空行
                                continue
                            N += 1
                        # 读取频道
                        label_line = of.readline()  # 字符串
                        label_line = label_line.strip()  # 仍为字符串
                        label_list = label_line.split('	')  # 列表list
                        AVG_index_list = []
                        for i in range(len(label_list)):
                            if label_list[i].find("LI Demod 1 X [AVG]") != -1:
                                AVG_index_list.append(i)
                        # 读取数据
                        data_lines = of.readlines()  # a list of strings
                        raw_data_matrix = []
                        for line in data_lines:
                            line_new = line.strip('\n')
                            line_new = line_new.split('\t')
                            raw_data_matrix.append(line_new)  # a list of lists of strings
                        raw_data_matrix = np.array(raw_data_matrix, dtype="float64")  # ndarray
                        dim = raw_data_matrix.shape  # 例如(301,45)
                        # 获取偏压
                        bias_array = np.zeros((1, dim[0]), dtype='float64')
                        bias_array += (raw_data_matrix.T)[0]  # shape例如:(1,301)
                        # 获取原始DIDV数据的平均
                        raw_avg_DIDV_data = np.zeros((1, dim[0]), dtype='float64')  # shape例如:(1,301)
                        for index in AVG_index_list:
                            raw_avg_DIDV_data = np.add(raw_avg_DIDV_data, (raw_data_matrix.T)[index])
                        raw_avg_DIDV_data = np.true_divide(raw_avg_DIDV_data,
                                                           np.mean(raw_avg_DIDV_data[0][0:10]))  # 高能量处的DIDV进行归一化,为什么不是全部??

                        self.filename_present = fileName_list[-1]
                        origin_bias_present = (bias_array[0])[::-1]  # 这样才是真正的一维数组
                        raw_data_present = (raw_avg_DIDV_data[0])[::-1]  # 这样才是真正的一维数组
                        self.raw_bias_dict[self.filename_present] = origin_bias_present
                        self.raw_data_dict[self.filename_present] = raw_data_present
                        self.process_bias_present = origin_bias_present  # mV为单位
                        self.process_data_present = raw_data_present  # DIDV数据不改动
                        self.if_rm_bkgd = False

                if "add to list and set some QWidget_Texts":
                    row = 0
                    while row < self.data_QList.count():
                        if self.data_QList.item(row).text() == self.filename_present:
                            break
                        row += 1
                    if row != self.data_QList.count():  # 前面已经问过相同的问题了,因此这里直接删除
                        self.data_QList.takeItem(row)

                    item = self._getItemFromName(parent=self.data_QList, name=self.filename_present)
                    bias_data_ori = np.vstack((origin_bias_present, raw_data_present))
                    bias_data_pro = np.vstack((self.process_bias_present, self.process_data_present))
                    item.setData(-1, bias_data_ori)
                    item.setData(1, bias_data_pro)
                    item.setFlags(self.listwidgetitemflags)
                    self.data_QList.addItem(item)
                    self.data_QList.sortItems()
                    self.data_QList.setCurrentItem(item)
                    E_min = format(self.process_bias_present[0], '.3f')
                    E_max = format(self.process_bias_present[-1], '.3f')
                    delta_E = format(self.process_bias_present[2] - self.process_bias_present[1], '.3f')  # 偏压是等差数列
                    self.label_parameter_exp_Energy.setText(
                        "ENE_SCALE:" + str(E_min) + ":" + str(delta_E) + ":" + str(E_max))
                    self.lineEdit_parameter_new_delta_Energy.setText(str(delta_E))
                    # self.lineEdit_parameter_offset.setText(str(format(DIDV_min2bias, '.3f')))
            return True

    def _handleOn_ImportAllFiguresFromDataList(self):
        row = self.data_QList.count()
        if row == 0:
            self.informMsg("data_list中無item")
            return
        for i in range(row):
            item = self.data_QList.item(i)
            try:
                item.data(4).remove()
            except Exception as e:
                self.informMsg(str(e))
            try:
                bias = item.data(1)[0]
                data = item.data(1)[1]
                self.canvas_plot_from_data(bias, data, item)
            except Exception as e:
                self.informMsg(str(e) + ":import失败")

    def _handleOnImportFigureFromDataList(self):  # 用来作原始曲线
        item = self.data_QList.currentItem()
        try:
            item.data(4).remove()
        except Exception as e:
            self.informMsg(str(e))
        try:
            self.filename_present = item.text()
            self.process_bias_present = item.data(1)[0]
            self.process_data_present = item.data(1)[1]
            self.if_rm_bkgd = False  # 当前是否已经去掉背底倾斜
            bias = item.data(1)[0]
            data = item.data(1)[1]
            self.canvas_plot_from_data(bias, data, item)
        except Exception as e:
            self.informMsg(str(e) + ":import失败")

    def _handleOnImportFittingLinesFromLinesList(self):  # 都是self.filename_present下的
        item = self.lines_QList.currentItem()
        try:
            item.data(4).remove()
        except Exception as e:
            self.informMsg(str(e))
        try:
            para = item.data(-1)
            bias = self.process_bias_present
            data = item.data(1)
            self.canvas_plot_from_data(bias, data, item)
            self.current_fitting_type = item.text()
        except Exception as e:
            self.informMsg(str(e) + ":import失败")
        return

    def find_axis(self, bias, data):
        bias_max = max(bias)
        bias_min = min(bias)
        data_min = min(data)
        data_max = max(data)
        if self.axis == []:
            self.axis = [bias_min,bias_max,data_min,data_max]
        else:
            if self.axis[0] > bias_min:
                self.axis[0] = bias_min
            if self.axis[1] < bias_max:
                self.axis[1] = bias_max
            if self.axis[2] > data_min:
                self.axis[2] = data_min
            if self.axis[3] > data_max:
                self.axis[3] = data_max

    def canvas_plot_from_data(self, bias, data, item:QListWidgetItem):
        try:
            row = self.lines_QList.count()
            for i in range(row):
                if item == self.lines_QList.item(i):
                    try:
                        linewidth = float(self.lineEdit_parameter_linewidth.text())
                    except:
                        linewidth = 1.5
                    line = self.axes.plot(bias, data, linestyle='-', linewidth=linewidth)[0]
                    item.setData(4, line)
                    item_ori = self.data_QList.currentItem()
                    try:
                        self.axes_error.cla()
                        self.axes_error.plot(bias, (data - item_ori.data(-1)[1]) / item_ori.data(-1)[1], linestyle='-',linewidth=1.5)
                    except:
                        self.axes_error.plot(bias, (data - item_ori.data(-1)[1]) / item_ori.data(-1)[1], linestyle='-', linewidth=1.5)
                    self.canvas.draw()
            row = self.data_QList.count()
            for i in range(row):
                if item == self.data_QList.item(i):
                    line = self.axes.plot(bias, data, linestyle='-', linewidth=2)[0]
                    item.setData(4, line)
                    self.canvas.draw()
        except Exception as e:
            self.informMsg(str(e) + ":axis fail")

    def _handleOn_Insert(self):
        # ax_inset = self.fig.subplot()
        # ax_inset = self.axes.inset_axes([0.71, 0.01, 0.28, 0.3])
        try:
            self.ax_inset.cla()
        except Exception as e:
            self.ax_inset = self.fig.add_axes([0.71, 0.30, 0.2, 0.2], polar=True)
        try:
            fit_type = self.current_fitting_type  # 只有进行拟合之后self.current_fitting_type才会有具体的值
            item = self._getItemFromName(self.lines_QList, fit_type)
            para_dict = item.data(-1)
        except Exception as e:
            self.informMsg(str(e))
            return
        if re.findall("2s", fit_type) != [] or re.findall("two", fit_type) != []:
            try:
                theta,gap1,gap2 = self.fitting_methods.generate_two_band_gap_func(fit_type=fit_type, kwargs=para_dict)
                self.ax_inset.plot(theta, gap1, color="tab:red", ls="solid", lw=1.5, label="legendre3")
                self.ax_inset.plot(theta, gap2, color="tab:blue", ls="solid", lw=1.5, label="legendre4")
                self.canvas.draw()
            except Exception as e:
                self.informMsg(str(e))
                return
        else:
            try:
                theta,gap = self.fitting_methods.generate_one_band_gap_func(fit_type=fit_type, kwargs=para_dict)
                self.ax_inset.plot(theta, gap, color="tab:blue", ls="solid", lw=1.5, label="legendre3")
                self.ax_inset.tick_params(grid_color="palegoldenrod")
                self.canvas.draw()
            except Exception as e:
                self.informMsg(str(e))
                return

    def _handleOnDeleteDataFromDataList(self):  # 图像连同数据全部销毁 # 如果该函数被调用,则一定是选中对象了,这说明List中至少有一个item
        item = self.data_QList.currentItem()
        filename = item.text()
        if filename in self.raw_bias_dict.keys():
            del self.raw_bias_dict[filename]
            del self.raw_data_dict[filename]

        if filename == self.filename_present:
            self._handleOn_Remove_All()
            self.filename_present = None
            self.fit_lines_dict = None
            self.parameters_dict = None
            self.process_data_present = None
            self.process_bias_present = None
            self.lines_QList.clear()  # 删除所有的item

        row = self.data_QList.row(item)
        self.data_QList.takeItem(row)

    def _handleOn_DeleteAllFromDataList(self):
        reply = self.questionMsg("是否要全部刪除")
        if reply == True:
            row = self.data_QList.count()
            for i in range(row):
                self._handleOnDeleteDataFromDataList()

    def _handleOnDeleteDataFromLineList(self):  # 如果该函数被调用,则一定是选中对象了,这说明List中至少有一个item
        item = self.lines_QList.currentItem()
        row = self.lines_QList.row(item)
        self.lines_QList.takeItem(row)
        try:
            item.data(4).remove()
            self.canvas.draw()
        except Exception as e:
            self.informMsg(str(e) + ":拟合曲线删除失败")

    def _SaveLineListByItem(self,item:QListWidgetItem) -> bool:
        if item is None:
            self.informMsg("未选中atom_list中的item")
            return False
        try:
            fileName = item.text() + ".dat"
            bias_data_pro = np.vstack((self.process_bias_present, item.data(1)))
            bias_data_pro = bias_data_pro.T
            fileName, filetype = QFileDialog.getSaveFileName(self.main_frame, "文件保存", "./" + fileName, "DAT Files (.dat)")
            with open(fileName, 'w') as f:
                try:
                    tplt = '{:<30}{:<20}'
                    str_1 = tplt.format("bias","dIdv") + "\n"
                    f.writelines(str_1)
                    for row in bias_data_pro:
                        str_1 = tplt.format(row[0], row[1]) + "\n"
                        f.writelines(str_1)
                except Exception as e:
                    self.informMsg("保存失败")
                    return False
            self.informMsg("已经保存")
            return True
        except Exception as e:
            self.informMsg(str(e))
            return True

    def _handleOnSaveLineList(self):
        item_now = self.lines_QList.currentItem()
        self._SaveLineListByItem(item_now)

    def _handleOnSave_All_LineListItem(self):
        row = self.lines_QList.count()
        for i in range(row):
            item_now = self.lines_QList.item(i)
            if_true = self._SaveLineListByItem(item_now)
            if if_true == False:
                self.informMsg("第%d文件的保存失败，请重新保存"%i)

    def _handleOn_Remove_All(self):  # self.button_remove_all的调用函数
        try:
            row = self.data_QList.count()
            for i in range(row):
                self.data_QList.item(i).data(4).remove()
        except Exception as e:
            self.informMsg(str(e) + ":實驗數據抹除失败")
        try:
            row = self.lines_QList.count()
            for i in range(row):
                self.lines_QList.item(i).data(4).remove()
        except Exception as e:
            self.informMsg(str(e) + ":擬合數據抹除失败")
        try:
            self.axes.cla()
            self.axes_error.cla()
            self.canvas.draw()
        except:
            self.informMsg("something wrong")
            return

    # 数据预处理:偏压范围重调
    def _handleOn_Symmetry_process(self):
        return

    def _handleOn_Amplification(self):
        item = self._getItemFromName(self.data_QList, self.filename_present)
        bias = item.data(-1)[0]
        self.process_bias_present = self.Amplification(bias)
        bias_data_pro = np.vstack((self.process_bias_present, self.process_data_present))
        item.setData(1, bias_data_pro)

    def Amplification(self,bias):
        try:
            amp_factor = float(self.text_parameter_amplifier.text())
            bias_new = np.multiply(bias, amp_factor)
            return bias_new
        except Exception as e:
            self.informMsg(str(e) + ":放大失败")
            return False

    # 去掉背底,通过线性拟合方法,因为信号通过放大得到,设原始信号为f(v),放大因子为a(v),则实际接收到的信号强度为a(v)f(v)
    # 而由于粒子-空穴对称性,DIDV谱应该关于v=0是对称的,如果GAP的一侧比另一侧要高,那么可能是由于放大因子a(v)随着偏压变化而变化导致,
    # 这时可以通过线性拟合(或更高阶)的方式矫正放大因子不为常数的情况,保证GAP两侧具有对称性
    def _handleOn_Remove_BackGround(self):
        if self.if_rm_bkgd == True:
            self.informMsg("已经消除了背底倾斜")
        else:
            item = self._getItemFromName(self.data_QList, self.filename_present)
            bias = item.data(-1)[0]
            data = item.data(-1)[1]
            self.process_data_present = self.Remove_BackGround(bias,data)
            bias_data_pro = np.vstack((self.process_bias_present, self.process_data_present))
            item.setData(1, bias_data_pro)

    def Remove_BackGround(self,bias,data):
        # 拟合点:GAP两侧趋于平行的部分对应的偏压和DIDV,一般取为总点数的1/10
        try:
            tot_Num = len(bias)
            fit_Num = int(np.floor(tot_Num / 10))
            bias_test = np.hstack((bias[0:fit_Num], bias[tot_Num - fit_Num:tot_Num]))  # x:y并不会取到y
            data_test = np.hstack((data[0:fit_Num], data[tot_Num - fit_Num:tot_Num]))  # x:y并不会取到y
            try:
                def f_lin(x: float, k: float, b: float):  # 定义线性函数,xdata对应的是第一个形参
                    return k * x + b

                k, b = optimize.curve_fit(f_lin, xdata=bias_test, ydata=data_test,
                                          bounds=([-np.inf, -np.inf], [np.inf, np.inf]))[0]
                self.informMsg("背底倾斜的斜率:" + str(k) + ", 截距:" + str(b))
                factor = k * bias + b
                data_new = data / factor
                return data_new
            except:
                self.informMsg("背底去除失败")
                return False
        except Exception as e:
            self.informMsg(str(e) + ":当前没有数据")
            return False

    def _handleOn_Data_Process_All(self):
        reply = self.questionMsg("是否要全部處理？請先檢查一下放大倍數是否寫了？？")
        if reply == False:
            return
        row = self.data_QList.count()
        for i in range(row):
            item = self.data_QList.item(i)
            bias = item.data(-1)[0]
            data = item.data(-1)[1]
            data_new = self.Remove_BackGround(bias, data)
            bias_new = self.Amplification(bias)
            try:  # 如果不加try就會報錯，因爲格式不同
                if data_new == False:
                    self.informMsg("扣除背底時出錯")
                    return
                if bias_new == False:
                    self.informMsg("放大偏壓時出錯")
                    return
            except:
                try:
                    bias_data_new = np.vstack((bias_new, data_new))
                    item.setData(1, bias_data_new)
                    if item.text() == self.filename_present:
                        self.process_bias_present = bias_new
                        self.process_data_present = data_new
                        self.if_rm_bkgd = True
                except Exception as e:
                    self.informMsg(str(e) + ":處理失敗")

    def get_fixed_para(self, para_name):
        para_val = self.lineEdit_widget[para_name].text()
        try:
            para_val = float(para_val)
        except:
            para_val = None
        return para_val

    def _handleOn_Fitting_Fixed_Para(self):
        reply = self.questionMsg("確定擬合？？")
        if reply == False:
            return
        for key in self.lineEdit_widget.keys():
            self.para_fixed[key] = self.get_fixed_para(key)
        try:
            self.fitting_family(self.para_fixed)
        except Exception as e:
            self.informMsg(str(e))

    def fitting_family(self,kwargs:dict):
        fit_type = self.current_fitting_type  # 只有进行拟合之后self.current_fitting_type才会有具体的值
        item = self._getItemFromName(self.lines_QList, fit_type)
        try:
            para = item.data(-1)
            try:
                for key in kwargs.keys():
                    para[key] = kwargs[key]
                item.setData(-1, para)
                para["bias"] = self.process_bias_present
                import re
                if re.findall("2s",fit_type) != [] or re.findall("two",fit_type) != []:
                    fit_data = self.fitting_methods.cal_two_band_with_given_parameters(fit_type=fit_type, kwargs=para)
                else:
                    fit_data = self.fitting_methods.cal_one_band_with_given_parameters(fit_type=fit_type, kwargs=para)
                item.setData(1, fit_data)
                item.data(4).remove()
                try:
                    offset = float(self.lineEdit_widget["offset"].text())
                except:
                    offset = 0
                bias = self.process_bias_present + offset
                self.canvas_plot_from_data(bias, fit_data, item)
            except Exception as e:
                self.informMsg(str(e))
        except:
            self.informMsg("请先拟合一次")
            return

    def _handleOn_Fitting_by_Methods(self):
        parameters = self.getParameters_bound_FromInput()
        parameters["bias"] = self.process_bias_present
        parameters["didv"] = self.process_data_present
        if self.checkBox_one_band.isChecked() == False and self.checkBox_two_band.isChecked() == False:
            self.informMsg("请选择其中一个checkbox")
            return
        elif self.checkBox_one_band.isChecked() == True and self.checkBox_two_band.isChecked() == False:
            fit_type = self.comboBox_one_band.currentText()
            para, para_dict, fit_data = self.fitting_methods.fitting_one_band_methods_manager(fit_type=fit_type, kwargs=parameters)
            self.informMsg(self.fitting_methods.info_one_band[fit_type]%para)
            self.setParameters(para_dict)
        elif self.checkBox_one_band.isChecked() == False and self.checkBox_two_band.isChecked() == True:
            fit_type = self.comboBox_two_band.currentText()
            para, para_dict, fit_data = self.fitting_methods.fitting_two_band_methods_manager(fit_type=fit_type, kwargs=parameters)
            self.informMsg(self.fitting_methods.info_two_band[fit_type] % para)
            self.setParameters(para_dict)
        else:
            self.informMsg("两个checkbox不能同时选择")
            return

        if "item":
            item = self._getItemFromName(self.lines_QList, fit_type)
            item.setData(-1, para_dict)
            item.setData(1, fit_data)
            self.lines_QList.addItem(item)
            self.lines_QList.sortItems()
            self.lines_QList.setCurrentItem(item)
            self.current_fitting_type = fit_type
            try:
                item.data(4).remove()
            except Exception as e:
                self.informMsg(str(e) + ":第一次用该函数拟合")
            try:
                offset = float(self.lineEdit_widget["offset"].text())
            except:
                offset = 0
            bias = self.process_bias_present + offset
            self.canvas_plot_from_data(bias, fit_data, item)

    def getpara_bound(self,para_name:str):
        lower_bound = self.lower_bound[para_name].text()
        upper_bound = self.upper_bound[para_name].text()
        try:
            lower_bound = float(lower_bound)
            if para_name == "x":
                if lower_bound < 0:
                    self.informMsg("x的下限是0")
                    lower_bound = 0
                    self.lower_bound[para_name].setText("0")
        except:
            if para_name == "offset":
                lower_bound = -0.1
            elif para_name == "phi":
                lower_bound = 0
            elif para_name == "x":
                lower_bound = 0
            else:
                lower_bound = 0.12
        try:
            upper_bound = float(upper_bound)
            if para_name == "x":
                if upper_bound > 1:
                    self.informMsg("x的上限是0")
                    upper_bound = 0
                    self.upper_bound[para_name].setText("1")
        except:
            if para_name == "offset":
                upper_bound = 0.1
            elif para_name == "phi":
                upper_bound = 2
            elif para_name == "x":
                upper_bound = 1
            else:
                upper_bound = 0.16
        self.para_dict[para_name] = [lower_bound, upper_bound]

    def getParameters_bound_FromInput(self):
        if self.check() == False:
            return
        for key in self.lower_bound.keys():
            self.getpara_bound(key)
        return self.para_dict

    def setParameters(self, kwargs:dict):
        for key in kwargs.keys():
            self.lineEdit_widget[key].setText(str(kwargs[key]))

    def check(self) -> bool:
        try:
            if self.process_bias_present == None:
                self.informMsg("目前无数据")
                return False
            if self.process_data_present == None:
                self.informMsg("目前无数据")
                return False
        except Exception as e:
            # self.informMsg(str(e))  # the truth value of an array with more than one element is ambiguous.
            return True

    # 两个信息提示框
    def informMsg(self, msg: str):
        msgBox = QMessageBox()
        msgBox.setWindowTitle("inform")
        msgBox.setText(msg)
        msgBox.exec_()

    def questionMsg(cls, msg: str):
        msgBox = QMessageBox()
        msgBox.setWindowTitle("确认框")
        reply = QMessageBox.information(msgBox,"标题",msg,QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            return True
        if reply == QMessageBox.No:
            return False
        msgBox.exec_()

# to do list:如何鼠标触碰后显示QWidgetList中的item的名称？？
if __name__ == "__main__":
    app = QApplication(sys.argv)
    form = PlotFrame_SC_GAP()
    form.mainWindow.show()
    app.exec_()