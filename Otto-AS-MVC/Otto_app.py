from PyQt5 import QtWidgets as qtw
from PyQt5.QtWidgets import QApplication
import sys
from Otto_GUI import Ui_Form
from Otto import ottoCycleController
from DieselController import DieselController

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.setWindowTitle("Thermodynamic Cycle Calculator")
        self.calculated = False

        # Add Matplotlib canvas to GUI
        self.figure = Figure(figsize=(8, 8), tight_layout=True, facecolor='none')
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot()
        self.main_VerticalLayout.addWidget(self.canvas)

        # Instantiate both controllers
        self.ottoController = ottoCycleController()
        self.dieselController = DieselController()

        # Hook up combo box for switching cycles
        self.cmb_Cycle.currentIndexChanged.connect(self.setCycle)

        # Hook up the Calculate button
        self.btn_Calculate.clicked.connect(self.calculateCycle)

        # Hook up log/axis plots
        self.cmb_Abcissa.currentIndexChanged.connect(self.updatePlot)
        self.cmb_Ordinate.currentIndexChanged.connect(self.updatePlot)
        self.chk_LogAbcissa.stateChanged.connect(self.updatePlot)
        self.chk_LogOrdinate.stateChanged.connect(self.updatePlot)

        # Connect Metric toggle
        self.rdo_Metric.toggled.connect(self.setUnits)

        # Widgets used by both controllers
        self.shared_widgets = {
            'lbl_THigh': self.lbl_THigh, 'lbl_TLow': self.lbl_TLow,
            'lbl_P0': self.lbl_P0, 'lbl_V0': self.lbl_V0, 'lbl_CR': self.lbl_CR,
            'le_THigh': self.le_THigh, 'le_TLow': self.le_TLow, 'le_P0': self.le_P0,
            'le_V0': self.le_V0, 'le_CR': self.le_CR,
            'le_T1': self.le_T1, 'le_T2': self.le_T2, 'le_T3': self.le_T3, 'le_T4': self.le_T4,
            'lbl_T1Units': self.lbl_T1Units, 'lbl_T2Units': self.lbl_T2Units,
            'lbl_T3Units': self.lbl_T3Units, 'lbl_T4Units': self.lbl_T4Units,
            'le_PowerStroke': self.le_PowerStroke, 'le_CompressionStroke': self.le_CompressionStroke,
            'le_HeatAdded': self.le_HeatAdded, 'le_Efficiency': self.le_Efficiency,
            'lbl_PowerStrokeUnits': self.lbl_PowerStrokeUnits,
            'lbl_CompressionStrokeUnits': self.lbl_CompressionStrokeUnits,
            'lbl_HeatInUnits': self.lbl_HeatInUnits,
            'rdo_Metric': self.rdo_Metric, 'cmb_Abcissa': self.cmb_Abcissa,
            'cmb_Ordinate': self.cmb_Ordinate,
            'chk_LogAbcissa': self.chk_LogAbcissa, 'chk_LogOrdinate': self.chk_LogOrdinate,
            'ax': self.ax, 'canvas': self.canvas
        }

        self.ottoController.setWidgets(self.shared_widgets)
        self.dieselController.setWidgets(self.shared_widgets)

        # Default cycle
        self.current_cycle = 'Otto'
        self.cmb_Cycle.setCurrentText("Otto")

        self.show()

    def clamp(self, val, low, high):
        try:
            f = float(val)
            return max(min(f, high), low)
        except ValueError:
            return float(low)

    def isfloat(self, value):
        try:
            float(value)
            return True
        except:
            return False

    def setUnits(self):
        if self.current_cycle == 'Otto':
            self.ottoController.updateView()

    def setCycle(self):
        self.current_cycle = self.cmb_Cycle.currentText()

    def updatePlot(self):
        if self.current_cycle == 'Otto':
            self.ottoController.updateView()
        else:
            # Diesel doesn't use cmb_Abcissa/cmb_Ordinate currently
            pass

    def calculateCycle(self):
        inputs = {
            'T_high': self.le_THigh.text(),
            'T_low': self.le_TLow.text(),
            'P0': self.le_P0.text(),
            'V0': self.le_V0.text(),
            'CR': self.le_CR.text(),
            'rc': "2"  # default for Diesel; optionally add a cutoff field
        }

        if self.current_cycle == 'Otto':
            self.ottoController.calc()
        elif self.current_cycle == 'Diesel':
            self.dieselController.calculate(inputs)

# Run the app
if __name__ == '__main__':
    app = QApplication(sys.argv)
    mw = MainWindow()
    sys.exit(app.exec())
