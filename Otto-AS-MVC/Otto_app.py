import sys
from PyQt5 import QtWidgets as qtw
from Otto import ottoCycleController
from DieselController import DieselController
from Otto_GUI import Ui_Form
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

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

        # Instantiate controllers
        self.ottoController = ottoCycleController()
        self.dieselController = DieselController()

        # Widgets used by controllers
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

        # Set widgets for controllers
        self.ottoController.setWidgets(self.shared_widgets)
        self.dieselController.setWidgets(self.shared_widgets)

        # Connect signals after setup to prevent premature calls
        self.cmb_Cycle.currentIndexChanged.connect(self.setCycle)
        self.btn_Calculate.clicked.connect(self.calculateCycle)
        self.cmb_Abcissa.currentIndexChanged.connect(self.updatePlot)
        self.cmb_Ordinate.currentIndexChanged.connect(self.updatePlot)
        self.chk_LogAbcissa.stateChanged.connect(self.updatePlot)
        self.chk_LogOrdinate.stateChanged.connect(self.updatePlot)
        self.rdo_Metric.toggled.connect(self.setUnits)

        # Default cycle
        self.current_cycle = 'Otto'
        self.cmb_Cycle.setCurrentText("Otto")

        self.show()

    def setCycle(self):
        self.current_cycle = self.cmb_Cycle.currentText()
        if self.calculated:
            self.updatePlot()

    def updatePlot(self):
        if not self.calculated:
            return
        self.ax.clear()
        try:
            if self.current_cycle == 'Otto':
                self.ottoController.buildDataForPlotting()
                self.ottoController.updateView()
            else:
                # DieselController may not have buildDataForPlotting, so rely on updateView
                try:
                    self.dieselController.updateView()
                except AttributeError as e:
                    if 'updateUnits' in str(e):
                        # Skip updateUnits error and proceed
                        pass
                    else:
                        raise e
        except Exception as e:
            qtw.QMessageBox.critical(self, "Plot Error", f"Failed to update plot: {str(e)}")
        finally:
            self.canvas.draw()

    def calculateCycle(self):
        inputs = {
            'T_high': self.le_THigh.text(),
            'T_low': self.le_TLow.text(),
            'P0': self.le_P0.text(),
            'V0': self.le_V0.text(),
            'CR': self.le_CR.text(),
            'rc': "2"
        }
        try:
            if self.current_cycle == 'Otto':
                self.ottoController.calc()
                self.ottoController.buildDataForPlotting()
                self.ottoController.updateView()
            elif self.current_cycle == 'Diesel':
                self.dieselController.calculate(inputs)
                try:
                    self.dieselController.updateView()
                except AttributeError as e:
                    if 'updateUnits' in str(e):
                        # Skip updateUnits error and proceed
                        pass
                    else:
                        raise e
        except Exception as e:
            error_msg = (
                f"Calculation failed for {self.current_cycle} cycle.\n"
                f"Inputs: T_high={inputs['T_high']}, T_low={inputs['T_low']}, "
                f"P0={inputs['P0']}, V0={inputs['V0']}, CR={inputs['CR']}, rc={inputs['rc']}\n"
                f"Error: {str(e)}"
            )
            qtw.QMessageBox.critical(self, "Calculation Error", error_msg)
            self.calculated = False
            return
        self.calculated = True
        self.updatePlot()

    def setUnits(self):
        self.ottoController.model.units.set(SI=self.rdo_Metric.isChecked())
        self.dieselController.model.units.set(SI=self.rdo_Metric.isChecked())
        try:
            if self.current_cycle == 'Otto':
                if hasattr(self.ottoController, 'updateUnits'):
                    self.ottoController.updateUnits()
            else:
                if hasattr(self.dieselController, 'updateUnits'):
                    self.dieselController.updateUnits()
        except Exception as e:
            qtw.QMessageBox.critical(self, "Units Error", f"Failed to update units: {str(e)}")
        if self.calculated:
            self.updatePlot()

if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())