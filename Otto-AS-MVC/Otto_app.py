import sys
from PyQt5 import QtWidgets as qtw
from Otto import ottoCycleController
from DieselController import DieselController
from Otto_GUI import Ui_Form
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg


class MainWindow(qtw.QWidget, Ui_Form):
    """Main GUI window for the Thermodynamic Cycle Calculator.

    Integrates Otto and Diesel cycle controllers to perform calculations and display
    results and plots. Manages user inputs, unit systems, and plot updates through a
    PyQt5 interface.

    Attributes:
        calculated (bool): Tracks whether a calculation has been successfully performed.
        figure (matplotlib.figure.Figure): Matplotlib figure for plotting.
        canvas (FigureCanvasQTAgg): Matplotlib canvas for rendering plots.
        ax (matplotlib.axes.Axes): Matplotlib axes for plotting.
        ottoController (ottoCycleController): Controller for Otto cycle calculations.
        dieselController (DieselController): Controller for Diesel cycle calculations.
        shared_widgets (dict): Dictionary mapping widget names to Qt objects.
        current_cycle (str): Current cycle type ('Otto' or 'Diesel').
    """

    def __init__(self):
        """Initializes the main window, sets up the GUI, and configures controllers."""
        super().__init__()
        self.setupUi(self)  # Set up UI from Otto_GUI
        self.setWindowTitle("Thermodynamic Cycle Calculator")
        self.calculated = False  # Track calculation status

        # Initialize Matplotlib canvas for plotting
        self.figure = Figure(figsize=(8, 8), tight_layout=True, facecolor='none')
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot()
        self.main_VerticalLayout.addWidget(self.canvas)  # Add canvas to layout

        # Instantiate controllers
        self.ottoController = ottoCycleController(ax=self.ax)
        self.dieselController = DieselController()

        # Define shared widgets for both controllers
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

        # Assign widgets to controllers
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

        # Set default cycle
        self.current_cycle = 'Otto'
        self.cmb_Cycle.setCurrentText("Otto")

        self.show()  # Display the window

    def setCycle(self):
        """Updates the current cycle type based on user selection.

        Triggers plot update if a calculation has been performed.
        """
        self.current_cycle = self.cmb_Cycle.currentText()  # Update cycle type
        if self.calculated:
            self.updatePlot()  # Refresh plot for new cycle

    def updatePlot(self):
        """Updates the plot based on the current cycle and plot settings.

        Rebuilds plotting data and updates the view for the selected cycle.
        """
        if not self.calculated:
            return  # Skip if no calculation has been performed

        self.ax.clear()  # Clear previous plot
        try:
            if self.current_cycle == 'Otto':
                self.ottoController.buildDataForPlotting()  # Generate plotting data
                self.ottoController.updateView()  # Update Otto plot
            else:
                self.dieselController.updateView()  # Update Diesel plot (no buildDataForPlotting)
        except Exception as e:
            qtw.QMessageBox.critical(self, "Plot Error", f"Failed to update plot: {str(e)}")
        finally:
            self.canvas.draw()  # Redraw canvas

    def calculateCycle(self):
        """Performs calculations for the selected cycle based on user inputs.

        Validates inputs, triggers the appropriate controller's calculation, and updates the view.
        """
        # Collect inputs from GUI
        inputs = {
            'T_high': self.le_THigh.text(),
            'T_low': self.le_TLow.text(),
            'P0': self.le_P0.text(),
            'V0': self.le_V0.text(),
            'CR': self.le_CR.text(),
            'rc': "2"  # Default cutoff ratio for Diesel
        }
        try:
            if self.current_cycle == 'Otto':
                self.ottoController.calc()  # Run Otto calculation
                self.ottoController.buildDataForPlotting()  # Prepare plotting data
                self.ottoController.updateView()  # Update view
            elif self.current_cycle == 'Diesel':
                self.dieselController.calculate(inputs)  # Run Diesel calculation
                self.dieselController.updateView()  # Update view
            self.calculated = True  # Mark calculation as successful
            self.updatePlot()  # Refresh plot
        except Exception as e:
            # Display detailed error message with inputs
            error_msg = (
                f"Calculation failed for {self.current_cycle} cycle.\n"
                f"Inputs: T_high={inputs['T_high']}, T_low={inputs['T_low']}, "
                f"P0={inputs['P0']}, V0={inputs['V0']}, CR={inputs['CR']}, rc={inputs['rc']}\n"
                f"Error: {str(e)}"
            )
            qtw.QMessageBox.critical(self, "Calculation Error", error_msg)
            self.calculated = False  # Reset calculation status
            return

    def setUnits(self):
        """Updates the unit system for both controllers and refreshes the view.

        Sets SI or English units based on the metric radio button and updates the GUI.
        """
        is_metric = self.rdo_Metric.isChecked()
        try:
            # Update unit systems for both controllers
            self.ottoController.model.units.set(SI=is_metric)
            # DieselController's model may not have a units object, so skip
            if hasattr(self.dieselController.model, 'units'):
                self.dieselController.model.units.set(SI=is_metric)

            # Update view units for the current cycle
            if self.current_cycle == 'Otto':
                self.ottoController.view.updateDisplayWidgets(Model=self.ottoController.model)
            else:
                self.dieselController.view.updateUnits(is_metric)
        except Exception as e:
            qtw.QMessageBox.critical(self, "Units Error", f"Failed to update units: {str(e)}")

        # Refresh plot if a calculation has been performed
        if self.calculated:
            self.updatePlot()


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)  # Initialize Qt application
    window = MainWindow()  # Create main window
    sys.exit(app.exec_())  # Run Qt event loop