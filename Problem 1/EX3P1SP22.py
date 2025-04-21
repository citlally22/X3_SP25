from PyQt5.QtWidgets import QApplication, QWidget, QLineEdit, QPushButton, QVBoxLayout, QScrollArea, QMessageBox
import PyQt5.QtWidgets as qtw
import PyQt5.QtGui as qtg
from Problem1 import Ui_Form
import sys
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar


class main_window(Ui_Form, QWidget):
    """Main application window for RLC circuit simulation with PyQt5 GUI.

    Inherits from Ui_Form (generated UI) and QWidget to create a window
    that allows users to input circuit parameters, simulate an RLC circuit,
    and visualize the transient response.
    """

    def __init__(self):
        """Initialize the main window, set up UI, and configure components."""
        super().__init__()
        self.setupUi(self)  # Load UI from Problem1.Ui_Form
        self.resize(600, 900)  # Set window size (width: 600px, height: 900px)
        self.setupImageLabel()  # Configure circuit image display
        self.setupInputFields()  # Set up input fields and simulate button
        self.setupPlot()  # Initialize matplotlib plot canvas
        self.show()  # Display the window

    def setupImageLabel(self):
        """Set up the QLabel to display the RLC circuit diagram."""
        self.pixMap = qtg.QPixmap()  # Create QPixmap for image
        self.pixMap.load(" Circuit1.png")  # Load circuit diagram image
        self.image_label = qtw.QLabel()  # Create QLabel for image display
        self.image_label.setPixmap(self.pixMap)  # Set image to label
        # Add image label to grid layout (row 0, spans 2 columns)
        self.layout_GridInput.addWidget(self.image_label, 0, 0, 1, 2)

    def setupInputFields(self):
        """Configure input fields and simulate button for user inputs."""
        # Create labels and input fields for circuit parameters
        self.L_label = qtw.QLabel("L (H):")
        self.L_input = QLineEdit("20")  # Default inductance: 20 H
        self.R_label = qtw.QLabel("R (Ω):")
        self.R_input = QLineEdit("10")  # Default resistance: 10 Ω
        self.C_label = qtw.QLabel("C (F):")
        self.C_input = QLineEdit("0.05")  # Default capacitance: 0.05 F
        self.Vm_label = qtw.QLabel("Voltage Magnitude (V):")
        self.Vm_input = QLineEdit("20")  # Default voltage magnitude: 20 V
        self.freq_label = qtw.QLabel("Frequency (rad/s):")
        self.freq_input = QLineEdit("20")  # Default frequency: 20 rad/s
        self.phase_label = qtw.QLabel("Phase (rad):")
        self.phase_input = QLineEdit("0")  # Default phase: 0 rad

        # Create and configure simulate button
        self.simulate_button = QPushButton("Simulate")
        self.simulate_button.clicked.connect(self.simulate)  # Connect to simulate method

        # Add widgets to grid layout (label in column 0, input in column 1)
        self.layout_GridInput.addWidget(self.L_label, 1, 0)
        self.layout_GridInput.addWidget(self.L_input, 1, 1)
        self.layout_GridInput.addWidget(self.R_label, 2, 0)
        self.layout_GridInput.addWidget(self.R_input, 2, 1)
        self.layout_GridInput.addWidget(self.C_label, 3, 0)
        self.layout_GridInput.addWidget(self.C_input, 3, 1)
        self.layout_GridInput.addWidget(self.Vm_label, 4, 0)
        self.layout_GridInput.addWidget(self.Vm_input, 4, 1)
        self.layout_GridInput.addWidget(self.freq_label, 5, 0)
        self.layout_GridInput.addWidget(self.freq_input, 5, 1)
        self.layout_GridInput.addWidget(self.phase_label, 6, 0)
        self.layout_GridInput.addWidget(self.phase_input, 6, 1)
        # Add simulate button (row 7, spans 2 columns)
        self.layout_GridInput.addWidget(self.simulate_button, 7, 0, 1, 2)

    def setupPlot(self):
        """Initialize matplotlib figure, canvas, and toolbar for plotting."""
        # Create figure and axis for plotting (5x4 inches)
        self.figure, self.ax = plt.subplots(figsize=(5, 4))
        self.canvas = FigureCanvas(self.figure)  # Create canvas for figure
        self.toolbar = NavigationToolbar(self.canvas, self)  # Add navigation toolbar
        plot_layout = QVBoxLayout()  # Create vertical layout for plot
        plot_layout.addWidget(self.toolbar)  # Add toolbar to layout
        plot_layout.addWidget(self.canvas)  # Add canvas to layout
        self.verticalLayout.addLayout(plot_layout)  # Add plot layout to main UI

    def simulate(self):
        """Simulate RLC circuit and plot transient response.

        Reads user inputs, validates them, solves the circuit's differential
        equations using odeint, and plots the results (i1, i2, vC) on the canvas.
        """
        try:
            # Parse user inputs as floats
            L = float(self.L_input.text())  # Inductance (H)
            R = float(self.R_input.text())  # Resistance (Ω)
            C = float(self.C_input.text())  # Capacitance (F)
            Vm = float(self.Vm_input.text())  # Voltage magnitude (V)
            freq = float(self.freq_input.text())  # Frequency (rad/s)
            phase = float(self.phase_input.text())  # Phase (rad)

            # Validate inputs
            if L <= 0 or R <= 0 or C <= 0:
                qtw.QMessageBox.warning(
                    self, "Invalid Input", "L, R, and C must be positive."
                )
                return
            if freq < 0:
                qtw.QMessageBox.warning(
                    self, "Invalid Input", "Frequency must be non-negative."
                )
                return

            # Define ODE system for RLC circuit
            def circuit_model(y, t, L, R, C, Vm, freq, phase):
                """Model RLC circuit dynamics as a system of ODEs.

                Args:
                    y (list): State variables [i1, vC]
                    t (float): Time (s)
                    L (float): Inductance (H)
                    R (float): Resistance (Ω)
                    C (float): Capacitance (F)
                    Vm (float): Voltage magnitude (V)
                    freq (float): Frequency (rad/s)
                    phase (float): Phase (rad)

                Returns:
                    list: Derivatives [di1_dt, dvC_dt]
                """
                i1, vC = y  # Current through inductor, voltage across capacitor
                v_t = Vm * np.sin(freq * t + phase)  # Source voltage
                di1_dt = (v_t - vC) / L  # Rate of change of i1
                dvC_dt = (i1 - vC / R) / C  # Rate of change of vC
                return [di1_dt, dvC_dt]

            # Set up simulation parameters
            t = np.linspace(0, 2, 2000)  # Time array: 0 to 2s, 2000 points
            y0 = [0, 0]  # Initial conditions: i1(0) = 0, vC(0) = 0
            # Solve ODEs
            sol = odeint(circuit_model, y0, t, args=(L, R, C, Vm, freq, phase))
            i1 = sol[:, 0]  # Extract i1 (current through inductor)
            vC = sol[:, 1]  # Extract vC (voltage across capacitor)
            i2 = vC / R  # Calculate i2 (current through resistor)

            # Clear previous plot and create new one
            self.ax.clear()
            self.ax.plot(t, i1, label="i1 (A)", linewidth=2, color="blue")
            self.ax.plot(t, i2, label="i2 (A)", linewidth=2, color="orange")
            self.ax.plot(t, vC, label="vC (V)", linewidth=2, color="green")
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Current (A) / Voltage (V)")
            self.ax.set_title("RLC Circuit Transient Response")
            self.ax.legend()  # Show legend
            self.ax.grid(True)  # Enable grid
            self.canvas.draw()  # Update canvas with new plot

        except ValueError:
            # Handle invalid numeric inputs
            qtw.QMessageBox.warning(
                self, "Invalid Input", "Please enter valid numeric values."
            )
        except Exception as e:
            # Handle other simulation errors
            qtw.QMessageBox.critical(self, "Error", f"Simulation failed: {e}")


if __name__ == "__main__":
    """Entry point for the application."""
    app = QApplication.instance()  # Check for existing QApplication
    if not app:
        app = QApplication(sys.argv)  # Create new QApplication if none exists
    app.aboutToQuit.connect(app.deleteLater)  # Clean up on quit
    main_win = main_window()  # Create main window
    sys.exit(app.exec_())  # Start event loop and exit with app status