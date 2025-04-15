# EX3P1SP22.py
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
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.resize(600, 900)  # Increase height for toolbar visibility
        self.setupImageLabel()
        self.setupInputFields()
        self.setupPlot()
        self.show()

    def setupImageLabel(self):
        self.pixMap = qtg.QPixmap()
        self.pixMap.load("Circuit1.png")
        self.image_label = qtw.QLabel()
        self.image_label.setPixmap(self.pixMap)
        self.layout_GridInput.addWidget(self.image_label, 0, 0, 1, 2)

    def setupInputFields(self):
        self.L_label = qtw.QLabel("L (H):")
        self.L_input = QLineEdit("20")
        self.R_label = qtw.QLabel("R (Ω):")
        self.R_input = QLineEdit("10")
        self.C_label = qtw.QLabel("C (F):")
        self.C_input = QLineEdit("0.05")
        self.Vm_label = qtw.QLabel("Voltage Magnitude (V):")
        self.Vm_input = QLineEdit("20")
        self.freq_label = qtw.QLabel("Frequency (rad/s):")
        self.freq_input = QLineEdit("20")
        self.phase_label = qtw.QLabel("Phase (rad):")
        self.phase_input = QLineEdit("0")

        self.simulate_button = QPushButton("Simulate")
        self.simulate_button.clicked.connect(self.simulate)

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
        self.layout_GridInput.addWidget(self.simulate_button, 7, 0, 1, 2)

    def setupPlot(self):
        self.figure, self.ax = plt.subplots(figsize=(5, 4))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        self.verticalLayout.addLayout(plot_layout)

    def simulate(self):
        try:
            # Get user inputs
            L = float(self.L_input.text())
            R = float(self.R_input.text())
            C = float(self.C_input.text())
            Vm = float(self.Vm_input.text())
            freq = float(self.freq_input.text())
            phase = float(self.phase_input.text())

            # Validate inputs
            if L <= 0 or R <= 0 or C <= 0:
                qtw.QMessageBox.warning(self, "Invalid Input", "L, R, and C must be positive.")
                return
            if freq < 0:
                qtw.QMessageBox.warning(self, "Invalid Input", "Frequency must be non-negative.")
                return

            # Define the system of ODEs
            def circuit_model(y, t, L, R, C, Vm, freq, phase):
                i1, vC = y
                v_t = Vm * np.sin(freq * t + phase)
                di1_dt = (v_t - vC) / L
                dvC_dt = (i1 - vC / R) / C
                return [di1_dt, dvC_dt]

            # Time points
            t = np.linspace(0, 2, 2000)  # Extended to 2 seconds
            y0 = [0, 0]  # Initial conditions: i1(0) = 0, vC(0) = 0
            sol = odeint(circuit_model, y0, t, args=(L, R, C, Vm, freq, phase))
            i1 = sol[:, 0]
            vC = sol[:, 1]
            i2 = vC / R

            # Plot results
            self.ax.clear()
            self.ax.plot(t, i1, label="i1 (A)", linewidth=2, color="blue")
            self.ax.plot(t, i2, label="i2 (A)", linewidth=2, color="orange")
            self.ax.plot(t, vC, label="vC (V)", linewidth=2, color="green")
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Current (A) / Voltage (V)")
            self.ax.set_title("RLC Circuit Transient Response")
            self.ax.legend()
            self.ax.grid(True)
            self.canvas.draw()
        except ValueError:
            qtw.QMessageBox.warning(self, "Invalid Input", "Please enter valid numeric values.")
        except Exception as e:
            qtw.QMessageBox.critical(self, "Error", f"Simulation failed: {e}")

if __name__ == "__main__":
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = main_window()
    sys.exit(app.exec_())