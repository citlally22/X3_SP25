# DieselView.py
from PyQt5 import QtWidgets as qtw

class DieselView:
    def __init__(self):
        self.widgets = {}
        self.ax = None
        self.canvas = None

    def setWidgets(self, widget_dict):
        self.widgets = widget_dict
        self.ax = widget_dict['ax']
        self.canvas = widget_dict['canvas']

    def updateResults(self, results):
        w = self.widgets
        w['le_T1'].setText(f"{results['T1']:.2f}")
        w['le_T2'].setText(f"{results['T2']:.2f}")
        w['le_T3'].setText(f"{results['T3']:.2f}")
        w['le_T4'].setText(f"{results['T4']:.2f}")
        w['le_HeatAdded'].setText(f"{results['q_in']:.2f}")
        w['le_CompressionStroke'].setText(f"{results['q_out']:.2f}")
        w['le_PowerStroke'].setText(f"{results['w_net']:.2f}")
        w['le_Efficiency'].setText(f"{results['eff']:.2f}")

    def plotPV(self, pv_data):
        self.ax.clear()
        v, p = zip(*pv_data)
        self.ax.plot(v, p, 'r-', label='Diesel Cycle')
        self.ax.set_xlabel('Volume (m³)')
        self.ax.set_ylabel('Pressure (Pa)')
        self.ax.set_title('Diesel PV Diagram')
        self.ax.grid(True)
        self.ax.legend()
        self.canvas.draw()