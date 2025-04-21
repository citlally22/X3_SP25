from PyQt5 import QtWidgets as qtw


class DieselView:
    def __init__(self):
        self.widgets = {}
        self.ax = None
        self.canvas = None
        self.results = None  # Store results for redrawing

    def setWidgets(self, widget_dict):
        self.widgets = widget_dict
        self.ax = widget_dict['ax']
        self.canvas = widget_dict['canvas']

    def updateResults(self, results):
        self.results = results  # Store results
        w = self.widgets
        is_metric = w['rdo_Metric'].isChecked()
        # Convert temperatures from Kelvin to Rankine if in English units
        T1 = results['T1'] * 9 / 5 if not is_metric else results['T1']
        T2 = results['T2'] * 9 / 5 if not is_metric else results['T2']
        T3 = results['T3'] * 9 / 5 if not is_metric else results['T3']
        T4 = results['T4'] * 9 / 5 if not is_metric else results['T4']
        # Work and heat are in kJ/kg (SI units); convert to BTU/lbm if English
        factor = 1.0 if is_metric else 0.4299  # kJ/kg to BTU/lbm
        w_exp = results['w_exp'] * factor
        w_comp = results['w_comp'] * factor
        q_in = results['q_in'] * factor

        w['le_T1'].setText(f"{T1:.2f}")
        w['le_T2'].setText(f"{T2:.2f}")
        w['le_T3'].setText(f"{T3:.2f}")
        w['le_T4'].setText(f"{T4:.2f}")
        w['le_HeatAdded'].setText(f"{q_in:.2f}")
        w['le_CompressionStroke'].setText(f"{w_comp:.2f}")
        w['le_PowerStroke'].setText(f"{w_exp:.2f}")
        w['le_Efficiency'].setText(f"{results['eff']:.2f}")

    def updateUnits(self, is_metric):
        w = self.widgets
        if is_metric:
            # Metric units: Kelvin, bar, m³
            w['lbl_T1Units'].setText("K")
            w['lbl_T2Units'].setText("K")
            w['lbl_T3Units'].setText("K")
            w['lbl_T4Units'].setText("K")
            w['lbl_PowerStrokeUnits'].setText("kJ/kg")
            w['lbl_CompressionStrokeUnits'].setText("kJ/kg")
            w['lbl_HeatInUnits'].setText("kJ/kg")
        else:
            # English units: Rankine, atm, ft³
            w['lbl_T1Units'].setText("R")
            w['lbl_T2Units'].setText("R")
            w['lbl_T3Units'].setText("R")
            w['lbl_T4Units'].setText("R")
            w['lbl_PowerStrokeUnits'].setText("BTU/lbm")
            w['lbl_CompressionStrokeUnits'].setText("BTU/lbm")
            w['lbl_HeatInUnits'].setText("BTU/lbm")
        # Redraw results and plot with new units
        if self.results:
            self.updateResults(self.results)
            self.updatePlot(
                self.results,
                w['cmb_Abcissa'].currentText(),
                w['cmb_Ordinate'].currentText(),
                w['chk_LogAbcissa'].isChecked(),
                w['chk_LogOrdinate'].isChecked(),
                self.ax,
                self.canvas
            )

    def updatePlot(self, results, x_type, y_type, log_x, log_y, ax, canvas):
        self.ax = ax
        self.canvas = canvas
        self.plotPV(results['pv'], log_x, log_y)

    def plotPV(self, pv_data, log_x=False, log_y=False):
        self.ax.clear()
        v, p = zip(*pv_data)
        if self.widgets['rdo_Metric'].isChecked():
            # Metric: m³, Pa
            v_label = 'Volume (m³)'
            p_label = 'Pressure (Pa)'
        else:
            # English: ft³, atm
            v = [vol / 0.0283168 for vol in v]  # m³ to ft³
            p = [press / 101325 for press in p]  # Pa to atm
            v_label = 'Volume (ft³)'
            p_label = 'Pressure (atm)'
        self.ax.plot(v, p, 'r-', label='Diesel Cycle')
        self.ax.set_xlabel(v_label)
        self.ax.set_ylabel(p_label)
        self.ax.set_title('Diesel PV Diagram')
        if log_x:
            self.ax.set_xscale('log')
        if log_y:
            self.ax.set_yscale('log')
        self.ax.grid(True)
        self.ax.legend()
        self.canvas.draw()