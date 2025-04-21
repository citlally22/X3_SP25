from DieselModel import DieselModel
from DieselView import DieselView

class DieselController:
    """Controller for Diesel cycle calculations and GUI updates."""

    def __init__(self):
        """Initialize Diesel model and view."""
        self.model = DieselModel()
        self.view = DieselView()

    def setWidgets(self, widget_dict):
        """Set GUI widgets for the view."""
        self.view.setWidgets(widget_dict)

    def calculate(self, inputs):
        try:
            # Validate inputs
            for key in ['CR', 'rc', 'T_high', 'T_low', 'P0', 'V0']:
                if not isinstance(inputs[key], str) or not inputs[key].strip() or not float(inputs[key]) > 0:
                    raise ValueError(f"Invalid {key}: must be a positive number")

            # Convert units based on rdo_Metric
            is_metric = self.view.widgets['rdo_Metric'].isChecked()
            T_low = float(inputs['T_low'])
            T_high = float(inputs['T_high'])
            V0 = float(inputs['V0'])
            P0 = float(inputs['P0'])

            # Model expects SI units (Kelvin, Pa, m³)
            if not is_metric:  # English units (Rankine, atm, ft³)
                T_low = T_low * 5 / 9  # Rankine to Kelvin
                T_high = T_high * 5 / 9  # Rankine to Kelvin
                V0 = V0 * 0.0283168  # ft³ to m³
                P0 = P0 * 101325  # atm to Pa
            else:  # Metric units (Kelvin, bar, m³)
                P0 = P0 * 1e5  # bar to Pa

            self.model.r = float(inputs['CR'])
            self.model.rc = float(inputs['rc'])
            self.model.T1 = T_low
            self.model.T_high = T_high  # Pass T_high to the model
            self.model.P1 = P0
            self.model.V1 = V0
            self.model.calculate()
            self.view.updateResults(self.model.results)
            self.view.plotPV(self.model.results['pv'])
        except ValueError as e:
            raise ValueError(f"Input error: {str(e)}")
        except Exception as e:
            raise Exception(f"Calculation error: {str(e)}")

    def updateView(self):
        """Update GUI results, units, and plot."""
        widgets = self.view.widgets
        self.view.updateResults(self.model.results)
        self.view.updateUnits(widgets['rdo_Metric'].isChecked())
        self.view.updatePlot(
            self.model.results,
            widgets['cmb_Abcissa'].currentText(),
            widgets['cmb_Ordinate'].currentText(),
            widgets['chk_LogAbcissa'].isChecked(),
            widgets['chk_LogOrdinate'].isChecked(),
            widgets['ax'],
            widgets['canvas']
        )