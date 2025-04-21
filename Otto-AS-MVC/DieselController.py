from DieselModel import DieselModel
from DieselView import DieselView


class DieselController:
    """Controller for managing Diesel cycle calculations and GUI updates.

    Coordinates between the DieselModel for thermodynamic calculations and
    DieselView for GUI interactions, handling input validation, unit conversions,
    and result display.

    Attributes:
        model (DieselModel): Instance of the Diesel cycle model.
        view (DieselView): Instance of the view for GUI interactions.
    """

    def __init__(self):
        """Initializes the controller with a Diesel model and view."""
        self.model = DieselModel()  # Create model instance
        self.view = DieselView()   # Create view instance

    def setWidgets(self, widget_dict):
        """Sets the GUI widgets for the view.

        Args:
            widget_dict (dict): Dictionary mapping widget names to their objects.
        """
        self.view.setWidgets(widget_dict)  # Pass widget dictionary to view

    def calculate(self, inputs):
        """Performs Diesel cycle calculations based on user inputs.

        Validates inputs, converts units to SI (Kelvin, Pa, m³) as needed, triggers
        model calculations, and updates the view with results.

        Args:
            inputs (dict): Dictionary with keys 'CR', 'rc', 'T_high', 'T_low', 'P0', 'V0'
                          containing string representations of input values.

        Raises:
            ValueError: If any input is invalid (non-numeric, empty, or non-positive).
            Exception: If model calculations fail.
        """
        try:
            # Validate inputs: ensure they are non-empty strings and positive numbers
            for key in ['CR', 'rc', 'T_high', 'T_low', 'P0', 'V0']:
                if not isinstance(inputs[key], str) or not inputs[key].strip() or not float(inputs[key]) > 0:
                    raise ValueError(f"Invalid {key}: must be a positive number")

            # Determine unit system from radio button
            is_metric = self.view.widgets['rdo_Metric'].isChecked()
            T_low = float(inputs['T_low'])
            T_high = float(inputs['T_high'])
            V0 = float(inputs['V0'])
            P0 = float(inputs['P0'])

            # Convert to SI units (Kelvin, Pa, m³) based on unit system
            if not is_metric:  # English units: Rankine, atm, ft³
                T_low = T_low * 5 / 9      # Convert Rankine to Kelvin
                T_high = T_high * 5 / 9    # Convert Rankine to Kelvin
                V0 = V0 * 0.0283168        # Convert ft³ to m³
                P0 = P0 * 101325           # Convert atm to Pa
            else:  # Metric units: Kelvin, bar, m³
                P0 = P0 * 1e5              # Convert bar to Pa

            # Set model parameters
            self.model.r = float(inputs['CR'])      # Compression ratio
            self.model.rc = float(inputs['rc'])     # Cutoff ratio
            self.model.T1 = T_low                   # Initial temperature
            self.model.T_high = T_high              # Maximum temperature
            self.model.P1 = P0                      # Initial pressure
            self.model.V1 = V0                      # Initial volume

            # Perform calculations and update view
            self.model.calculate()
            self.view.updateResults(self.model.results)
            self.view.plotPV(self.model.results['pv'])  # Plot P-V diagram

        except ValueError as e:
            raise ValueError(f"Input error: {str(e)}")
        except Exception as e:
            raise Exception(f"Calculation error: {str(e)}")

    def updateView(self):
        """Updates the GUI with calculation results, units, and plots.

        Refreshes the results display, updates unit labels based on the selected unit
        system, and redraws the plot based on user-selected axes and log scale options.
        """
        widgets = self.view.widgets
        # Update displayed results
        self.view.updateResults(self.model.results)
        # Update unit labels based on metric/English selection
        self.view.updateUnits(widgets['rdo_Metric'].isChecked())
        # Redraw plot with current settings
        self.view.updatePlot(
            self.model.results,
            widgets['cmb_Abcissa'].currentText(),    # X-axis property
            widgets['cmb_Ordinate'].currentText(),  # Y-axis property
            widgets['chk_LogAbcissa'].isChecked(),  # Log scale for x-axis
            widgets['chk_LogOrdinate'].isChecked(), # Log scale for y-axis
            widgets['ax'],                          # Matplotlib axes
            widgets['canvas']                       # Matplotlib canvas
        )