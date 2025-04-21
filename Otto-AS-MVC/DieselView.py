from PyQt5 import QtWidgets as qtw


class DieselView:
    """Manages the GUI for displaying Diesel cycle results and plots.

    Handles widget updates, unit conversions for display, and P-V diagram plotting
    based on calculation results from the DieselModel.

    Attributes:
        widgets (dict): Dictionary mapping widget names to their Qt objects.
        ax (matplotlib.axes.Axes): Matplotlib axes for plotting.
        canvas (matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg): Matplotlib canvas.
        results (dict): Stored results from the DieselModel for redrawing.
    """

    def __init__(self):
        """Initializes the view with empty widgets and results."""
        self.widgets = {}        # Dictionary to store GUI widgets
        self.ax = None          # Matplotlib axes for plotting
        self.canvas = None      # Matplotlib canvas for rendering
        self.results = None     # Store results for redrawing

    def setWidgets(self, widget_dict):
        """Sets the GUI widgets and initializes plotting components.

        Args:
            widget_dict (dict): Dictionary mapping widget names to their Qt objects.
        """
        self.widgets = widget_dict              # Store widget dictionary
        self.ax = widget_dict['ax']             # Set matplotlib axes
        self.canvas = widget_dict['canvas']     # Set matplotlib canvas

    def updateResults(self, results):
        """Updates GUI widgets with Diesel cycle calculation results.

        Converts units based on the selected system (metric or English) and updates
        text fields with temperatures, work, heat, and efficiency.

        Args:
            results (dict): Dictionary containing cycle results (T1, T2, P1, P2, etc.).
        """
        self.results = results  # Store results for redrawing
        w = self.widgets
        is_metric = w['rdo_Metric'].isChecked()  # Check unit system

        # Convert temperatures from Kelvin to Rankine if English units
        T1 = results['T1'] * 9 / 5 if not is_metric else results['T1']
        T2 = results['T2'] * 9 / 5 if not is_metric else results['T2']
        T3 = results['T3'] * 9 / 5 if not is_metric else results['T3']
        T4 = results['T4'] * 9 / 5 if not is_metric else results['T4']

        # Convert work and heat from kJ/kg to BTU/lbm if English units
        factor = 1.0 if is_metric else 0.4299  # Conversion factor: kJ/kg to BTU/lbm
        w_exp = results['w_exp'] * factor      # Expansion work
        w_comp = results['w_comp'] * factor    # Compression work
        q_in = results['q_in'] * factor        # Heat input

        # Update text fields with formatted values
        w['le_T1'].setText(f"{T1:.2f}")
        w['le_T2'].setText(f"{T2:.2f}")
        w['le_T3'].setText(f"{T3:.2f}")
        w['le_T4'].setText(f"{T4:.2f}")
        w['le_HeatAdded'].setText(f"{q_in:.2f}")
        w['le_CompressionStroke'].setText(f"{w_comp:.2f}")
        w['le_PowerStroke'].setText(f"{w_exp:.2f}")
        w['le_Efficiency'].setText(f"{results['eff']:.2f}")

    def updateUnits(self, is_metric):
        """Updates unit labels in the GUI based on the selected unit system.

        Also refreshes results and plot to reflect unit changes.

        Args:
            is_metric (bool): True for metric units (K, bar, m³), False for English (R, atm, ft³).
        """
        w = self.widgets
        if is_metric:
            # Set metric unit labels: Kelvin, bar, m³
            w['lbl_T1Units'].setText("K")
            w['lbl_T2Units'].setText("K")
            w['lbl_T3Units'].setText("K")
            w['lbl_T4Units'].setText("K")
            w['lbl_PowerStrokeUnits'].setText("kJ/kg")
            w['lbl_CompressionStrokeUnits'].setText("kJ/kg")
            w['lbl_HeatInUnits'].setText("kJ/kg")
        else:
            # Set English unit labels: Rankine, atm, ft³
            w['lbl_T1Units'].setText("R")
            w['lbl_T2Units'].setText("R")
            w['lbl_T3Units'].setText("R")
            w['lbl_T4Units'].setText("R")
            w['lbl_PowerStrokeUnits'].setText("BTU/lbm")
            w['lbl_CompressionStrokeUnits'].setText("BTU/lbm")
            w['lbl_HeatInUnits'].setText("BTU/lbm")

        # Refresh results and plot with updated units if results exist
        if self.results:
            self.updateResults(self.results)
            self.updatePlot(
                self.results,
                w['cmb_Abcissa'].currentText(),    # X-axis property
                w['cmb_Ordinate'].currentText(),   # Y-axis property
                w['chk_LogAbcissa'].isChecked(),   # Log scale for x-axis
                w['chk_LogOrdinate'].isChecked(),  # Log scale for y-axis
                self.ax,
                self.canvas
            )

    def updatePlot(self, results, x_type, y_type, log_x, log_y, ax, canvas):
        """Updates the P-V diagram based on results and plot settings.

        Args:
            results (dict): Dictionary containing cycle results, including 'pv' data.
            x_type (str): X-axis property (not used in current P-V plot).
            y_type (str): Y-axis property (not used in current P-V plot).
            log_x (bool): True to use logarithmic x-axis.
            log_y (bool): True to use logarithmic y-axis.
            ax (matplotlib.axes.Axes): Matplotlib axes for plotting.
            canvas (matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg): Matplotlib canvas.
        """
        self.ax = ax          # Update axes
        self.canvas = canvas  # Update canvas
        self.plotPV(results['pv'], log_x, log_y)  # Plot P-V diagram

    def plotPV(self, pv_data, log_x=False, log_y=False):
        """Plots the P-V diagram for the Diesel cycle.

        Args:
            pv_data (list): List of (volume, pressure) tuples for the cycle.
            log_x (bool, optional): True to use logarithmic x-axis. Defaults to False.
            log_y (bool, optional): True to use logarithmic y-axis. Defaults to False.
        """
        self.ax.clear()  # Clear previous plot
        v, p = zip(*pv_data)  # Unpack volume and pressure data
        is_metric = self.widgets['rdo_Metric'].isChecked()  # Check unit system

        if is_metric:
            # Metric units: m³, Pa
            v_label = 'Volume (m³)'
            p_label = 'Pressure (Pa)'
        else:
            # English units: ft³, atm
            v = [vol / 0.0283168 for vol in v]  # Convert m³ to ft³
            p = [press / 101325 for press in p]  # Convert Pa to atm
            v_label = 'Volume (ft³)'
            p_label = 'Pressure (atm)'

        # Plot the Diesel cycle
        self.ax.plot(v, p, 'r-', label='Diesel Cycle')
        self.ax.set_xlabel(v_label)
        self.ax.set_ylabel(p_label)
        self.ax.set_title('Diesel PV Diagram')

        # Apply logarithmic scales if selected
        if log_x:
            self.ax.set_xscale('log')
        if log_y:
            self.ax.set_yscale('log')

        self.ax.grid(True)    # Add grid
        self.ax.legend()      # Show legend
        self.canvas.draw()    # Redraw canvas