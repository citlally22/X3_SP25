from Air import *
from matplotlib import pyplot as plt
from PyQt5 import QtWidgets as qtw
import sys
import numpy as np


class ottoCycleModel:
    """Models an air-standard Otto cycle for thermodynamic analysis.

    Represents the Otto cycle with four processes: isentropic compression, constant-volume
    heat addition, isentropic expansion, and constant-volume heat rejection. Uses the air
    class for thermodynamic properties.

    Attributes:
        units (units): Units object for handling SI or English units.
        air (air): Air model for thermodynamic calculations.
        p_initial (float): Initial pressure in Pa.
        T_initial (float): Initial temperature in Kelvin.
        T_high (float): Maximum temperature (state 3) in Kelvin.
        Ratio (float): Compression ratio (V1/V2).
        V_Cylinder (float): Cylinder volume at bottom dead center in m³.
        State1 (stateProps): Thermodynamic state at start of compression.
        State2 (stateProps): State after isentropic compression.
        State3 (stateProps): State after constant-volume heat addition.
        State4 (stateProps): State after isentropic expansion.
        W_Compression (float): Compression work in J.
        W_Power (float): Power stroke work in J.
        Q_In (float): Heat input in J.
        Q_Out (float): Heat rejection in J.
        W_Cycle (float): Net cycle work in J.
        Eff (float): Thermal efficiency (fraction).
        upperCurve (StateDataForPlotting): Data for upper cycle curve (states 2-3-4-1).
        lowerCurve (StateDataForPlotting): Data for lower cycle curve (states 1-2).
    """

    def __init__(self, p_initial=1000.0, v_cylinder=1.0, t_initial=298, t_high=1500.0, ratio=6.0, name='Air Standard Otto Cycle'):
        """Initializes the Otto cycle model with specified parameters.

        Args:
            p_initial (float, optional): Initial pressure in Pa. Defaults to 1000.0.
            v_cylinder (float, optional): Cylinder volume in m³. Defaults to 1.0.
            t_initial (float, optional): Initial temperature in Kelvin. Defaults to 298.
            t_high (float, optional): Maximum temperature in Kelvin. Defaults to 1500.0.
            ratio (float, optional): Compression ratio. Defaults to 6.0.
            name (str, optional): Cycle name. Defaults to 'Air Standard Otto Cycle'.
        """
        self.units = units()  # Initialize units object
        self.air = air()      # Initialize air model
        self.air.set(P=p_initial, T=t_initial)  # Set initial state
        self.p_initial = p_initial
        self.T_initial = t_initial
        self.T_high = t_high
        self.Ratio = ratio
        self.V_Cylinder = v_cylinder
        self.air.n = self.V_Cylinder / self.air.State.v  # Calculate moles
        self.air.m = self.air.n * self.air.MW / 1000.0   # Calculate mass in kg

        # Define the four states of the Otto cycle
        self.State1 = self.air.set(P=self.p_initial, T=self.T_initial)  # State 1: Start
        self.State2 = self.air.set(v=self.State1.v/self.Ratio, s=self.State1.s)  # State 2: Isentropic compression
        self.State3 = self.air.set(T=self.T_high, v=self.State2.v)  # State 3: Constant-volume heat addition
        self.State4 = self.air.set(v=self.State1.v, s=self.State3.s)  # State 4: Isentropic expansion

        # Calculate work and heat (extensive, in J)
        self.W_Compression = self.air.n * (self.State2.u - self.State1.u)  # Compression work
        self.W_Power = self.air.n * (self.State3.u - self.State4.u)  # Power stroke work
        self.Q_In = self.air.n * (self.State3.u - self.State2.u)  # Heat input
        self.Q_Out = self.air.n * (self.State4.u - self.State1.u)  # Heat rejection

        # Calculate net work and efficiency
        self.W_Cycle = self.W_Power - self.W_Compression
        self.Eff = self.W_Cycle / self.Q_In if self.Q_In != 0 else 0.0

        # Initialize data for plotting
        self.upperCurve = StateDataForPlotting()  # States 2-3-4-1
        self.lowerCurve = StateDataForPlotting()  # States 1-2

    def getSI(self):
        """Returns whether SI units are used.

        Returns:
            bool: True if SI units, False if English units.
        """
        return self.units.SI


class ottoCycleController:
    """Controller for managing Otto cycle calculations and GUI updates.

    Coordinates between the OttoCycleModel for calculations and OttoCycleView for
    GUI interactions, handling input validation, plotting, and result display.

    Attributes:
        model (ottoCycleModel): Otto cycle model instance.
        view (ottoCycleView): View instance for GUI updates.
        widgets (dict): Dictionary mapping widget names to their Qt objects.
    """

    def __init__(self, model=None, ax=None):
        """Initializes the controller with a model and view.

        Args:
            model (ottoCycleModel, optional): Otto cycle model. Defaults to None (creates new).
            ax (matplotlib.axes.Axes, optional): Matplotlib axes for plotting. Defaults to None.
        """
        self.model = ottoCycleModel() if model is None else model  # Initialize model
        self.view = ottoCycleView()  # Initialize view
        self.view.ax = ax  # Set axes for plotting
        self.widgets = None  # Initialize widgets dictionary

    def calc(self, T_0=None, P_0=None, V_0=None, T_High=None, CR=None, SI=True):
        """Performs Otto cycle calculations based on user inputs.

        Validates inputs, updates the model, and triggers calculations. Supports both
        GUI inputs (via widgets) and direct parameter inputs.

        Args:
            T_0 (float, optional): Initial temperature (K or R). Defaults to None.
            P_0 (float, optional): Initial pressure (Pa or atm). Defaults to None.
            V_0 (float, optional): Cylinder volume (m³ or ft³). Defaults to None.
            T_High (float, optional): Maximum temperature (K or R). Defaults to None.
            CR (float, optional): Compression ratio. Defaults to None.
            SI (bool, optional): True for SI units, False for English. Defaults to True.

        Raises:
            ValueError: If inputs are invalid (non-numeric, non-positive).
            Exception: If calculations fail.
        """
        if self.widgets and all(x is None for x in [T_0, P_0, V_0, T_High, CR]):
            # Use GUI inputs from widgets
            w = self.widgets
            try:
                T0 = float(w['le_TLow'].text())
                P0 = float(w['le_P0'].text())
                V0 = float(w['le_V0'].text())
                TH = float(w['le_THigh'].text())
                CR = float(w['le_CR'].text())
                metric = w['rdo_Metric'].isChecked()
            except ValueError as e:
                raise ValueError(f"Invalid input: {str(e)}")
        else:
            # Use provided parameters
            if any(x is None for x in [T_0, P_0, V_0, T_High, CR]):
                raise ValueError("All parameters (T_0, P_0, V_0, T_High, CR) must be provided")
            T0, P0, V0, TH, CR = T_0, P_0, V_0, T_High, CR
            metric = SI

        # Convert units to SI if necessary
        try:
            if not metric:
                T0 = T0 * 5 / 9  # Rankine to Kelvin
                TH = TH * 5 / 9  # Rankine to Kelvin
                P0 = P0 * 101325  # atm to Pa
                V0 = V0 * 0.0283168  # ft³ to m³
            self.model = ottoCycleModel(p_initial=P0, v_cylinder=V0, t_initial=T0, t_high=TH, ratio=CR)
            self.model.units.set(SI=metric)  # Set unit system
            self.buildDataForPlotting()  # Generate plotting data
        except Exception as e:
            raise Exception(f"Calculation error: {str(e)}")

    def setWidgets(self, widgets=None):
        """Sets the GUI widgets for the controller and view.

        Args:
            widgets (dict, optional): Dictionary mapping widget names to Qt objects. Defaults to None.
        """
        self.widgets = widgets
        if widgets:
            # Assign widgets to view for GUI updates
            self.view.lbl_THigh = widgets['lbl_THigh']
            self.view.lbl_TLow = widgets['lbl_TLow']
            self.view.lbl_P0 = widgets['lbl_P0']
            self.view.lbl_V0 = widgets['lbl_V0']
            self.view.lbl_CR = widgets['lbl_CR']
            self.view.le_THigh = widgets['le_THigh']
            self.view.le_TLow = widgets['le_TLow']
            self.view.le_P0 = widgets['le_P0']
            self.view.le_V0 = widgets['le_V0']
            self.view.le_CR = widgets['le_CR']
            self.view.le_T1 = widgets['le_T1']
            self.view.le_T2 = widgets['le_T2']
            self.view.le_T3 = widgets['le_T3']
            self.view.le_T4 = widgets['le_T4']
            self.view.lbl_T1Units = widgets['lbl_T1Units']
            self.view.lbl_T2Units = widgets['lbl_T2Units']
            self.view.lbl_T3Units = widgets['lbl_T3Units']
            self.view.lbl_T4Units = widgets['lbl_T4Units']
            self.view.le_PowerStroke = widgets['le_PowerStroke']
            self.view.le_CompressionStroke = widgets['le_CompressionStroke']
            self.view.le_HeatAdded = widgets['le_HeatAdded']
            self.view.le_Efficiency = widgets['le_Efficiency']
            self.view.lbl_PowerStrokeUnits = widgets['lbl_PowerStrokeUnits']
            self.view.lbl_CompressionStrokeUnits = widgets['lbl_CompressionStrokeUnits']
            self.view.lbl_HeatInUnits = widgets['lbl_HeatInUnits']
            self.view.rdo_Metric = widgets['rdo_Metric']
            self.view.cmb_Abcissa = widgets['cmb_Abcissa']
            self.view.cmb_Ordinate = widgets['cmb_Ordinate']
            self.view.chk_LogAbcissa = widgets['chk_LogAbcissa']
            self.view.chk_LogOrdinate = widgets['chk_LogOrdinate']
            self.view.ax = widgets['ax']
            self.view.canvas = widgets['canvas']

    def buildDataForPlotting(self):
        """Generates state data for plotting the Otto cycle.

        Creates data for the upper curve (states 2-3-4-1) and lower curve (states 1-2)
        to represent the cycle on a thermodynamic diagram (e.g., P-V, T-s).
        """
        self.model.upperCurve.clear()  # Clear previous upper curve data
        self.model.lowerCurve.clear()  # Clear previous lower curve data
        a = air()  # Temporary air instance for state calculations

        # Region: States from 2-3 (constant volume, T from T2 to T3)
        DeltaT = np.linspace(self.model.State2.T, self.model.State3.T, 30)
        for T in DeltaT:
            state = a.set(T=T, v=self.model.State2.v)  # Constant volume
            self.model.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))

        # Region: States from 3-4 (isentropic expansion, v from v3 to v4)
        DeltaV = np.linspace(self.model.State3.v, self.model.State4.v, 30)
        for v in DeltaV:
            state = a.set(v=v, s=self.model.State3.s)  # Constant entropy
            self.model.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))

        # Region: States from 4-1 (constant volume, T from T4 to T1)
        DeltaT = np.linspace(self.model.State4.T, self.model.State1.T, 30)
        for T in DeltaT:
            state = a.set(T=T, v=self.model.State4.v)  # Constant volume
            self.model.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))

        # Region: States from 1-2 (isentropic compression, v from v1 to v2)
        DeltaV = np.linspace(self.model.State1.v, self.model.State2.v, 30)
        for v in DeltaV:
            state = a.set(v=v, s=self.model.State1.s)  # Constant entropy
            self.model.lowerCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))

    def plot_cycle_XY(self, X='s', Y='T', logx=False, logy=False, mass=False, total=False):
        """Plots the Otto cycle with specified thermodynamic properties.

        Args:
            X (str, optional): X-axis property ('T', 'P', 'v', 's', 'h', 'u'). Defaults to 's'.
            Y (str, optional): Y-axis property ('T', 'P', 'v', 's', 'h', 'u'). Defaults to 'T'.
            logx (bool, optional): True for logarithmic x-axis. Defaults to False.
            logy (bool, optional): True for logarithmic y-axis. Defaults to False.
            mass (bool, optional): True for mass-based units. Defaults to False.
            total (bool, optional): True for total (extensive) properties. Defaults to False.
        """
        self.view.plot_cycle_XY(self.model, X=X, Y=Y, logx=logx, logy=logy, mass=mass, total=total)

    def print_summary(self):
        """Prints a summary of the Otto cycle results."""
        self.view.print_summary(self.model)

    def get_summary(self):
        """Returns a formatted string summary of the Otto cycle.

        Returns:
            str: Summary string with efficiency, work, and heat values.
        """
        return self.view.get_summary(self.model)

    def updateView(self):
        """Updates the GUI with the current model results and plot."""
        self.view.updateView(cycle=self.model)


class ottoCycleView:
    """Manages the GUI for displaying Otto cycle results and plots.

    Updates widgets with cycle results, handles unit conversions, and plots
    thermodynamic diagrams (e.g., P-V, T-s) using matplotlib.

    Attributes:
        lbl_THigh (QLabel): Label for maximum temperature input.
        lbl_TLow (QLabel): Label for initial temperature input.
        lbl_P0 (QLabel): Label for initial pressure input.
        lbl_V0 (QLabel): Label for cylinder volume input.
        lbl_CR (QLabel): Label for compression ratio input.
        le_THigh (QLineEdit): Input field for maximum temperature.
        le_TLow (QLineEdit): Input field for initial temperature.
        le_P0 (QLineEdit): Input field for initial pressure.
        le_V0 (QLineEdit): Input field for cylinder volume.
        le_CR (QLineEdit): Input field for compression ratio.
        le_T1 (QLineEdit): Display field for state 1 temperature.
        le_T2 (QLineEdit): Display field for state 2 temperature.
        le_T3 (QLineEdit): Display field for state 3 temperature.
        le_T4 (QLineEdit): Display field for state 4 temperature.
        lbl_T1Units (QLabel): Unit label for state 1 temperature.
        lbl_T2Units (QLabel): Unit label for state 2 temperature.
        lbl_T3Units (QLabel): Unit label for state 3 temperature.
        lbl_T4Units (QLabel): Unit label for state 4 temperature.
        le_Efficiency (QLineEdit): Display field for cycle efficiency.
        le_PowerStroke (QLineEdit): Display field for power stroke work.
        le_CompressionStroke (QLineEdit): Display field for compression work.
        le_HeatAdded (QLineEdit): Display field for heat input.
        lbl_PowerStrokeUnits (QLabel): Unit label for power stroke work.
        lbl_CompressionStrokeUnits (QLabel): Unit label for compression work.
        lbl_HeatInUnits (QLabel): Unit label for heat input.
        rdo_Metric (QRadioButton): Radio button for metric units.
        cmb_Abcissa (QComboBox): Combo box for x-axis property.
        cmb_Ordinate (QComboBox): Combo box for y-axis property.
        chk_LogAbcissa (QCheckBox): Checkbox for logarithmic x-axis.
        chk_LogOrdinate (QCheckBox): Checkbox for logarithmic y-axis.
        canvas (FigureCanvasQTAgg): Matplotlib canvas for plotting.
        ax (matplotlib.axes.Axes): Matplotlib axes for plotting.
    """

    def __init__(self):
        """Initializes the view with Qt widgets and plotting components."""
        # Initialize input labels
        self.lbl_THigh = qtw.QLabel()
        self.lbl_TLow = qtw.QLabel()
        self.lbl_P0 = qtw.QLabel()
        self.lbl_V0 = qtw.QLabel()
        self.lbl_CR = qtw.QLabel()
        # Initialize input fields
        self.le_THigh = qtw.QLineEdit()
        self.le_TLow = qtw.QLineEdit()
        self.le_P0 = qtw.QLineEdit()
        self.le_V0 = qtw.QLineEdit()
        self.le_CR = qtw.QLineEdit()
        # Initialize result display fields
        self.le_T1 = qtw.QLineEdit()
        self.le_T2 = qtw.QLineEdit()
        self.le_T3 = qtw.QLineEdit()
        self.le_T4 = qtw.QLineEdit()
        # Initialize unit labels
        self.lbl_T1Units = qtw.QLabel()
        self.lbl_T2Units = qtw.QLabel()
        self.lbl_T3Units = qtw.QLabel()
        self.lbl_T4Units = qtw.QLabel()
        # Initialize result fields
        self.le_Efficiency = qtw.QLineEdit()
        self.le_PowerStroke = qtw.QLineEdit()
        self.le_CompressionStroke = qtw.QLineEdit()
        self.le_HeatAdded = qtw.QLineEdit()
        # Initialize unit labels for results
        self.lbl_PowerStrokeUnits = qtw.QLabel()
        self.lbl_CompressionStrokeUnits = qtw.QLabel()
        self.lbl_HeatInUnits = qtw.QLabel()
        # Initialize unit selection and plot controls
        self.rdo_Metric = qtw.QRadioButton()
        self.cmb_Abcissa = qtw.QComboBox()
        self.cmb_Ordinate = qtw.QComboBox()
        self.chk_LogAbcissa = qtw.QCheckBox()
        self.chk_LogOrdinate = qtw.QCheckBox()
        # Initialize plotting components
        self.canvas = None
        self.ax = None

    def updateView(self, cycle):
        """Updates the GUI with cycle results and plot.

        Args:
            cycle (ottoCycleModel): Otto cycle model with calculation results.
        """
        cycle.units.SI = self.rdo_Metric.isChecked()  # Set unit system
        logx = self.chk_LogAbcissa.isChecked()       # Check x-axis log scale
        logy = self.chk_LogOrdinate.isChecked()      # Check y-axis log scale
        xvar = self.cmb_Abcissa.currentText()        # Get x-axis property
        yvar = self.cmb_Ordinate.currentText()       # Get y-axis property
        # Plot cycle with selected properties
        self.plot_cycle_XY(cycle, X=xvar, Y=yvar, logx=logx, logy=logy, mass=False, total=True)
        # Update display widgets with results
        self.updateDisplayWidgets(Model=cycle)

    def print_summary(self, cycle):
        """Prints a summary of the Otto cycle results to the console.

        Args:
            cycle (ottoCycleModel): Otto cycle model with calculation results.
        """
        print('Cycle Summary for:', cycle.name)
        print(f'\tEfficiency: {cycle.Eff*100:.3f}%')  # Show as percentage
        print(f'\tPower Stroke: {cycle.W_Power/1000:.3f} kJ/kg')  # Convert J to kJ
        print(f'\tCompression Stroke: {cycle.W_Compression/1000:.3f} kJ/kg')  # Convert J to kJ
        print(f'\tHeat Added: {cycle.Q_In/1000:.3f} kJ/kg')  # Convert J to kJ
        # Print state details
        cycle.State1.print()
        cycle.State2.print()
        cycle.State3.print()
        cycle.State4.print()

    def get_summary(self, cycle):
        """Returns a formatted string summary for the Otto cycle.

        Args:
            cycle (ottoCycleModel): Otto cycle model with calculation results.

        Returns:
            str: Formatted string with efficiency, work, and heat values.
        """
        s = r'Summary:'
        s += f'\nEfficiency: {cycle.Eff*100:.1f}%'  # Show as percentage
        s += f'\nPower Stroke: {cycle.W_Power*cycle.air.n/1000:.1f} kJ'  # Extensive, in kJ
        s += f'\nCompression Stroke: {cycle.W_Compression*cycle.air.n/1000:.1f} kJ'  # Extensive, in kJ
        s += f'\nHeat Added: {cycle.Q_In*cycle.air.n/1000:.1f} kJ'  # Extensive, in kJ
        return s

    def convertDataCol(self, cycle, data=None, colName='T', mass=False, total=False):
        """Converts thermodynamic data to the desired unit system and basis.

        Args:
            cycle (ottoCycleModel): Otto cycle model with units and air properties.
            data (list, optional): List of data values to convert. Defaults to None.
            colName (str, optional): Property name ('T', 'P', 'v', 's', 'h', 'u'). Defaults to 'T'.
            mass (bool, optional): True for mass-based units. Defaults to False.
            total (bool, optional): True for total (extensive) properties. Defaults to False.

        Returns:
            list: Converted data values.
        """
        UC = cycle.units
        n = cycle.air.n
        MW = cycle.air.MW / 1000.0  # Convert kg/kmol to kg/mol
        # Define conversion factors
        TCF = 1.0 if UC.SI else UC.CF_T
        PCF = 1.0 if UC.SI else UC.CF_P
        hCF = 1.0 if UC.SI else UC.CF_e
        uCF = 1.0 if UC.SI else UC.CF_e
        sCF = 1.0 if UC.SI else UC.CF_s
        vCF = 1.0 if UC.SI else UC.CF_v
        nCF = 1.0 if UC.SI else UC.CF_n
        mCF = 1.0 if UC.SI else UC.CF_Mass
        if mass:
            hCF /= MW
            uCF /= MW
            sCF /= MW
            vCF /= MW
        elif total:
            hCF *= n * nCF
            uCF *= n * nCF
            sCF *= n * nCF
            vCF *= n * nCF
        # Convert data based on property
        w = colName.lower()
        if w == 't':
            return [T * TCF for T in data]
        if w == 'h':
            return [h * hCF for h in data]
        if w == 'u':
            return [u * uCF for h in data]
        if w == 's':
            return [s * sCF for s in data]
        if w == 'v':
            return [v * vCF for v in data]
        if w == 'p':
            return [P * PCF for P in data]
        return data  # Return unchanged if property not recognized

    def plot_cycle_XY(self, cycle, X='s', Y='T', logx=False, logy=False, mass=False, total=False):
        """Plots the Otto cycle with specified thermodynamic properties.

        Args:
            cycle (ottoCycleModel): Otto cycle model with calculation results.
            X (str, optional): X-axis property ('T', 'P', 'v', 's', 'h', 'u'). Defaults to 's'.
            Y (str, optional): Y-axis property ('T', 'P', 'v', 's', 'h', 'u'). Defaults to 'T'.
            logx (bool, optional): True for logarithmic x-axis. Defaults to False.
            logy (bool, optional): True for logarithmic y-axis. Defaults to False.
            mass (bool, optional): True for mass-based units. Defaults to False.
            total (bool, optional): True for total (extensive) properties. Defaults to False.
        """
        if X == Y:
            return  # Avoid plotting same property on both axes

        QTPlotting = True
        if self.ax is None:
            self.ax = plt.subplot()  # Create new axes if none provided
            QTPlotting = False

        ax = self.ax
        ax.clear()  # Clear previous plot

        # Set axis scales
        ax.set_xscale('log' if logx else 'linear')
        ax.set_yscale('log' if logy else 'linear')

        # Convert and plot upper and lower curves
        XdataLC = self.convertDataCol(cycle, colName=X, data=cycle.lowerCurve.getDataCol(X), mass=mass, total=total)
        YdataLC = self.convertDataCol(cycle, colName=Y, data=cycle.lowerCurve.getDataCol(Y), mass=mass, total=total)
        XdataUC = self.convertDataCol(cycle, colName=X, data=cycle.upperCurve.getDataCol(X), mass=mass, total=total)
        YdataUC = self.convertDataCol(cycle, colName=Y, data=cycle.upperCurve.getDataCol(Y), mass=mass, total=total)
        ax.plot(XdataLC, YdataLC, color='k', label='Compression (1-2)')  # Lower curve
        ax.plot(XdataUC, YdataUC, color='g', label='Upper Cycle (2-3-4-1)')  # Upper curve

        # Set axis labels
        cycle.units.setPlotUnits(SI=cycle.units.SI, mass=mass, total=total)
        ax.set_ylabel(cycle.lowerCurve.getAxisLabel(Y, Units=cycle.units), fontsize='large')
        ax.set_xlabel(cycle.lowerCurve.getAxisLabel(X, Units=cycle.units), fontsize='large')

        # Set plot title
        cycle.name = 'Otto Cycle'
        ax.set_title(cycle.name, fontsize='large')

        # Customize tick marks
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize='large')

        # Plot state points (1, 2, 3, 4)
        state1 = dc(cycle.State1)
        state1.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW/1000.0, mass=mass, total=total)
        state2 = dc(cycle.State2)
        state2.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW/1000.0, mass=mass, total=total)
        state3 = dc(cycle.State3)
        state3.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW/1000.0, mass=mass, total=total)
        state4 = dc(cycle.State4)
        state4.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW/1000.0, mass=mass, total=total)

        ax.plot(state1.getVal(X), state1.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k', label='State 1')
        ax.plot(state2.getVal(X), state2.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k', label='State 2')
        ax.plot(state3.getVal(X), state3.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k', label='State 3')
        ax.plot(state4.getVal(X), state4.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k', label='State 4')

        ax.legend()  # Show legend
        ax.grid(True)  # Add grid

        if QTPlotting:
            self.canvas.draw()  # Redraw canvas for Qt
        else:
            plt.show()  # Show plot for non-Qt

    def updateDisplayWidgets(self, Model=None):
        """Updates GUI widgets with Otto cycle results.

        Args:
            Model (ottoCycleModel, optional): Otto cycle model with results. Defaults to None.
        """
        if Model is None:
            return
        U = Model.units
        SI = U.SI
        CFT = 1.0 if SI else U.CF_T  # Temperature conversion factor
        CFE = 1.0 if SI else U.CF_E  # Energy conversion factor

        # Update input labels with units
        self.lbl_THigh.setText(f'T High ({U.TUnits})')
        self.lbl_TLow.setText(f'T Low ({U.TUnits})')
        self.lbl_P0.setText(f'P0 ({U.PUnits})')
        self.lbl_V0.setText(f'V0 ({U.VUnits})')
        # Update state temperatures
        self.le_T1.setText(f'{Model.State1.T * CFT:.2f}')
        self.le_T2.setText(f'{Model.State2.T * CFT:.2f}')
        self.le_T3.setText(f'{Model.State3.T * CFT:.2f}')
        self.le_T4.setText(f'{Model.State4.T * CFT:.2f}')
        # Update unit labels
        self.lbl_T1Units.setText(U.TUnits)
        self.lbl_T2Units.setText(U.TUnits)
        self.lbl_T3Units.setText(U.TUnits)
        self.lbl_T4Units.setText(U.TUnits)
        # Update results
        self.le_Efficiency.setText(f'{Model.Eff * 100:.3f}')  # Show as percentage
        self.le_PowerStroke.setText(f'{Model.air.n * Model.W_Power * CFE / 1000:.3f}')  # Convert to kJ
        self.le_CompressionStroke.setText(f'{Model.air.n * Model.W_Compression * CFE / 1000:.3f}')  # Convert to kJ
        self.le_HeatAdded.setText(f'{Model.air.n * Model.Q_In * CFE / 1000:.3f}')  # Convert to kJ
        # Update result unit labels
        self.lbl_PowerStrokeUnits.setText(U.EUnits)
        self.lbl_CompressionStrokeUnits.setText(U.EUnits)
        self.lbl_HeatInUnits.setText(U.EUnits)


def main():
    """Demonstrates the Otto cycle controller with a sample calculation and plot."""
    app = qtw.QApplication(sys.argv)  # Initialize Qt application
    oc = ottoCycleController()        # Create controller
    # Calculate cycle with sample parameters (English units)
    oc.calc(T_0=540.0, P_0=1.0, V_0=0.02, T_High=3600.0, CR=8.0, SI=False)
    # Plot P-V diagram with total properties
    oc.plot_cycle_XY(X='v', Y='P', total=True)
    sys.exit(app.exec_())  # Run Qt event loop


if __name__ == "__main__":
    main()