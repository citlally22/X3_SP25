from Air import *
from matplotlib import pyplot as plt
from PyQt5 import QtWidgets as qtw
import sys
import numpy as np

class ottoCycleModel:
    def __init__(self, p_initial=1000.0, v_cylinder=1.0, t_initial=298, t_high=1500.0, ratio=6.0, name='Air Standard Otto Cycle'):
        """
        Constructor for an air standard Otto cycle. The Otto has 4 primary states and consists of four thermodynamic
        processes:
        1. Isentropic compression from: v1, T1, P1 to v2, T2, P2 (Note v2=v1/C.R.)
        2. Constant volume heat addition: v3=v2
        3. Isentropic expansion (power stroke): v3=v1
        4. Constant volume heat rejection.
        Compression stroke work = (u2-u1)
        Power stroke work = (u3-u4)
        Heat in = (u3-u2)
        Heat out = (u4-u1)
        """
        self.units = units()
        self.air = air()
        self.air.set(P=p_initial, T=t_initial)
        self.p_initial = p_initial
        self.T_initial = t_initial
        self.T_high = t_high
        self.Ratio = ratio
        self.V_Cylinder = v_cylinder
        self.air.n = self.V_Cylinder / self.air.State.v
        self.air.m = self.air.n * self.air.MW

        self.State1 = self.air.set(P=self.p_initial, T=self.T_initial)
        self.State2 = self.air.set(v=self.State1.v/self.Ratio, s=self.State1.s)
        self.State3 = self.air.set(T=self.T_high, v=self.State2.v)
        self.State4 = self.air.set(v=self.State1.v, s=self.State3.s)

        self.W_Compression = self.air.n * (self.State2.u - self.State1.u)
        self.W_Power = self.air.n * (self.State3.u - self.State4.u)  # Fixed: Use State4 instead of State3
        self.Q_In = self.air.n * (self.State3.u - self.State2.u)
        self.Q_Out = self.air.n * (self.State4.u - self.State1.u)

        self.W_Cycle = self.W_Power - self.W_Compression
        self.Eff = self.W_Cycle / self.Q_In

        self.upperCurve = StateDataForPlotting()
        self.lowerCurve = StateDataForPlotting()

    def getSI(self):
        return self.units.SI

class ottoCycleController:
    def __init__(self, model=None, ax=None):
        self.model = ottoCycleModel() if model is None else model
        self.view = ottoCycleView()
        self.view.ax = ax

    def calc(self):
        w = self.widgets
        try:
            T0 = float(w['le_TLow'].text())
            P0 = float(w['le_P0'].text())
            V0 = float(w['le_V0'].text())
            TH = float(w['le_THigh'].text())
            CR = float(w['le_CR'].text())
            metric = w['rdo_Metric'].isChecked()
            self.model = ottoCycleModel(p_initial=P0, v_cylinder=V0, t_initial=T0, t_high=TH, ratio=CR)
            self.model.units.set(SI=metric)
        except ValueError as e:
            raise ValueError(f"Invalid input: {str(e)}")
        except Exception as e:
            raise Exception(f"Calculation error: {str(e)}")

    def setWidgets(self, widgets=None):
        self.widgets = widgets
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
        """
        I want to create state data between states 1-2, 2-3, 3-4, 4-1
        I'll piece together an upperCurve data set from 2-3, 3-4, 4-1
        The lowerCurve data set is 1-2
        """
        self.model.upperCurve.clear()
        self.model.lowerCurve.clear()
        a = air()
        # Region: states from 2-3 (v=const, T from T2->T3)
        DeltaT = np.linspace(self.model.State2.T, self.model.State3.T, 30)
        for T in DeltaT:
            state = a.set(T=T, v=self.model.State2.v)
            self.model.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))
        # Region: states from 3-4 (v=from TDC to BDC, s=const.)
        DeltaV = np.linspace(self.model.State3.v, self.model.State4.v, 30)
        for v in DeltaV:
            state = a.set(v=v, s=self.model.State3.s)
            self.model.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))
        # Region: states from 4-1 (v=const, T from T4->T1)
        DeltaT = np.linspace(self.model.State4.T, self.model.State1.T, 30)
        for T in DeltaT:
            state = a.set(T=T, v=self.model.State4.v)
            self.model.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))
        # Region: states from 1-2 (v=from BDC to TDC, s=const.)
        DeltaV = np.linspace(self.model.State1.v, self.model.State2.v, 30)
        for v in DeltaV:
            state = a.set(v=v, s=self.model.State1.s)
            self.model.lowerCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))

    def plot_cycle_XY(self, X='s', Y='T', logx=False, logy=False, mass=False, total=False):
        self.view.plot_cycle_XY(self.model, X=X, Y=Y, logx=logx, logy=logy, mass=mass, total=total)

    def print_summary(self):
        self.view.print_summary(self.model)

    def get_summary(self):
        return self.view.get_summary(self.model)

    def updateView(self):
        self.view.updateView(cycle=self.model)

class ottoCycleView:
    def __init__(self):
        self.lbl_THigh = qtw.QLabel()
        self.lbl_TLow = qtw.QLabel()
        self.lbl_P0 = qtw.QLabel()
        self.lbl_V0 = qtw.QLabel()
        self.lbl_CR = qtw.QLabel()
        self.le_THigh = qtw.QLineEdit()
        self.le_TLow = qtw.QLineEdit()
        self.le_P0 = qtw.QLineEdit()
        self.le_V0 = qtw.QLineEdit()
        self.le_CR = qtw.QLineEdit()
        self.le_T1 = qtw.QLineEdit()
        self.le_T2 = qtw.QLineEdit()
        self.le_T3 = qtw.QLineEdit()
        self.le_T4 = qtw.QLineEdit()
        self.lbl_T1Units = qtw.QLabel()
        self.lbl_T2Units = qtw.QLabel()
        self.lbl_T3Units = qtw.QLabel()
        self.lbl_T4Units = qtw.QLabel()
        self.le_Efficiency = qtw.QLineEdit()
        self.le_PowerStroke = qtw.QLineEdit()
        self.le_CompressionStroke = qtw.QLineEdit()
        self.le_HeatAdded = qtw.QLineEdit()
        self.lbl_PowerStrokeUnits = qtw.QLabel()
        self.lbl_CompressionStrokeUnits = qtw.QLabel()
        self.lbl_HeatInUnits = qtw.QLabel()
        self.rdo_Metric = qtw.QRadioButton()
        self.cmb_Abcissa = qtw.QComboBox()
        self.cmb_Ordinate = qtw.QComboBox()
        self.chk_LogAbcissa = qtw.QCheckBox()
        self.chk_LogOrdinate = qtw.QCheckBox()
        self.canvas = None
        self.ax = None

    def updateView(self, cycle):
        cycle.units.SI = self.rdo_Metric.isChecked()
        logx = self.chk_LogAbcissa.isChecked()
        logy = self.chk_LogOrdinate.isChecked()
        xvar = self.cmb_Abcissa.currentText()
        yvar = self.cmb_Ordinate.currentText()
        self.plot_cycle_XY(cycle, X=xvar, Y=yvar, logx=logx, logy=logy, mass=False, total=True)
        self.updateDisplayWidgets(Model=cycle)

    def print_summary(self, cycle):
        print('Cycle Summary for: ', cycle.name)
        print('\tEfficiency: {:0.3f}%'.format(cycle.Eff))
        print('\tPower Stroke: {:0.3f} kJ/kg'.format(cycle.W_Power))
        print('\tCompression Stroke: {:0.3f} kJ/kg'.format(cycle.W_Compression))
        print('\tHeat Added: {:0.3f} kJ/kg'.format(cycle.Q_In))
        cycle.State1.print()
        cycle.State2.print()
        cycle.State3.print()
        cycle.State4.print()

    def get_summary(self, cycle):
        '''
        This returns a formatted string to put on the plot of the Otto cycle.
        '''
        s = r'Summary:'
        s += '\nEfficiency: {:0.1f}%'.format(cycle.Eff * 100)  # Adjusted to show percentage
        s += '\nPower Stroke: {:0.1f} kJ'.format(cycle.W_Power * cycle.air.n)
        s += '\nCompression Stroke: {:0.1f} kJ'.format(cycle.W_Compression * cycle.air.n)
        s += '\nHeat Added: {:0.1f} kJ'.format(cycle.Q_In * cycle.air.n)
        return s

    def convertDataCol(self, cycle, data=None, colName='T', mass=False, total=False):
        UC = cycle.units
        n = cycle.air.n
        MW = cycle.air.MW
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
        w = colName.lower()
        if w == 't':
            return [T * TCF for T in data]
        if w == 'h':
            return [h * hCF for h in data]
        if w == 'u':
            return [u * uCF for u in data]
        if w == 's':
            return [s * sCF for s in data]
        if w == 'v':
            return [v * vCF for v in data]
        if w == 'p':
            return [P * PCF for P in data]

    def plot_cycle_XY(self, cycle, X='s', Y='T', logx=False, logy=False, mass=False, total=False):
        """
        I want to plot any two thermodynamic properties on X and Y
        Data is in molar metric units. I may need to convert it.
        """
        if X == Y:
            return
        QTPlotting = True
        if self.ax is None:
            self.ax = plt.subplot()
            QTPlotting = False

        ax = self.ax
        ax.clear()

        ax.set_xscale('log' if logx else 'linear')
        ax.set_yscale('log' if logy else 'linear')

        # Plot the upper and lower curves
        XdataLC = self.convertDataCol(cycle, colName=X, data=cycle.lowerCurve.getDataCol(X), mass=mass, total=total)
        YdataLC = self.convertDataCol(cycle, colName=Y, data=cycle.lowerCurve.getDataCol(Y), mass=mass, total=total)
        XdataUC = self.convertDataCol(cycle, colName=X, data=cycle.upperCurve.getDataCol(X), mass=mass, total=total)
        YdataUC = self.convertDataCol(cycle, colName=Y, data=cycle.upperCurve.getDataCol(Y), mass=mass, total=total)
        ax.plot(XdataLC, YdataLC, color='k')
        ax.plot(XdataUC, YdataUC, color='g')

        # Add axis labels
        cycle.units.setPlotUnits(SI=cycle.units.SI, mass=mass, total=total)
        ax.set_ylabel(cycle.lowerCurve.getAxisLabel(Y, Units=cycle.units), fontsize='large')
        ax.set_xlabel(cycle.lowerCurve.getAxisLabel(X, Units=cycle.units), fontsize='large')

        # Put a title on the plot
        cycle.name = 'Otto Cycle'
        ax.set_title(cycle.name, fontsize='large')

        # Modify the tick marks
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize='large')

        # Plot the circles for states 1, 2, 3, and 4
        state1 = dc(cycle.State1)
        state1.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW, mass=mass, total=total)
        state2 = dc(cycle.State2)
        state2.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW, mass=mass, total=total)
        state3 = dc(cycle.State3)
        state3.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW, mass=mass, total=total)
        state4 = dc(cycle.State4)
        state4.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW, mass=mass, total=total)

        ax.plot(state1.getVal(X), state1.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(state2.getVal(X), state2.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(state3.getVal(X), state3.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(state4.getVal(X), state4.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k')

        if QTPlotting:
            self.canvas.draw()
        else:
            plt.show()

    def updateDisplayWidgets(self, Model=None):
        U = Model.units
        SI = U.SI
        CFT = 1.0 if SI else U.CF_T
        CFE = 1.0 if SI else U.CF_E

        self.lbl_THigh.setText('T High ({})'.format(Model.units.TUnits))
        self.lbl_TLow.setText('T Low ({})'.format(Model.units.TUnits))
        self.lbl_P0.setText('P0 ({})'.format(Model.units.PUnits))
        self.lbl_V0.setText('V0 ({})'.format(Model.units.VUnits))
        self.le_T1.setText('{:0.2f}'.format(Model.State1.T * CFT))
        self.le_T2.setText('{:0.2f}'.format(Model.State2.T * CFT))
        self.le_T3.setText('{:0.2f}'.format(Model.State3.T * CFT))
        self.le_T4.setText('{:0.2f}'.format(Model.State4.T * CFT))
        self.lbl_T1Units.setText(Model.units.TUnits)
        self.lbl_T2Units.setText(Model.units.TUnits)
        self.lbl_T3Units.setText(Model.units.TUnits)
        self.lbl_T4Units.setText(Model.units.TUnits)

        self.le_Efficiency.setText('{:0.3f}'.format(Model.Eff * 100))  # Show as percentage
        self.le_PowerStroke.setText('{:0.3f}'.format(Model.air.n * Model.W_Power * CFE))
        self.le_CompressionStroke.setText('{:0.3f}'.format(Model.air.n * Model.W_Compression * CFE))
        self.le_HeatAdded.setText('{:0.3f}'.format(Model.air.n * Model.Q_In * CFE))
        self.lbl_PowerStrokeUnits.setText(Model.units.EUnits)
        self.lbl_CompressionStrokeUnits.setText(Model.units.EUnits)
        self.lbl_HeatInUnits.setText(Model.units.EUnits)

def main():
    oc = ottoCycleController()
    oc.set(T_0=540.0, P_0=1.0, T_High=3600.0, ratio=8.0, V_0=0.02, SI=False)
    oc.plot_cycle_XY(X='v', Y='P', total=True)

if __name__ == "__main__":
    app = qtw.QApplication(sys.argv)
    main()
    sys.exit(app.exec_())
