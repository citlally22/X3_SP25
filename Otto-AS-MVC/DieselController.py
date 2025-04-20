# DieselController.py
from DieselModel import DieselModel
from DieselView import DieselView

class DieselController:
    def __init__(self):
        self.model = DieselModel()
        self.view = DieselView()

    def setWidgets(self, widget_dict):
        self.view.setWidgets(widget_dict)

    def calculate(self, inputs):
        self.model.r = float(inputs['CR'])
        self.model.rc = float(inputs['rc'])
        self.model.T1 = float(inputs['T_low'])
        self.model.P1 = float(inputs['P0']) * 1e5  # from bar to Pa
        self.model.calculate()
        self.view.updateResults(self.model.results)
        self.view.plotPV(self.model.results['pv'])