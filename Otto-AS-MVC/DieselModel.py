# DieselModel.py
from Air import air
import numpy as np

class DieselModel:
    def __init__(self):
        self.r = 18
        self.rc = 2
        self.T1 = 300
        self.P1 = 0.1e6
        self.air = air()
        self.results = {}

    def calculate(self):
        self.air.set(T=self.T1, P=self.P1)
        s1 = self.air.State.s
        u1 = self.air.State.u
        v1 = self.air.State.v

        V1 = v1
        V2 = V1 / self.r
        self.air.set(s=s1, v=V2)
        T2 = self.air.State.T
        P2 = self.air.State.P
        u2 = self.air.State.u

        V3 = self.rc * V2
        self.air.set(P=P2, v=V3)
        T3 = self.air.State.T
        s3 = self.air.State.s
        u3 = self.air.State.u

        V4 = V1
        self.air.set(s=s3, v=V4)
        T4 = self.air.State.T
        P4 = self.air.State.P
        u4 = self.air.State.u

        q_in = u3 - u2
        q_out = u1 - u4
        w_net = q_in - q_out
        eff = 100.0 * w_net / q_in

        self.results = {
            'T1': self.T1, 'T2': T2, 'T3': T3, 'T4': T4,
            'P1': self.P1, 'P2': P2, 'P3': P2, 'P4': P4,
            'q_in': q_in, 'q_out': q_out, 'w_net': w_net, 'eff': eff,
            'pv': self._pv_points(V1, V2, V3, V4, s1, s3, P2)
        }

    def _pv_points(self, V1, V2, V3, V4, s1, s3, P2):
        pv = []
        v = np.linspace(V1, V2, 30)
        for vol in v:
            self.air.set(s=s1, v=vol)
            pv.append((vol, self.air.State.P))
        v = np.linspace(V2, V3, 30)
        for vol in v:
            self.air.set(P=P2, v=vol)
            pv.append((vol, self.air.State.P))
        v = np.linspace(V3, V4, 30)
        for vol in v:
            self.air.set(s=s3, v=vol)
            pv.append((vol, self.air.State.P))
        return pv