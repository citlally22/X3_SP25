from Air import air
import numpy as np

class DieselModel:
    def __init__(self):
        self.r = 18
        self.rc = 2
        self.T1 = 300
        self.T_high = 2000  # Default T_high in Kelvin
        self.P1 = 0.1e6
        self.V1 = 0.001
        self.air = air()
        self.results = {}
        self.R = 287  # Gas constant for air (J/kg·K)
        self.gamma = 1.4  # Specific heat ratio for air
        self.cp = 1005  # Specific heat at constant pressure (J/kg·K)
        self.cv = self.cp / self.gamma  # Specific heat at constant volume

    def calculate(self):
        # State 1: Use ideal gas law to find specific volume
        self.air.set(T=self.T1, P=self.P1)
        s1 = self.air.State.s
        v1 = (self.R * self.T1) / self.P1
        m = self.V1 / v1  # Mass of the gas
        u1 = self.cv * self.T1  # Internal energy

        v2 = v1 / self.r  # Specific volume at state 2
        # State 2: Isentropic compression (1-2)
        T2 = self.T1 * (self.r ** (self.gamma - 1))
        self.air.set(T=T2, v=v2)
        P2 = self.R * T2 / v2  # Ideal gas law
        u2 = self.cv * T2
        h2 = self.cp * T2
        V2 = v2 * m

        # State 3: Constant-pressure heat addition (2-3), set T3 = T_high
        T3 = self.T_high
        v3 = (self.R * T3) / P2
        self.air.set(T=T3, P=P2)
        s3 = self.air.State.s
        u3 = self.cv * T3
        h3 = self.cp * T3
        V3 = v3 * m
        self.rc = v3 / v2  # Update cutoff ratio based on T_high

        # State 4: Isentropic expansion (3-4)
        v4 = v1
        expansion_ratio = (v3 / v4)  # Should be r * rc
        T4 = T3 * (expansion_ratio) ** (self.gamma - 1)  # Correct isentropic relation
        self.air.set(T=T4, v=v4)
        P4 = self.R * T4 / v4
        u4 = self.cv * T4
        h4 = self.cp * T4
        V4 = v4 * m

        # Calculate heat and work
        q_in = h3 - h2  # Constant pressure: q_in = h3 - h2
        q_out = u4 - u1  # Constant volume: q_out = u4 - u1
        w_comp = u2 - u1  # Compression work (isentropic)
        w_exp = u3 - u4  # Expansion work (isentropic)
        w_net = w_exp - w_comp
        eff = 100.0 * w_net / q_in if q_in != 0 else 0

        self.results = {
            'T1': self.T1, 'T2': T2, 'T3': T3, 'T4': T4,
            'P1': self.P1, 'P2': P2, 'P3': P2, 'P4': P4,
            'q_in': q_in, 'q_out': q_out, 'w_net': w_net, 'w_comp': w_comp, 'w_exp': w_exp, 'eff': eff,
            'pv': self._pv_points(v1, v2, v3, v4, s1, s3, P2)
        }

    def _pv_points(self, v1, v2, v3, v4, s1, s3, P2):
        pv = []
        v = np.linspace(v1, v2, 30)
        for vol in v:
            self.air.set(s=s1, v=vol)
            pv.append((vol, self.air.State.P))
        v = np.linspace(v2, v3, 30)
        for vol in v:
            self.air.set(P=P2, v=vol)
            pv.append((vol, self.air.State.P))
        v = np.linspace(v3, v4, 30)
        for vol in v:
            self.air.set(s=s3, v=vol)
            pv.append((vol, self.air.State.P))
        return pv