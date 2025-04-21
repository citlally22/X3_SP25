from Air import air
import numpy as np


class DieselModel:
    """Models a Diesel cycle for thermodynamic analysis.

    Calculates the states, work, heat transfer, and efficiency of a Diesel cycle
    using air as an ideal gas with constant specific heats. Generates P-V diagram
    points for plotting.

    Attributes:
        r (float): Compression ratio (V1/V2).
        rc (float): Cutoff ratio (V3/V2), updated based on T_high.
        T1 (float): Temperature at state 1 in Kelvin.
        T_high (float): Maximum temperature (T3) in Kelvin.
        P1 (float): Pressure at state 1 in Pa.
        V1 (float): Volume at state 1 in m³.
        air (air): Instance of the air class for thermodynamic properties.
        results (dict): Dictionary containing cycle results (temperatures, pressures, etc.).
        R (float): Specific gas constant for air in J/kg·K.
        gamma (float): Specific heat ratio (cp/cv) for air.
        cp (float): Specific heat at constant pressure in J/kg·K.
        cv (float): Specific heat at constant volume in J/kg·K.
    """

    def __init__(self):
        """Initializes the Diesel model with default parameters."""
        self.r = 18                        # Default compression ratio
        self.rc = 2                        # Default cutoff ratio
        self.T1 = 300                      # Initial temperature in Kelvin
        self.T_high = 2000                 # Maximum temperature in Kelvin
        self.P1 = 0.1e6                    # Initial pressure in Pa
        self.V1 = 0.001                    # Initial volume in m³
        self.air = air()                   # Air model for thermodynamic calculations
        self.results = {}                  # Store cycle results
        self.R = 287                       # Gas constant for air in J/kg·K
        self.gamma = 1.4                   # Specific heat ratio
        self.cp = 1005                     # Specific heat at constant pressure in J/kg·K
        self.cv = self.cp / self.gamma     # Specific heat at constant volume in J/kg·K

    def calculate(self):
        """Calculates the Diesel cycle states, work, heat transfer, and efficiency.

        Computes the four states of the Diesel cycle (isentropic compression,
        constant-pressure heat addition, isentropic expansion, constant-volume
        heat rejection) using constant specific heats and ideal gas assumptions.
        Stores results in self.results.
        """
        # State 1: Initial state using ideal gas law
        self.air.set(T=self.T1, P=self.P1)  # Set air state
        s1 = self.air.State.s               # Entropy at state 1 in J/mol·K
        v1 = (self.R * self.T1) / self.P1   # Specific volume in m³/kg
        m = self.V1 / v1                    # Mass of air in kg
        u1 = self.cv * self.T1              # Internal energy in J/kg

        # State 2: Isentropic compression (1-2)
        v2 = v1 / self.r                    # Specific volume after compression
        T2 = self.T1 * (self.r ** (self.gamma - 1))  # Isentropic temperature relation
        self.air.set(T=T2, v=v2)            # Set air state
        P2 = self.R * T2 / v2               # Pressure using ideal gas law
        u2 = self.cv * T2                   # Internal energy in J/kg
        h2 = self.cp * T2                   # Enthalpy in J/kg
        V2 = v2 * m                         # Total volume in m³

        # State 3: Constant-pressure heat addition (2-3), T3 = T_high
        T3 = self.T_high                    # Set maximum temperature
        v3 = (self.R * T3) / P2             # Specific volume at constant pressure
        self.air.set(T=T3, P=P2)            # Set air state
        s3 = self.air.State.s               # Entropy at state 3 in J/mol·K
        u3 = self.cv * T3                   # Internal energy in J/kg
        h3 = self.cp * T3                   # Enthalpy in J/kg
        V3 = v3 * m                         # Total volume in m³
        self.rc = v3 / v2                   # Update cutoff ratio based on T_high

        # State 4: Isentropic expansion (3-4)
        v4 = v1                             # Specific volume returns to v1
        expansion_ratio = v3 / v4           # Should equal r * rc for consistency
        T4 = T3 / (expansion_ratio ** (self.gamma - 1))  # Isentropic temperature relation
        self.air.set(T=T4, v=v4)            # Set air state
        P4 = self.R * T4 / v4               # Pressure using ideal gas law
        u4 = self.cv * T4                   # Internal energy in J/kg
        h4 = self.cp * T4                   # Enthalpy in J/kg
        V4 = v4 * m                         # Total volume in m³

        # Calculate thermodynamic quantities (all per kg)
        q_in = h3 - h2                      # Heat input at constant pressure
        q_out = u4 - u1                     # Heat rejection at constant volume
        w_comp = u2 - u1                    # Compression work (isentropic)
        w_exp = u3 - u4                     # Expansion work (isentropic)
        w_net = w_exp - w_comp              # Net work output
        eff = 100.0 * w_net / q_in if q_in != 0 else 0  # Thermal efficiency in %

        # Store results in dictionary
        self.results = {
            'T1': self.T1, 'T2': T2, 'T3': T3, 'T4': T4,  # Temperatures in K
            'P1': self.P1, 'P2': P2, 'P3': P2, 'P4': P4,  # Pressures in Pa
            'q_in': q_in, 'q_out': q_out, 'w_net': w_net,  # Heat and work in J/kg
            'w_comp': w_comp, 'w_exp': w_exp, 'eff': eff,  # Work and efficiency
            'pv': self._pv_points(v1, v2, v3, v4, s1, s3, P2)  # P-V diagram points
        }

    def _pv_points(self, v1, v2, v3, v4, s1, s3, P2):
        """Generates points for the P-V diagram of the Diesel cycle.

        Creates points along the isentropic compression (1-2), constant-pressure
        heat addition (2-3), and isentropic expansion (3-4) processes.

        Args:
            v1 (float): Specific volume at state 1 in m³/kg.
            v2 (float): Specific volume at state 2 in m³/kg.
            v3 (float): Specific volume at state 3 in m³/kg.
            v4 (float): Specific volume at state 4 in m³/kg.
            s1 (float): Entropy at state 1 in J/mol·K.
            s3 (float): Entropy at state 3 in J/mol·K.
            P2 (float): Pressure at state 2 (and 3) in Pa.

        Returns:
            list: List of (volume, pressure) tuples for plotting.
        """
        pv = []
        # Isentropic compression (1-2)
        v = np.linspace(v1, v2, 30)  # Generate 30 points from v1 to v2
        for vol in v:
            self.air.set(s=s1, v=vol)  # Set state with constant entropy
            pv.append((vol, self.air.State.P))  # Store volume and pressure

        # Constant-pressure heat addition (2-3)
        v = np.linspace(v2, v3, 30)  # Generate 30 points from v2 to v3
        for vol in v:
            self.air.set(P=P2, v=vol)  # Set state with constant pressure
            pv.append((vol, self.air.State.P))  # Store volume and pressure

        # Isentropic expansion (3-4)
        v = np.linspace(v3, v4, 30)  # Generate 30 points from v3 to v4
        for vol in v:
            self.air.set(s=s3, v=vol)  # Set state with constant entropy
            pv.append((vol, self.air.State.P))  # Store volume and pressure

        return pv