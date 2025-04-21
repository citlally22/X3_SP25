import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve, minimize_scalar
from copy import deepcopy as dc


class StateDataForPlotting:
    """Stores thermodynamic state data for plotting.

    Attributes:
        T (list): Temperatures in Kelvin.
        P (list): Pressures in Pascals.
        h (list): Enthalpies in J/mol.
        u (list): Internal energies in J/mol.
        s (list): Entropies in J/mol*K.
        v (list): Specific volumes in m^3/mol.
    """

    def __init__(self):
        """Initializes empty lists for thermodynamic properties."""
        self.T = []
        self.P = []
        self.h = []
        self.u = []
        self.s = []
        self.v = []

    def clear(self):
        """Clears all stored data lists."""
        self.T.clear()
        self.P.clear()
        self.h.clear()
        self.u.clear()
        self.s.clear()
        self.v.clear()

    def add(self, vals):
        """Adds a set of thermodynamic state values to the respective lists.

        Args:
            vals (tuple): Tuple containing (T, P, u, h, s, v) values.
        """
        T, P, u, h, s, v = vals
        self.T.append(T)
        self.P.append(P)
        self.h.append(h)
        self.u.append(u)
        self.s.append(s)
        self.v.append(v)

    def getAxisLabel(self, W='T', Units=None):
        """Returns the plot axis label for a specified property.

        Args:
            W (str, optional): Property identifier ('T', 'P', 'h', 'u', 's', 'v'). Defaults to 'T'.
            Units (units, optional): Units object for formatting labels. Defaults to None.

        Returns:
            str: Formatted axis label with units.
        """
        Units = Units if Units is not None else units()
        w = W.lower()
        if w == 't':
            return Units.TPlotUnits
        if w == 'h':
            return Units.hPlotUnits
        if w == 'u':
            return Units.uPlotUnits
        if w == 's':
            return Units.sPlotUnits
        if w == 'v':
            return Units.vPlotUnits
        if w == 'p':
            return Units.PPlotUnits
        return ""  # Return empty string for invalid input

    def getDataCol(self, W='T'):
        """Retrieves the data list for a specified property.

        Args:
            W (str, optional): Property identifier ('T', 'P', 'h', 'u', 's', 'v'). Defaults to 'T'.

        Returns:
            list: List of values for the specified property.
        """
        w = W.lower()
        if w == 't':
            return self.T
        if w == 'h':
            return self.h
        if w == 'u':
            return self.u
        if w == 's':
            return self.s
        if w == 'v':
            return self.v
        if w == 'p':
            return self.P
        return []  # Return empty list for invalid input


class stateProps:
    """Stores and manages thermodynamic state properties.

    Attributes:
        name (str): Optional name for the state.
        T (float): Temperature in Kelvin.
        P (float): Pressure in Pascals.
        h (float): Enthalpy in J/mol.
        u (float): Internal energy in J/mol.
        s (float): Entropy in J/mol*K.
        v (float): Specific volume in m^3/mol.
    """

    def __init__(self):
        """Initializes an empty thermodynamic state."""
        self.name = None
        self.T = None
        self.P = None
        self.h = None
        self.u = None
        self.s = None
        self.v = None

    def __mul__(self, other):
        """Multiplies thermodynamic properties by a scalar.

        Args:
            other (float or int): Scalar to multiply by.

        Returns:
            stateProps: New state with scaled properties (h, u, s, v).
        """
        if isinstance(other, (float, int)):
            b = stateProps()
            b.T = self.T
            b.P = self.P
            b.h = self.h * other if self.h is not None else None
            b.u = self.u * other if self.u is not None else None
            b.s = self.s * other if self.s is not None else None
            b.v = self.v * other if self.v is not None else None
            return b
        return NotImplemented

    def __rmul__(self, other):
        """Handles right multiplication (scalar * stateProps).

        Args:
            other (float or int): Scalar to multiply by.

        Returns:
            stateProps: Result of multiplication.
        """
        return self * other

    def __truediv__(self, other):
        """Divides thermodynamic properties by a scalar.

        Args:
            other (float or int): Scalar to divide by.

        Returns:
            stateProps: New state with scaled properties (h, u, s, v).
        """
        if isinstance(other, (float, int)):
            b = stateProps()
            b.T = self.T
            b.P = self.P
            b.h = self.h / other if self.h is not None else None
            b.u = self.u / other if self.u is not None else None
            b.s = self.s / other if self.s is not None else None
            b.v = self.v / other if self.v is not None else None
            return b
        return NotImplemented

    def ConvertStateData(self, SI=True, mass=False, total=False, n=1.0, MW=1.0, Units=None):
        """Converts state data to specified units (SI or US, molar/mass/total basis).

        Args:
            SI (bool, optional): Use SI units if True, else US customary. Defaults to True.
            mass (bool, optional): Use mass basis if True. Defaults to False.
            total (bool, optional): Use total properties if True. Defaults to False.
            n (float, optional): Number of moles. Defaults to 1.0.
            MW (float, optional): Molecular weight in kg/kmol. Defaults to 1.0.
            Units (units, optional): Units object for conversion factors. Defaults to None.
        """
        UC = Units if Units is not None else units()
        UC.set(SI=SI, mass=mass, total=total)
        mCF = 1.0 if SI else UC.CF_Mass
        TCF = 1.0 if SI else UC.CF_T
        PCF = 1.0 if SI else UC.CF_P
        vCF = 1.0 if SI else UC.CF_v
        uCF = 1.0 if SI else UC.CF_e
        hCF = 1.0 if SI else UC.CF_e
        sCF = 1.0 if SI else UC.CF_s
        nCF = 1.0 if SI else UC.CF_n
        if mass:
            mCF /= MW
            vCF /= MW
            uCF /= MW
            hCF /= MW
            sCF /= MW
        elif total:
            vCF *= n * nCF
            uCF *= n * nCF
            hCF *= n * nCF
            sCF *= n * nCF

        # Apply conversion factors to non-None properties
        self.P = self.P * PCF if self.P is not None else None
        self.T = self.T * TCF if self.T is not None else None
        self.h = self.h * hCF if self.h is not None else None
        self.u = self.u * uCF if self.u is not None else None
        self.v = self.v * vCF if self.v is not None else None
        self.s = self.s * sCF if self.s is not None else None

    def getVal(self, name='T'):
        """Retrieves the value of a specified thermodynamic property.

        Args:
            name (str, optional): Property identifier ('T', 'P', 'h', 'u', 's', 'v'). Defaults to 'T'.

        Returns:
            float: Value of the specified property or None if not set.
        """
        n = name.lower()
        if n == 't':
            return self.T
        if n == 'h':
            return self.h
        if n == 'u':
            return self.u
        if n == 's':
            return self.s
        if n == 'v':
            return self.v
        if n == 'p':
            return self.P
        return None

    def print(self, units=None):
        """Prints thermodynamic state properties with units.

        Args:
            units (units, optional): Units object for formatting. Defaults to None.
        """
        U = units if units is not None else units()
        if self.name is not None:
            print(self.name)
        print(f'v={self.v:.4f} {U.vUnits}.')
        print(f'u={self.u:.4f} {U.uUnits}')
        print(f'h={self.h:.4f} {U.hUnits}')
        print(f's={self.s:.4f} {U.sUnits}')


class units:
    """Manages unit conversions and labels for thermodynamic properties.

    Attributes:
        SI (bool): True for SI units, False for US customary units.
        sUnits (str): Entropy units.
        uUnits (str): Internal energy units.
        vUnits (str): Specific volume units.
        VUnits (str): Total volume units.
        hUnits (str): Enthalpy units.
        mUnits (str): Mass units.
        TUnits (str): Temperature units.
        PUnits (str): Pressure units.
        EUnits (str): Energy units.
        CF_* (float): Conversion factors for various properties.
    """

    def __init__(self):
        """Initializes unit settings with default SI units on a molar basis."""
        self.SI = True
        self.sUnits = 'J/mol*K'
        self.uUnits = 'J/mol'
        self.vUnits = 'm^3/mol'
        self.VUnits = 'm^3'
        self.hUnits = self.uUnits
        self.mUnits = 'kg'
        self.TUnits = 'K'
        self.PUnits = 'Pa'
        self.EUnits = 'J'

        # Define conversion factors from SI to US customary units
        self.CF_E = 1.0 / 1055.06  # J to Btu
        self.CF_Length = 3.28084  # m to ft
        self.CF_V = self.CF_Length ** 3.0  # m^3 to ft^3
        self.CF_P = 1.0 / 101325  # Pa to atm
        self.CF_Mass = 2.20462  # kg to lb
        self.CF_T = 9.0 / 5.0  # K to R
        self.CF_n = 1 / 453.59  # mol to lbmol
        self.CF_v = self.CF_V / self.CF_n  # m^3/mol to ft^3/lbmol
        self.CF_e = self.CF_E / self.CF_n  # J/mol to Btu/lbmol
        self.CF_s = self.CF_e / (self.CF_n * self.CF_T)  # J/mol*K to Btu/lbmol*R

        self.setPlotUnits()

    def set(self, SI=True, mass=False, total=False):
        """Configures units based on system, basis, and extent.

        Args:
            SI (bool, optional): Use SI units if True, else US customary. Defaults to True.
            mass (bool, optional): Use mass basis if True, else molar. Defaults to False.
            total (bool, optional): Use total properties if True. Defaults to False.
        """
        self.SI = SI
        if SI:
            self.sUnits = f'J/{"" if total else ("kg*" if mass else "mol*")}K'
            self.uUnits = f'J{"" if total else ("/kg" if mass else "/mol")}'
            self.vUnits = f'm^3{"" if total else ("/kg" if mass else "/mol")}'
            self.VUnits = 'm^3'
            self.hUnits = self.uUnits
            self.mUnits = 'kg'
            self.TUnits = 'K'
            self.PUnits = 'Pa'
            self.EUnits = 'J'
        else:
            self.sUnits = f'BTU/{"" if total else ("lb*" if mass else "lbmol*")}R'
            self.uUnits = f'BTU{"" if total else ("/lb" if mass else "/lbmol")}'
            self.vUnits = f'ft^3{"" if total else ("/lb" if mass else "/lbmol")}'
            self.VUnits = 'ft^3'
            self.hUnits = self.uUnits
            self.mUnits = 'lb'
            self.TUnits = 'R'
            self.PUnits = 'Atm'
            self.EUnits = 'Btu'

        self.setPlotUnits(SI=SI, mass=mass, total=total)

    def setPlotUnits(self, SI=True, mass=True, total=False):
        """Sets plot unit labels for thermodynamic properties.

        Args:
            SI (bool, optional): Use SI units if True. Defaults to True.
            mass (bool, optional): Use mass basis if True. Defaults to True.
            total (bool, optional): Use total properties if True. Defaults to False.
        """
        if SI:
            self.PPlotUnits = r'P $\left(Pa\right)$'
            self.TPlotUnits = r'T $\left(K\right)$'
            if total:
                self.sPlotUnits = r'S $\left(\frac{J}{K}\right)$'
                self.uPlotUnits = r'U $\left(J\right)$'
                self.hPlotUnits = r'H $\left(J\right)$'
                self.vPlotUnits = r'V $\left(m^3\right)$'
            elif mass:
                self.sPlotUnits = r's $\left(\frac{J}{kg*K}\right)$'
                self.uPlotUnits = r'u $\left(\frac{J}{kg}\right)$'
                self.hPlotUnits = r'h $\left(\frac{J}{kg}\right)$'
                self.vPlotUnits = r'v $\left(\frac{m^3}{kg}\right)$'
            else:
                self.sPlotUnits = r'$\bar{s} \left(\frac{J}{mol*K}\right)$'
                self.uPlotUnits = r'$\bar{u} \left(\frac{J}{mol}\right)$'
                self.hPlotUnits = r'$\bar{h} \left(\frac{J}{mol}\right)$'
                self.vPlotUnits = r'$\bar{v} \left(\frac{m^3}{mol}\right)$'
        else:
            self.PPlotUnits = r'P $\left(atm\right)$'
            self.TPlotUnits = r'T $\left(^{o}R\right)$'
            if total:
                self.sPlotUnits = r'S $\left(\frac{Btu}{^{o}R}\right)$'
                self.uPlotUnits = r'U $\left(Btu\right)$'
                self.hPlotUnits = r'H $\left(Btu\right)$'
                self.vPlotUnits = r'V $\left(ft^3\right)$'
            elif mass:
                self.sPlotUnits = r's $\left(\frac{Btu}{lb\cdot^{o}R}\right)$'
                self.uPlotUnits = r'u $\left(\frac{Btu}{lb}\right)$'
                self.hPlotUnits = r'h $\left(\frac{Btu}{lb}\right)$'
                self.vPlotUnits = r'v $\left(\frac{ft^3}{lb}\right)$'
            else:
                self.sPlotUnits = r'$\bar{s} \left(\frac{Btu}{lb_{mol}\cdot^{o}R}\right)$'
                self.uPlotUnits = r'$\bar{u} \left(\frac{Btu}{lb_{mol}}\right)$'
                self.hPlotUnits = r'$\bar{h} \left(\frac{Btu}{lb_{mol}}\right)$'
                self.vPlotUnits = r'$\bar{v} \left(\frac{ft^3}{lb_{mol}}\right)$'

    def T_RtoK(self, T):
        """Converts temperature from Rankine to Kelvin.

        Args:
            T (float): Temperature in Rankine.

        Returns:
            float: Temperature in Kelvin.
        """
        return T * 5.0 / 9.0

    def T_FtoC(self, T):
        """Converts temperature from Fahrenheit to Celsius.

        Args:
            T (float): Temperature in Fahrenheit.

        Returns:
            float: Temperature in Celsius.
        """
        return (T - 32.0) * 5.0 / 9.0

    def T_RtoF(self, T):
        """Converts temperature from Rankine to Fahrenheit.

        Args:
            T (float): Temperature in Rankine.

        Returns:
            float: Temperature in Fahrenheit.
        """
        return T - 459.67

    def T_FtoK(self, T):
        """Converts temperature from Fahrenheit to Kelvin.

        Args:
            T (float): Temperature in Fahrenheit.

        Returns:
            float: Temperature in Kelvin.
        """
        return self.T_RtoK(self.T_FtoR(T))

    def T_CtoK(self, T):
        """Converts temperature from Celsius to Kelvin.

        Args:
            T (float): Temperature in Celsius.

        Returns:
            float: Temperature in Kelvin.
        """
        return T + 273.15

    def T_CtoF(self, T):
        """Converts temperature from Celsius to Fahrenheit.

        Args:
            T (float): Temperature in Celsius.

        Returns:
            float: Temperature in Fahrenheit.
        """
        return T * 9.0 / 5.0 + 32

    def T_KtoC(self, T):
        """Converts temperature from Kelvin to Celsius.

        Args:
            T (float): Temperature in Kelvin.

        Returns:
            float: Temperature in Celsius.
        """
        return T - 273.15

    def T_KtoR(self, T):
        """Converts temperature from Kelvin to Rankine.

        Args:
            T (float): Temperature in Kelvin.

        Returns:
            float: Temperature in Rankine.
        """
        return T * 9 / 5

    def T_FtoR(self, T):
        """Converts temperature from Fahrenheit to Rankine.

        Args:
            T (float): Temperature in Fahrenheit.

        Returns:
            float: Temperature in Rankine.
        """
        return T + 459.67


class air:
    """Models air as an ideal gas with thermodynamic properties.

    Attributes:
        RBar (float): Universal gas constant in J/mol*K.
        MW (float): Molecular weight of air in kg/kmol.
        R (float): Gas constant for air in kJ/kg*K.
        StandardState (stateProps): Reference state at T=0°C, P=1 atm.
        State (stateProps): Current thermodynamic state.
        n (float): Number of moles.
        m (float): Mass in kg.
    """

    def __init__(self):
        """Initializes air model with standard state at T=0°C, P=1 atm."""
        self.RBar = 8.3145  # Universal gas constant in J/mol*K
        self.MW = 28.97  # Molecular weight in kg/kmol
        self.R = self.RBar / self.MW  # Specific gas constant in kJ/kg*K
        self.StandardState = stateProps()
        self.StandardState.P = 101325.0  # Standard pressure in Pa
        self.StandardState.T = 273.15  # Standard temperature in K
        self.StandardState.v = self.RBar * self.StandardState.T / self.StandardState.P
        self.StandardState.u = 0  # Reference internal energy
        self.StandardState.h = 0  # Reference enthalpy
        self.StandardState.s = 0  # Reference entropy
        self.State = stateProps()
        self.n = 1.0  # Default number of moles
        self.m = self.n * self.MW / 1000.0  # Mass in kg

    def cv(self, T):
        """Calculates molar specific heat at constant volume.

        Args:
            T (float): Temperature in Kelvin.

        Returns:
            float: Specific heat at constant volume in J/mol*K.
        """
        return self.cp(T) - self.RBar

    def cp(self, T):
        """Calculates molar specific heat at constant pressure using a polynomial.

        Args:
            T (float): Temperature in Kelvin.

        Returns:
            float: Specific heat at constant pressure in J/mol*K.
        """
        TLowRange = 1630.0  # Threshold for low-temperature coefficients
        # Polynomial coefficients for cp = Rbar * (a + b*T + c*T^2 + d*T^3 + e*T^4)
        a = 3.653 if T < TLowRange else 2.753
        b = -1.337e-3 if T < TLowRange else 0.002
        c = 3.294e-6 if T < TLowRange else -1.0e-6
        d = -1.913e-9 if T < TLowRange else 3.0e-10
        e = 0.2763e-12 if T < TLowRange else -3.0e-14
        return self.RBar * (a + b * T + c * T**2 + d * T**3 + e * T**4)

    def deltau(self, T1=None, T2=None):
        """Calculates change in molar internal energy between two temperatures.

        Args:
            T1 (float, optional): Initial temperature in Kelvin. Defaults to standard state.
            T2 (float, optional): Final temperature in Kelvin. Defaults to standard state.

        Returns:
            float: Change in internal energy in J/mol.

        Raises:
            ValueError: If temperatures are non-positive.
        """
        T1 = self.StandardState.T if T1 is None else T1
        T2 = self.StandardState.T if T2 is None else T2
        if T1 <= 0 or T2 <= 0:
            raise ValueError(f"Invalid temperatures in deltau: T1={T1}, T2={T2}")
        return quad(self.cv, T1, T2)[0]  # Integrate cv over temperature range

    def deltah(self, T1=None, T2=None):
        """Calculates change in molar enthalpy between two temperatures.

        Args:
            T1 (float, optional): Initial temperature in Kelvin. Defaults to standard state.
            T2 (float, optional): Final temperature in Kelvin. Defaults to standard state.

        Returns:
            float: Change in enthalpy in J/mol.

        Raises:
            ValueError: If temperatures are non-positive.
        """
        T1 = self.StandardState.T if T1 is None else T1
        T2 = self.StandardState.T if T2 is None else T2
        if T1 <= 0 or T2 <= 0:
            raise ValueError(f"Invalid temperatures in deltah: T1={T1}, T2={T2}")
        return quad(self.cp, T1, T2)[0]  # Integrate cp over temperature range

    def deltas_tv(self, T1=None, T2=None, V1=None, V2=None):
        """Calculates change in molar entropy as a function of temperature and volume.

        Args:
            T1 (float, optional): Initial temperature in Kelvin. Defaults to standard state.
            T2 (float, optional): Final temperature in Kelvin. Defaults to standard state.
            V1 (float, optional): Initial specific volume in m^3/mol. Defaults to standard state.
            V2 (float, optional): Final specific volume in m^3/mol. Defaults to standard state.

        Returns:
            float: Change in entropy in J/mol*K.

        Raises:
            ValueError: If temperatures or volumes are non-positive.
            Exception: If integration or logarithm calculation fails.
        """
        T1 = self.StandardState.T if T1 is None else T1
        T2 = self.StandardState.T if T2 is None else T2
        V1 = self.StandardState.v if V1 is None else V1
        V2 = self.StandardState.v if V2 is None else V2
        if T1 <= 0 or T2 <= 0:
            raise ValueError(f"Invalid temperatures in deltas_tv: T1={T1}, T2={T2}")
        if V1 <= 0 or V2 <= 0:
            raise ValueError(f"Invalid volumes in deltas_tv: V1={V1}, V2={V2}")
        fn = lambda T: 0 if T == 0 else self.cv(T) / T
        try:
            deltaS = quad(fn, T1, T2, limit=100)[0]
        except Exception as e:
            raise Exception(f"Integration failed in deltas_tv: T1={T1}, T2={T2}, error={str(e)}")
        try:
            deltaS += self.RBar * math.log(V2 / V1)  # Volume contribution to entropy
        except Exception as e:
            raise Exception(f"Log calculation failed in deltas_tv: V1={V1}, V2={V2}, error={str(e)}")
        return deltaS

    def deltas_tp(self, T1=None, T2=None, P1=None, P2=None):
        """Calculates change in molar entropy as a function of temperature and pressure.

        Args:
            T1 (float, optional): Initial temperature in Kelvin. Defaults to standard state.
            T2 (float, optional): Final temperature in Kelvin. Defaults to standard state.
            P1 (float, optional): Initial pressure in Pa. Defaults to standard state.
            P2 (float, optional): Final pressure in Pa. Defaults to standard state.

        Returns:
            float: Change in entropy in J/mol*K.

        Raises:
            ValueError: If temperatures or pressures are non-positive.
            Exception: If integration or logarithm calculation fails.
        """
        T1 = self.StandardState.T if T1 is None else T1
        T2 = self.StandardState.T if T2 is None else T2
        P1 = self.StandardState.P if P1 is None else P1
        P2 = self.StandardState.P if P2 is None else P2
        if T1 <= 0 or T2 <= 0:
            raise ValueError(f"Invalid temperatures in deltas_tp: T1={T1}, T2={T2}")
        if P1 <= 0 or P2 <= 0:
            raise ValueError(f"Invalid pressures in deltas_tp: P1={P1}, P2={P2}")
        fn = lambda T: 0 if T == 0.0 else self.cp(T) / T
        try:
            deltaS = quad(fn, T1, T2, limit=100)[0]
        except Exception as e:
            raise Exception(f"Integration failed in deltas_tp: T1={T1}, T2={T2}, error={str(e)}")
        try:
            deltaS += self.RBar * math.log(P1 / P2)  # Pressure contribution to entropy
        except Exception as e:
            raise Exception(f"Log calculation failed in deltas_tp: P1={P1}, P2={P2}, error={str(e)}")
        return deltaS

    def set(self, P=None, T=None, v=None, h=None, u=None, s=None, name=None):
        """Sets two thermodynamic properties and calculates the full state.

        Args:
            P (float, optional): Pressure in Pa. Defaults to None.
            T (float, optional): Temperature in Kelvin. Defaults to None.
            v (float, optional): Specific volume in m^3/mol. Defaults to None.
            h (float, optional): Enthalpy in J/mol. Defaults to None.
            u (float, optional): Internal energy in J/mol. Defaults to None.
            s (float, optional): Entropy in J/mol*K. Defaults to None.
            name (str, optional): Name of the state. Defaults to None.

        Returns:
            stateProps: Deep copy of the calculated state.
        """
        # Assign provided properties to the state
        self.State.P = P
        self.State.T = T
        self.State.v = v
        self.State.h = h
        self.State.u = u
        self.State.s = s
        self.State.name = name
        if all(x is None for x in [T, P, u, v, h, s]):
            return  # No properties provided, exit
        self.calc()  # Calculate remaining properties
        return dc(self.State)

    def calc(self):
        """Calculates the complete thermodynamic state based on two known properties.

        Uses the ideal gas law and specific heat functions relative to the standard
        state (T=0°C, P=101.325 kPa, u=0, h=0, s=0).

        Raises:
            ValueError: If input properties are invalid (e.g., non-positive).
            Exception: If numerical solver fails to converge.
        """
        # Define physical bounds for temperature and pressure
        T_MIN, T_MAX = 1e-6, 5000  # K
        P_MIN, P_MAX = 1e-6, 1e8   # Pa

        # Case 1: P, T - Direct calculation using ideal gas law
        if self.State.P is not None and self.State.T is not None:
            if self.State.P <= 0 or self.State.T <= 0:
                raise ValueError(f"Invalid P,T: P={self.State.P}, T={self.State.T}")
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.u = self.deltau(T2=self.State.T)
            self.State.h = self.deltah(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 2: P, u - Solve for T numerically
        elif self.State.P is not None and self.State.u is not None:
            if self.State.P <= 0:
                raise ValueError(f"Invalid P in P,u: P={self.State.P}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltau(T2=T) - self.State.u) ** 2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in P,u case: u={self.State.u}, P={self.State.P}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.h = self.deltah(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 3: P, v - Calculate T using ideal gas law
        elif self.State.P is not None and self.State.v is not None:
            if self.State.P <= 0 or self.State.v <= 0:
                raise ValueError(f"Invalid P,v: P={self.State.P}, v={self.State.v}")
            self.State.T = self.State.v * self.State.P / self.RBar
            self.State.u = self.deltau(T2=self.State.T)
            self.State.h = self.deltah(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 4: P, h - Solve for T numerically
        elif self.State.P is not None and self.State.h is not None:
            if self.State.P <= 0:
                raise ValueError(f"Invalid P in P,h: P={self.State.P}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltah(T2=T) - self.State.h) ** 2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in P,h case: h={self.State.h}, P={self.State.P}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.u = self.deltau(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 5: P, s - Solve for T numerically
        elif self.State.P is not None and self.State.s is not None:
            if self.State.P <= 0:
                raise ValueError(f"Invalid P in P,s: P={self.State.P}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltas_tp(T2=T, P2=self.State.P) - self.State.s) ** 2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in P,s case: s={self.State.s}, P={self.State.P}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.u = self.deltau(T2=self.State.T)
            self.State.h = self.deltah(T2=self.State.T)
        # Case 7: T, v - Calculate P using ideal gas law
        elif self.State.T is not None and self.State.v is not None:
            if self.State.T <= 0 or self.State.v <= 0:
                raise ValueError(f"Invalid T,v: T={self.State.T}, v={self.State.v}")
            self.State.P = self.State.T * self.RBar / self.State.v
            self.State.u = self.deltau(T2=self.State.T)
            self.State.h = self.deltah(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 9: T, s - Solve for P numerically
        elif self.State.T is not None and self.State.s is not None:
            if self.State.T <= 0:
                raise ValueError(f"Invalid T in T,s: T={self.State.T}")
            def objective(P):
                P = max(P, P_MIN)
                return (self.deltas_tp(T2=self.State.T, P2=P) - self.State.s) ** 2
            result = minimize_scalar(objective, bounds=(P_MIN, P_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for P in T,s case: s={self.State.s}, T={self.State.T}, error={result.message}")
            self.State.P = max(result.x, P_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.u = self.deltau(T2=self.State.T)
            self.State.h = self.deltah(T2=self.State.T)
        # Case 10: u, v - Solve for T numerically
        elif self.State.u is not None and self.State.v is not None:
            if self.State.v <= 0:
                raise ValueError(f"Invalid v in u,v: v={self.State.v}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltau(T2=T) - self.State.u) ** 2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in u,v case: u={self.State.u}, v={self.State.v}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.P = self.State.T * self.RBar / self.State.v
            self.State.h = self.deltah(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 12: u, s - Solve for T and P numerically
        elif self.State.u is not None and self.State.s is not None:
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltau(T2=T) - self.State.u) ** 2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in u,s case: u={self.State.u}, s={self.State.s}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            def objective(P):
                P = max(P, P_MIN)
                return (self.deltas_tp(T2=self.State.T, P2=P) - self.State.s) ** 2
            result = minimize_scalar(objective, bounds=(P_MIN, P_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for P in u,s case: s={self.State.s}, T={self.State.T}, error={result.message}")
            self.State.P = max(result.x, P_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.h = self.deltah(T2=self.State.T)
        # Case 13: v, h - Solve for T numerically
        elif self.State.v is not None and self.State.h is not None:
            if self.State.v <= 0:
                raise ValueError(f"Invalid v in v,h: v={self.State.v}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltah(T2=T) - self.State.h) ** 2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in v,h case: h={self.State.h}, v={self.State.v}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.P = self.State.T * self.RBar / self.State.v
            self.State.u = self.deltau(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 14: v, s - Solve for T numerically
        elif self.State.v is not None and self.State.s is not None:
            if self.State.v <= 0:
                raise ValueError(f"Invalid v in v,s: v={self.State.v}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltas_tv(T2=T, V2=self.State.v) - self.State.s) ** 2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in v,s case: s={self.State.s}, v={self.State.v}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.P = self.RBar * self.State.T / self.State.v
            self.State.h = self.deltah(T2=self.State.T)
            self.State.u = self.deltau(T2=self.State.T)
        # Case 15: h, s - Solve for T and P numerically
        elif self.State.h is not None and self.State.s is not None:
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltah(T2=T) - self.State.h) ** 2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in h,s case: h={self.State.h}, s={self.State.s}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            def objective(P):
                P = max(P, P_MIN)
                return (self.deltas_tp(T2=self.State.T, P2=P) - self.State.s) ** 2
            result = minimize_scalar(objective, bounds=(P_MIN, P_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for P in h,s case: s={self.State.s}, T={self.State.T}, error={result.message}")
            self.State.P = max(result.x, P_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.u = self.deltau(T2=self.State.T)

    def getSummary_MassBasis(self, units=None):
        """Returns a string summary of the state on a mass basis.

        Args:
            units (units, optional): Units object for formatting. Defaults to None.

        Returns:
            str: Formatted string with thermodynamic properties.
        """
        UC = units if units is not None else units()
        mCF = 1.0 if UC.SI else UC.CF_Mass
        TCF = 1.0 if UC.SI else UC.CF_T
        PCF = 1.0 if UC.SI else UC.CF_P
        vCF = 1.0 if UC.SI else UC.CF_V
        uCF = 1.0 if UC.SI else UC.CF_E
        hCF = 1.0 if UC.SI else UC.CF_E
        sCF = 1.0 if UC.SI else UC.CF_s

        # Convert properties to mass basis and format
        stTmp = ''
        stTmp += f'T={self.State.T * TCF:.2f} {UC.TUnits}\n'
        stTmp += f'P={self.State.P * PCF / 1000.0:.3f} {UC.PUnits}\n'
        stTmp += f'v={self.State.v * vCF * 1000.0 / self.MW:.4f} {UC.vUnits}\n'
        stTmp += f'u={self.State.u * uCF / self.MW:.4f} {UC.uUnits}\n'
        stTmp += f'h={self.State.h * hCF / self.MW:.4f} {UC.hUnits}\n'
        stTmp += f's={self.State.s * sCF / self.MW:.4f} {UC.sUnits}'
        return stTmp

    def print_MassBasis(self):
        """Prints the state summary on a mass basis."""
        print(self.getSummary_MassBasis())

    def getSummary_Extensive(self, units=None):
        """Returns a string summary of extensive state properties.

        Args:
            units (units, optional): Units object for formatting. Defaults to None.

        Returns:
            str: Formatted string with extensive properties.
        """
        UC = units if units is not None else units()
        mCF = 1.0 if UC.SI else UC.CF_Mass
        TCF = 1.0 if UC.SI else UC.CF_T
        PCF = 1.0 if UC.SI else UC.CF_P
        vCF = 1.0 if UC.SI else UC.CF_V
        uCF = 1.0 if UC.SI else UC.CF_E
        hCF = 1.0 if UC.SI else UC.CF_E
        sCF = 1.0 if UC.SI else UC.CF_s

        # Scale properties by number of moles
        stTmp = ''
        stTmp += f'T={self.n * self.State.T * TCF:.2f} {UC.TUnits}\n'
        stTmp += f'P={self.n * self.State.P * PCF / 1000.0:.3f} {UC.PUnits}\n'
        stTmp += f'v={self.n * self.State.v * vCF * 1000.0:.4f} {UC.vUnits}\n'
        stTmp += f'u={self.n * self.State.u * uCF:.4f} {UC.uUnits}\n'
        stTmp += f'h={self.n * self.State.h * hCF:.4f} {UC.hUnits}\n'
        stTmp += f's={self.n * self.State.s * sCF:.4f} {UC.sUnits}'
        return stTmp

    def print_Extensive(self):
        """Prints extensive thermodynamic properties."""
        ext = self.State * self.n
        print(f'T={ext.T:.2f} K')
        print(f'P={ext.P / 1000.0:.3f} kPa')
        print(f'v={ext.v:.4f} m^3.')
        print(f'u={ext.u:.4f} kJ')
        print(f'h={ext.h:.4f} kJ')
        print(f's={ext.s:.4f} kJ/K')


def main():
    """Demonstrates usage of the air class by setting and printing a state."""
    a = air()
    a.set(P=a.StandardState.P, T=200)
    a.print_Extensive()


if __name__ == "__main__":
    main()