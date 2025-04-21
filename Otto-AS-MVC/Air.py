import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve, minimize_scalar
from copy import deepcopy as dc

class StateDataForPlotting:
    """
    I'm making this class for easy storage of data for plotting.
    """
    def __init__(self):
        self.T = []
        self.P = []
        self.h = []
        self.u = []
        self.s = []
        self.v = []

    def clear(self):
        self.T.clear()
        self.P.clear()
        self.h.clear()
        self.u.clear()
        self.s.clear()
        self.v.clear()

    def add(self, vals):
        T, P, u, h, s, v = vals
        self.T.append(T)
        self.P.append(P)
        self.h.append(h)
        self.u.append(u)
        self.s.append(s)
        self.v.append(v)

    def getAxisLabel(self, W='T', Units=None):
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

    def getDataCol(self, W='T'):
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

class stateProps:
    """
    For storage and retrieval of a thermodynamic state
    T, P, u, h, s, v
    """
    def __init__(self):
        self.name = None
        self.T = None
        self.P = None
        self.h = None
        self.u = None
        self.s = None
        self.v = None

    def __mul__(self, other):
        if type(other) in (float, int):
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
        return self * other

    def __truediv__(self, other):
        if type(other) in (float, int):
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

        self.P = self.P * PCF if self.P is not None else None
        self.T = self.T * TCF if self.T is not None else None
        self.h = self.h * hCF if self.h is not None else None
        self.u = self.u * uCF if self.u is not None else None
        self.v = self.v * vCF if self.v is not None else None
        self.s = self.s * sCF if self.s is not None else None

    def getVal(self, name='T'):
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

    def print(self, units=None):
        U = units if units is not None else units()
        if self.name is not None:
            print(self.name)
        print('v={:0.4f} {}.'.format(self.v, U.vUnits))
        print('u={:0.4f} {}'.format(self.u, U.uUnits))
        print('h={:0.4f} {}'.format(self.h, U.hUnits))
        print('s={:0.4f} {}'.format(self.s, U.sUnits))

class units:
    """
    For air, I'm assuming the default units are on a molar basis.
    """
    def __init__(self):
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

        self.CF_E = 1.0 / 1055.06  # J to Btu
        self.CF_Length = 3.28084  # m to ft
        self.CF_V = self.CF_Length**3.0  # m^3 to ft^3
        self.CF_P = 1.0 / 101325  # Pa to atm
        self.CF_Mass = 2.20462  # kg to lb
        self.CF_T = 9.0 / 5.0  # K to R
        self.CF_n = 1 / 453.59  # mol to lbmol
        self.CF_v = self.CF_V / self.CF_n  # m^3/mol to ft^3/lbmol
        self.CF_e = self.CF_E / self.CF_n  # J/mol to Btu/lbmol
        self.CF_s = self.CF_e / (self.CF_n * self.CF_T)  # J/mol*K to Btu/lbmol*R

        self.setPlotUnits()

    def set(self, SI=True, mass=False, total=False):
        self.SI = SI
        if SI:
            self.sUnits = 'J/{}K'.format('' if total else ('kg*' if mass else 'mol*'))
            self.uUnits = 'J{}'.format('' if total else ('/kg' if mass else '/mol'))
            self.vUnits = 'm^3{}'.format('' if total else ('/kg' if mass else '/mol'))
            self.VUnits = 'm^3'
            self.hUnits = self.uUnits
            self.mUnits = 'kg'
            self.TUnits = 'K'
            self.PUnits = 'Pa'
            self.EUnits = 'J'
        else:
            self.sUnits = 'BTU/{}R'.format('' if total else ('lb*' if mass else 'lbmol*'))
            self.uUnits = 'BTU{}'.format('' if total else ('/lb' if mass else '/lbmol'))
            self.vUnits = 'ft^3{}'.format('' if total else ('/lb' if mass else '/lbmol'))
            self.VUnits = 'ft^3'
            self.hUnits = self.uUnits
            self.mUnits = 'lb'
            self.TUnits = 'R'
            self.PUnits = 'Atm'
            self.EUnits = 'Btu'

        self.setPlotUnits(SI=SI, mass=mass, total=total)

    def setPlotUnits(self, SI=True, mass=True, total=False):
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
        return T * 5.0 / 9.0

    def T_FtoC(self, T):
        return (T - 32.0) * 5.0 / 9.0

    def T_RtoF(self, T):
        return T - 459.67

    def T_FtoK(self, T):
        return self.T_RtoK(self.T_FtoR(T))

    def T_CtoK(self, T):
        return T + 273.15

    def T_CtoF(self, T):
        return T * 9.0 / 5.0 + 32

    def T_KtoC(self, T):
        return T - 273.15

    def T_KtoR(self, T):
        return T * 9 / 5

    def T_FtoR(self, T):
        return T + 459.67

class air:
    def __init__(self):
        """
        Air as an ideal gas.
        I choose to always specify air in molar metric units.
        The standard state is T=0C and P=1 atm (101.325kPa)
        So u0=0, h0=0, s0=0, v0=RT/P
        """
        self.RBar = 8.3145  # J/mol*K or kJ/kmol*K
        self.MW = 28.97  # kg/kmol or g/mol or lb/lbmol
        self.R = self.RBar / self.MW  # kJ/kg*K or J/g*K
        self.StandardState = stateProps()
        self.StandardState.P = 101325.0  # P in Pa
        self.StandardState.T = 273.15  # T in K
        self.StandardState.v = self.RBar * self.StandardState.T / self.StandardState.P
        self.StandardState.u = 0
        self.StandardState.h = 0
        self.StandardState.s = 0
        self.State = stateProps()
        self.n = 1.0  # moles
        self.m = self.n * self.MW / 1000.0  # mass in kg

    def cv(self, T):
        return self.cp(T) - self.RBar

    def cp(self, T):
        """
        For air as an ideal gas, cp is a function of temperature as given by:
        cp=Rbar(a+b*T+c*T**2+d*T**3+e*T**4)
        """
        TLowRange = 1630.0
        a = 3.653 if T < TLowRange else 2.753
        b = -1.337E-3 if T < TLowRange else 0.002
        c = 3.294E-6 if T < TLowRange else -1.0E-6
        d = -1.913E-9 if T < TLowRange else 3.0E-10
        e = 0.2763E-12 if T < TLowRange else -3.0E-14
        return self.RBar * (a + b*T + c*T**2 + d*T**3 + e*T**4)

    def deltau(self, T1=None, T2=None):
        """
        To calculate changes in molar internal energy for air as an ideal gas u=u(T)
        cv=du/dT|v -> delta u=int((cv)dT, T1, T2)
        """
        if T1 is None:
            T1 = self.StandardState.T
        if T2 is None:
            T2 = self.StandardState.T
        if T1 <= 0 or T2 <= 0:
            raise ValueError(f"Invalid temperatures in deltau: T1={T1}, T2={T2}")
        return quad(self.cv, T1, T2)[0]

    def deltah(self, T1=None, T2=None):
        """
        To calculate changes in molar internal energy for air as an ideal gas u=u(T)
        cp=dh/dT|p -> delta h=int((cp)dT, T1, T2)
        """
        if T1 is None:
            T1 = self.StandardState.T
        if T2 is None:
            T2 = self.StandardState.T
        if T1 <= 0 or T2 <= 0:
            raise ValueError(f"Invalid temperatures in deltah: T1={T1}, T2={T2}")
        return quad(self.cp, T1, T2)[0]

    def deltas_tv(self, T1=None, T2=None, V1=None, V2=None):
        """
        For calculating changes in molar entropy for air as an ideal gas s=s(T,V)
        Tds=du+Pdv -> delta s = int(cv/T*dT, T1, T2)+R ln(V2/V1)
        """
        if T1 is None:
            T1 = self.StandardState.T
        if T2 is None:
            T2 = self.StandardState.T
        if V1 is None:
            V1 = self.StandardState.v
        if V2 is None:
            V2 = self.StandardState.v
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
            deltaS += self.RBar * math.log(V2 / V1)
        except Exception as e:
            raise Exception(f"Log calculation failed in deltas_tv: V1={V1}, V2={V2}, error={str(e)}")
        return deltaS

    def deltas_tp(self, T1=None, T2=None, P1=None, P2=None):
        """
        For calculating changes in molar entropy for air as an ideal gas s=s(T,V)
        Tds=dh-vdP -> delta s = int(cp/T*dT, T1, T2)-R ln(P2/P1)
        """
        if T1 is None:
            T1 = self.StandardState.T
        if T2 is None:
            T2 = self.StandardState.T
        if P1 is None:
            P1 = self.StandardState.P
        if P2 is None:
            P2 = self.StandardState.P
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
            deltaS += self.RBar * math.log(P1 / P2)
        except Exception as e:
            raise Exception(f"Log calculation failed in deltas_tp: P1={P1}, P2={P2}, error={str(e)}")
        return deltaS

    def set(self, P=None, T=None, v=None, h=None, u=None, s=None, name=None):
        """
        This allows me to set two properties and calculate the state of the air
        """
        self.State.P = P
        self.State.T = T
        self.State.v = v
        self.State.h = h
        self.State.u = u
        self.State.s = s
        self.State.name = name
        if T is None and P is None and u is None and v is None and h is None and s is None:
            return
        else:
            self.calc()
        return dc(self.State)

    def calc(self):
        '''
        To calculate the state of ideal gas air, we use the ideal gas law and specific heat functions relative to
        the standard state of T=0C, P=101.325 kPa where u=0, h=0, s=0, v=vo by declaration
        '''
        # Define physical bounds
        T_MIN, T_MAX = 1e-6, 5000  # K
        P_MIN, P_MAX = 1e-6, 1e8   # Pa

        # Case 1: P, T
        if self.State.P is not None and self.State.T is not None:
            if self.State.P <= 0 or self.State.T <= 0:
                raise ValueError(f"Invalid P,T: P={self.State.P}, T={self.State.T}")
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.u = self.deltau(T2=self.State.T)
            self.State.h = self.deltah(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 2: P, u
        elif self.State.P is not None and self.State.u is not None:
            if self.State.P <= 0:
                raise ValueError(f"Invalid P in P,u: P={self.State.P}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltau(T2=T) - self.State.u)**2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in P,u case: u={self.State.u}, P={self.State.P}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.h = self.deltah(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 3: P, v
        elif self.State.P is not None and self.State.v is not None:
            if self.State.P <= 0 or self.State.v <= 0:
                raise ValueError(f"Invalid P,v: P={self.State.P}, v={self.State.v}")
            self.State.T = self.State.v * self.State.P / self.RBar
            self.State.u = self.deltau(T2=self.State.T)
            self.State.h = self.deltah(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 4: P, h
        elif self.State.P is not None and self.State.h is not None:
            if self.State.P <= 0:
                raise ValueError(f"Invalid P in P,h: P={self.State.P}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltah(T2=T) - self.State.h)**2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in P,h case: h={self.State.h}, P={self.State.P}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.u = self.deltau(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 5: P, s
        elif self.State.P is not None and self.State.s is not None:
            if self.State.P <= 0:
                raise ValueError(f"Invalid P in P,s: P={self.State.P}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltas_tp(T2=T, P2=self.State.P) - self.State.s)**2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in P,s case: s={self.State.s}, P={self.State.P}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.u = self.deltau(T2=self.State.T)
            self.State.h = self.deltah(T2=self.State.T)
        # Case 7: T, v
        elif self.State.T is not None and self.State.v is not None:
            if self.State.T <= 0 or self.State.v <= 0:
                raise ValueError(f"Invalid T,v: T={self.State.T}, v={self.State.v}")
            self.State.P = self.State.T * self.RBar / self.State.v
            self.State.u = self.deltau(T2=self.State.T)
            self.State.h = self.deltah(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 9: T, s
        elif self.State.T is not None and self.State.s is not None:
            if self.State.T <= 0:
                raise ValueError(f"Invalid T in T,s: T={self.State.T}")
            def objective(P):
                P = max(P, P_MIN)
                return (self.deltas_tp(T2=self.State.T, P2=P) - self.State.s)**2
            result = minimize_scalar(objective, bounds=(P_MIN, P_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for P in T,s case: s={self.State.s}, T={self.State.T}, error={result.message}")
            self.State.P = max(result.x, P_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.u = self.deltau(T2=self.State.T)
            self.State.h = self.deltah(T2=self.State.T)
        # Case 10: u, v
        elif self.State.u is not None and self.State.v is not None:
            if self.State.v <= 0:
                raise ValueError(f"Invalid v in u,v: v={self.State.v}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltau(T2=T) - self.State.u)**2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in u,v case: u={self.State.u}, v={self.State.v}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.P = self.State.T * self.RBar / self.State.v
            self.State.h = self.deltah(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 12: u, s
        elif self.State.u is not None and self.State.s is not None:
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltau(T2=T) - self.State.u)**2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in u,s case: u={self.State.u}, s={self.State.s}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            def objective(P):
                P = max(P, P_MIN)
                return (self.deltas_tp(T2=self.State.T, P2=P) - self.State.s)**2
            result = minimize_scalar(objective, bounds=(P_MIN, P_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for P in u,s case: s={self.State.s}, T={self.State.T}, error={result.message}")
            self.State.P = max(result.x, P_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.h = self.deltah(T2=self.State.T)
        # Case 13: v, h
        elif self.State.v is not None and self.State.h is not None:
            if self.State.v <= 0:
                raise ValueError(f"Invalid v in v,h: v={self.State.v}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltah(T2=T) - self.State.h)**2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in v,h case: h={self.State.h}, v={self.State.v}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.P = self.State.T * self.RBar / self.State.v
            self.State.u = self.deltau(T2=self.State.T)
            self.State.s = self.deltas_tp(T2=self.State.T, P2=self.State.P)
        # Case 14: v, s
        elif self.State.v is not None and self.State.s is not None:
            if self.State.v <= 0:
                raise ValueError(f"Invalid v in v,s: v={self.State.v}")
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltas_tv(T2=T, V2=self.State.v) - self.State.s)**2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in v,s case: s={self.State.s}, v={self.State.v}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            self.State.P = self.RBar * self.State.T / self.State.v
            self.State.h = self.deltah(T2=self.State.T)
            self.State.u = self.deltau(T2=self.State.T)
        # Case 15: h, s
        elif self.State.h is not None and self.State.s is not None:
            def objective(T):
                T = max(T, T_MIN)
                return (self.deltah(T2=T) - self.State.h)**2
            result = minimize_scalar(objective, bounds=(T_MIN, T_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for T in h,s case: h={self.State.h}, s={self.State.s}, error={result.message}")
            self.State.T = max(result.x, T_MIN)
            def objective(P):
                P = max(P, P_MIN)
                return (self.deltas_tp(T2=self.State.T, P2=P) - self.State.s)**2
            result = minisize_scalar(objective, bounds=(P_MIN, P_MAX), method='bounded')
            if not result.success:
                raise Exception(f"Failed to solve for P in h,s case: s={self.State.s}, T={self.State.T}, error={result.message}")
            self.State.P = max(result.x, P_MIN)
            self.State.v = self.RBar * self.State.T / self.State.P
            self.State.u = self.deltau(T2=self.State.T)

    def getSummary_MassBasis(self, units=None):
        UC = units if units is not None else units()
        mCF = 1.0 if UC.SI else UC.CF_Mass
        TCF = 1.0 if UC.SI else UC.CF_T
        PCF = 1.0 if UC.SI else UC.CF_P
        vCF = 1.0 if UC.SI else UC.CF_V
        uCF = 1.0 if UC.SI else UC.CF_E
        hCF = 1.0 if UC.SI else UC.CF_E
        sCF = 1.0 if UC.SI else UC.CF_S

        stTmp = ''
        stTmp += 'T={:0.2f} {}\n'.format(self.State.T * TCF, UC.TUnits)
        stTmp += 'P={:0.3f} {}\n'.format(self.State.P * PCF / 1000.0, UC.PUnits)
        stTmp += 'v={:0.4f} {}\n'.format(self.State.v * vCF * 1000.0 / self.MW, UC.vUnits)
        stTmp += 'u={:0.4f} {}\n'.format(self.State.u * uCF / self.MW, UC.uUnits)
        stTmp += 'h={:0.4f} {}\n'.format(self.State.h * hCF / self.MW, UC.hUnits)
        stTmp += 's={:0.4f} {}'.format(self.State.s * sCF / self.MW, UC.sUnits)
        return stTmp

    def print_MassBasis(self):
        print(self.getSummary_MassBasis())

    def getSummary_Extensive(self, units=None):
        UC = units if units is not None else units()
        mCF = 1.0 if UC.SI else UC.CF_Mass
        TCF = 1.0 if UC.SI else UC.CF_T
        PCF = 1.0 if UC.SI else UC.CF_P
        vCF = 1.0 if UC.SI else UC.CF_V
        uCF = 1.0 if UC.SI else UC.CF_E
        hCF = 1.0 if UC.SI else UC.CF_E
        sCF = 1.0 if UC.SI else UC.CF_S

        stTmp = ''
        stTmp += 'T={:0.2f} {}\n'.format(self.n * self.State.T * TCF, UC.TUnits)
        stTmp += 'P={:0.3f} {}\n'.format(self.n * self.State.P * PCF / 1000.0, UC.PUnits)
        stTmp += 'v={:0.4f} {}\n'.format(self.n * self.State.v * vCF * 1000.0, UC.vUnits)
        stTmp += 'u={:0.4f} {}\n'.format(self.n * self.State.u * uCF, UC.uUnits)
        stTmp += 'h={:0.4f} {}\n'.format(self.n * self.State.h * hCF, UC.hUnits)
        stTmp += 's={:0.4f} {}'.format(self.n * self.State.s * sCF, UC.sUnits)
        return stTmp

    def print_Extensive(self):
        ext = self.State * self.n
        print('T={:0.2f} {}'.format(ext.T, 'K'))
        print('P={:0.3f} {}'.format(ext.P / 1000.0, 'kPa'))
        print('v={:0.4f} {}.'.format(ext.v, 'm^3'))
        print('u={:0.4f} {}'.format(ext.u, 'kJ'))
        print('h={:0.4f} {}'.format(ext.h, 'kJ'))
        print('s={:0.4f} {}'.format(ext.s, 'kJ/K'))

def main():
    a = air()
    a.set(P=a.StandardState.P, T=200)
    a.print_Extensive()

if __name__ == "__main__":
    main()