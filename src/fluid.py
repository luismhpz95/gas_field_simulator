import numpy as np
import pandas as pd
from src.interpolator import LinearInterpolator

class Fluid:
    
    Pstd = 1.0       # стандартное давление [атм]
    Tstd = 293.15    # стандартная температура [K]
    R = 8.314        # универсальная газовая постоянная [Дж/(моль·К)]
    
    def __init__(self, M: float, rho_c: float, xa: float, xy: float, T: float):
        """
        Параметры
        ----------
        M : float
            Молярная масса [кг/моль]
        rho_c : float
            Плотность в стандартных условиях [кг/м³]
        xa : float
            Молярная доля N₂
        xy : float
            Молярная доля CO₂
        """
        self.M = M
        self.rho_c = rho_c
        self.xa = xa
        self.xy = xy
        self.T = T

    def z(self, P: float) -> float:
        """
        Рассчитать коэффициент сверхсжимаемости Z по методике GERG-91 мод.

        Параметры
        ----------
        P : float
            Давление, атм.
        T : float
            Температура, К.
        """
        xe = 1 -  self.xa -  self.xy    # xe = Xэ

        B_as = 0.72 + 1.875e-5 * (320 - self.T)**2       #B_as = B*

        z_c = 1 - (0.0741 *  self.rho_c - 0.006 - 0.063 *  self.xa - 0.0575 *  self.xy)**2

        Me = (24.05525 * z_c *  self.rho_c - 28.0135 *  self.xa - 44.01 *  self.xy) / xe    # Me = Mэ

        H = 128.64 + 47.479 * Me
        
        B1 = (-0.425468 + 2.865e-3 * self.T - 4.62073e-6 * self.T**2
              + (8.77118e-4 - 5.56281e-6 * self.T + 8.8151e-9 * self.T**2) * H 
              + (-8.24747e-7 + 4.31436e-9 * self.T - 6.08319e-12 * self.T**2) * H**2)

        B2 = -0.1446 + 7.4091e-4 * self.T - 9.1195e-7 * self.T**2

        B23 = -0.339693 + 1.61176e-3 * self.T - 2.04429e-6 * self.T**2

        B3 = -0.86834 + 4.0376e-3 * self.T - 5.1657e-6 * self.T**2
        
        Bm = (xe**2 * B1 + xe *  self.xa * B_as * (B1 + B2) 
              - 1.73 * xe *  self.xy * (B1 * B3)**0.5 
              +  self.xa**2 * B2 + 2 *  self.xa *  self.xy  * B23 
              +  self.xy**2 * B3)

        C1 = (-0.302488 + 1.95861e-3 * self.T - 3.16302e-6 * self.T**2 
              + (6.46422e-4 - 4.22876e-6 * self.T + 6.88157e-9 * self.T**2) * H 
              + (-3.32805e-7 + 2.2316e-9 * self.T - 3.67713e-12 * self.T**2) * H**2)
        
        C2 = 7.8498e-3 - 3.9895e-5 * self.T + 6.1187e-8 * self.T**2

        C3 = 2.0513e-3 + 3.4888e-5 * self.T - 8.3703e-8 * self.T**2

        C223 = 5.52066e-3 - 1.68609e-5 * self.T + 1.57169e-8 * self.T**2

        C233 = 3.58783e-3 + 8.06674e-6 * self.T - 3.25798e-8 * self.T**2

        C_as = 0.92 + 0.0013 * (self.T - 270)    # C_as = C*
        
        Cm = (xe**3 * C1 + 3 * xe**2 *  self.xa * C_as * (C1**2 * C2)**(1/3) 
              + 2.76 * xe**2 *  self.xy * (C1**2 * C3)**(1/3) 
              + 3 * xe *  self.xa**2 * C_as * (C1 * C2**2)**(1/3) 
              + 6.6 * xe *  self.xa *  self.xy * (C1 * C2 * C3)**(1/3) 
              + 2.76 * xe *  self.xy**2 * (C1 * C3**2)**(1/3) 
              +  self.xa**3 * C2 + 3 *  self.xa**2 *  self.xy * C223 
              + 3 *  self.xa *  self.xy**2 * C233 +  self.xy**3 * C3)

        P_MPa = P * 0.101325
        
        b = 1000 * P_MPa / (2.7715 * self.T)
        
        B0 = b * Bm

        C0 = b**2 * Cm
        
        A0 = 1 + 1.5 * (B0 + C0)

        A1 = 1 + B0
        
        A2 = np.cbrt((A0 - (A0**2 - A1**3)**0.5)) 
        
        Z = (1 + A2 + A1/A2) / 3

        return Z

    def ro(self, P: float) -> float:
        """
        Плотность газа [кг/м³] при давлении P.
        """
        P_pa = P * 101325
        Z = self.z(P)
        Ro = (P_pa * self.M) / (Z * self.R * self.T)

        return Ro

    def bg(self, P: float) -> float:
        """
        Рассчитать объёмный коэффициент расширения газа Bg.
        """
    
        Z = self.z(P)
        
        Bg = (self.Pstd * Z * self.T) / (P * self.Tstd)

        return Bg

    def mu(self, P: float) -> float:
        """
        Рассчитать вязкость газа Mu.
        Mu = интерполяция по табличным данным зависимости вязкости от давления
          
        """
        
        df = pd.read_csv('interp_data.csv', sep=';')

        LI_viscosity = LinearInterpolator(df['pressure, atm'].tolist(),
                                          df['viscosity, cP'].tolist())
        Mu = LI_viscosity.predict(P)

        return Mu
