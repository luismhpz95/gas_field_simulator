import numpy as np
from src.state import NodeState
from src.fluid import Fluid

class Pipe:
    
    g = 9.81
    
    def __init__(self, L: float, D: float, roughness: float, fluid: Fluid, vertical_depth: float = 0.0):

        self.L = L                  # Длина трубы [м]
        self.D = D                  # Внутренний диаметр [м]
        self.roughness = roughness  # Абсолютная шероховатость стенки δ [м]
        self.fluid = fluid
        self.H = vertical_depth     # Вертикальная глубина H [м] (0 для горизонтальной трубы)
    
    def f_lambda(self, Re: float):
       # Коэффициент гидравлического сопротивления λ
       if Re <= 0:
           return 0.0
       
       # Ламинарный: λ = 64/Re
       if Re < 2300:
          return 64/Re
       
       # Турбулентный: уравнение Колбрука–Уайта
       lambda_0 = 0.02
       for _ in range(1000):
           lambda_new = (-2 * np.log10(self.roughness/(3.7*self.D)
                                       + 2.51/(Re*np.sqrt(lambda_0)))) ** -2

            if abs(lambda_new - lambda_0) <= 1e-6:
               return lambda_new

            lambda_0 = lambda_new
        return lambda_new   

    def dp(self, P_in: float, q: float) -> NodeState:  # состояние элемента, q [ст.м³/сут]
        """
        Рассчитать перепад давления в трубе.
        
        P_in - Давление на входе в трубу [атм]
        q_std - Коммерческий расход [ст.м³/сут]
        
        """
        # Свойства газа при давлении на входе
        bg = self.fluid.bg(P_in)
        rho = self.fluid.ro(P_in)
        mu_Pas = self.fluid.mu(P_in)/1000

        v = 4 * q * bg / (np.pi * self.D**2 * 86400)
    
        Re = rho * v * self.D / mu_Pas

        f_lambda = self.f_lambda(Re)        
        
        dP_Pa = (f_lambda * (self.L / self.D * rho * v**2 / 2) 
                  + rho * self.g * self.H)
        dP = dP_Pa/101325 # [атм]
        
        P_out = P_in - dP

        # Состояние элемента трубопровода.
        return NodeState(
            name = f"Pipe(L={self.L}, D={self.D})",
            P_in = P_in,
            P_out = P_out,
            dP = deltaP,
            q_std = q,
            q_res = q * bg,
            v = v,
            rho = rho
        )
