import numpy as np
from src.fluid import Fluid
from src.pipe import Pipe

class Well:
  
    beta = 0.00852702
    
    def __init__(self, fluid: Fluid, k: float, h: float, re: float, rw: float, pipe: Pipe = None):

        self.fluid = fluid
        self.k = k          # Проницаемость [мД]
        self.h = h          # Эффективная мощность [м]
        self.re = re        # Радиус контура питания [м]
        self.rw = rw        # Радиус скважины [м]
        self.pipe = pipe    #
    
    def productivity_coef(self, P_res: float, P_bhp: float) -> float:

        P_avg = (P_res + P_bhp)/2 
        mu = self.fluid.mu(P_avg)  # [cP]
        C = (self.beta * self.k * self.h) / (mu * np.log(self.re / self.rw))
        return C
    
    def q(self, P_res: float, P_bhp: float) -> float:
        
        if P_bhp >= P_res:
            return 0.0
        
        C = self.productivity_coef(P_res, P_bhp)
        q_std = C * (P_res - P_bhp) / self.fluid.bg(P_res)
        return q_std
