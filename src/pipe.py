import numpy as np
from src.state import NodeState

class Pipe:
    
    g = 9.81
    
    def __init__(self, L: float, D: float, roughness: float, fluid: Fluid, vertical_depth: float = 0.0):
        # L [м], D [м], roughness (δ) [м], vertical_depth (H) [м]

        self.L = L
        self.D = D 
        self.roughness = roughness
        self.fluid = fluid
        self.H = vertical_depth
    
    def f_lambda(self, Re: float):
       # Коэффициент трения λ рассчитывать по уравнению Колбрука–Уайта 
       
       if Re < 2300:
            
        return 64/Re
       
       lambda_new = 0
       lambda_0 = 0.02

       while True:
        
        lambda_new = (-2 * np.log10(self.roughness/(3.7*self.D)
                                         +2.51/(Re*np.sqrt(lambda_0))))**-2
        
        if abs(lambda_new - lambda_0) <= 1e-6:
            return lambda_new
        
        lambda_0 = lambda_new

    def dp(self, P_in: float, q: float) -> NodeState:  # состояние элемента, q [ст.м³/сут]
        
        bg = self.fluid.bg(P_in)
        rho = self.fluid.ro(P_in)
        mu_Pas = self.fluid.mu(P_in)/1000

        v = 4 * q * bg / (np.pi * self.D**2 * 86400)
    
        Re = rho * v * self.D / mu_Pas

        f_lambda = self.f_lambda(Re)        
        
        deltaP = (f_lambda * (self.L/self.D*rho*v**2/2) + rho*self.g*self.H)/101325

        P_out = P_in - deltaP

            
        return NodeState(
            name=f"Pipe(L={self.L}, D={self.D})",
            P_in=P_in,
            P_out=P_out,
            dP=deltaP,
            q_std=q,
            q_res=q*bg,
            v=v,
            rho=rho
    )
