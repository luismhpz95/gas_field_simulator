import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from src.state import NodeState
from src.reservoir import Reservoir
from src.well import Well
from src.pipe import Pipe
from src.compressor import DCS


class FieldSimulator:
    """
    Симулятор газового куста из 3 скважин.
    """
    
    def __init__(self, reservoir: Reservoir, wells: list, shlyf: Pipe, dcs: DCS):
        self.reservoir = reservoir
        self.wells = wells
        self.shlyf = shlyf
        self.dcs = dcs
    
    def solve(self, P_res: float) -> dict[str, NodeState]:
        """
        Найти рабочую точку системы.
        """
        
        # Начальное приближение
        q0 = 500.0
        P_man0 = self.dcs.P_in() + 5.0
        x0 = [q0, q0, q0, P_man0]
        
        def equations(x):
            q1, q2, q3, P_man = x
            
            eq1 = self.well_equation(0, P_res, P_man, q1)
            eq2 = self.well_equation(1, P_res, P_man, q2)
            eq3 = self.well_equation(2, P_res, P_man, q3)
            
            q_total = q1 + q2 + q3 + self.dcs.q_ext
            P_man_calc = self.dcs.P_in() + self.shlyf.dp(P_man, q_total).dP
            eq4 = P_man - P_man_calc
            
            return [eq1, eq2, eq3, eq4]
        
        # Решение системы
        solution = fsolve(equations, x0)
        q1, q2, q3, P_man = solution
        
        # Если дебит отрицательный - закрываем скважину
        if q1 < 0:
            q1 = 0.0
        if q2 < 0:
            q2 = 0.0
        if q3 < 0:
            q3 = 0.0
        
        # Создаем словарь для состояний
        states = {}
        
        # Состояния скважин
        for i, (well, q) in enumerate(zip(self.wells, [q1, q2, q3])):
            if well.pipe is not None:
                result_nkt = well.pipe.dp(P_man, q)
                result_nkt.name = f'well_{i+1}'  # ← solo cambias el nombre
                states[f'well_{i+1}'] = result_nkt
            else:
                states[f'well_{i+1}'] = NodeState(
                    name=f'well_{i+1}',
                    P_in=P_man,
                    P_out=P_man,
                    dP=0.0,
                    q_std=q
                    )
        
        # Состояние шлейфа
        q_total = q1 + q2 + q3 + self.dcs.q_ext
        shlyf_result = self.shlyf.dp(P_man, q_total)
        states['shlyf'] = shlyf_result
        
        # Состояние ДКС
        P_in_dcs = self.dcs.P_in()
        states['dcs'] = NodeState(
            name='dcs',
            P_in=P_in_dcs,
            P_out=self.dcs.P_line,
            dP=self.dcs.P_line - P_in_dcs,
            q_std=q_total,
            q_res=None,
            v=None,
            rho=None
        )
        
        return states
    
    def well_equation(self, well_idx: int, P_res: float, P_man: float, q: float) -> float:
        well = self.wells[well_idx]
        
        if well.pipe is not None:
            result_nkt = well.pipe.dp(P_man, q)
            P_bhp = P_man + result_nkt.P_out
        else:
            P_bhp = P_man
        
        q_ipr = well.q(P_res, P_bhp)
        
        return q_ipr - q
    
    def run(self, N_days: int, dt: float = 1.0) -> pd.DataFrame:
        """
        Запустить динамическую симуляцию.
        """
        
        results = []
        P_res = self.reservoir.resprops.P
        Gp = 0.0
        
        print("=" * 60)
        print("ЗАПУСК СИМУЛЯЦИИ")
        print("=" * 60)
        
        for day in range(1, N_days + 1):
            states = self.solve(P_res)
            
            q1 = states['well_1'].q_std
            q2 = states['well_2'].q_std
            q3 = states['well_3'].q_std
            q_total = q1 + q2 + q3
            P_man = states['well_1'].P_in
            
            P_res_new = self.reservoir.p2(q_total, dt)
            Gp += q_total * dt / 1000
            
            results.append({
                't': day,
                'P_res': P_res,
                'P_man': P_man,
                'q1': q1,
                'q2': q2,
                'q3': q3,
                'q_total': q_total,
                'Gp': Gp
            })

            self.reservoir.resprops.P = P_res_new
            P_res = P_res_new
            
            if day % 10 == 0 or day == 1:
                print(f"День {day:3d} | P_res = {P_res:6.2f} атм | q_total = {q_total:6.2f} м³ст/сут | Gp = {Gp:.1f} тыс. м³")
        
        print("=" * 60)
        print(f"СИМУЛЯЦИЯ ЗАВЕРШЕНА. Дней: {N_days}")
        print("=" * 60)
        
        return pd.DataFrame(results)
