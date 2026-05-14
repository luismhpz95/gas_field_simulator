from src.fluid import Fluid

@dataclass
class ResProps:
    """
    Свойства пласта (текущее состояние)
    """
    P: float    # давление [атм]
    V: float    # объем [м³]
    T: float    # температура [К]


class Reservoir:
    """
    Газовый пласт с материальным балансом.
    """
    
    Pstd = 1.0       # стандартное давление [атм]
    Tstd = 293.15    # стандартная температура [К]
    Z_std = 1

    def __init__(self, resprops: ResProps, fluid: Fluid):
        """
        Начальное состояние пласта (P [атм], V [м³], T [К])
        """

        self.resprops = resprops
        self.fluid = fluid
        
        # Расчет плотности в стандартных условиях
        
        P_Pa = self.Pstd * 101325
        self.ro_std = (P_Pa  * fluid.M) / (self.Z_std * fluid.R * self.Tstd)
    
    def p2(self, q_total: float, dt: float = 1.0) -> float:
        """
        Рассчитать давление в пласте на следующем шаге по времени.
        
        Параметры
        ----------
        q_total - Суммарный дебит скважин [ст.м³/сут]
        dt -  Шаг по времени [сутки]
           
        P_res на следующем шаге [атм]
        
        """
        P_res = self.resprops.P
        V_res = self.resprops.V
        
        ro_res = self.fluid.ro(P_res)
        
        Z = self.fluid.z(P_res)
        
        # Формула материального баланса согласно разделу 2.3:
        P2 = P_res - (Z * self.ro_std / ro_res) * (q_total / V_res) * dt
        
        return P2
