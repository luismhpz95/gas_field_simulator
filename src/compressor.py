class DCS:
    """
    Дожимная компрессорная станция (ДКС).
    Повышает давление газа до давления магистрали.
    """
    
    def __init__(self, CR: float, P_line: float, q_ext: float = 0.0):
        # CR — степень сжатия (≥ 1.0), P_line [атм]
        # q_ext — расход стороннего газа, поступающего на манифолд [ст.м³/сут]
        
        self.CR = CR
        self.P_line = P_line
        self.q_ext = q_ext
    
    def P_in(self) -> float:
        """
        Давление на входе в ДКС.
        При CR = 1.0 станция отключена, возвращается P_line.
        
        Returns: Давление на входе [атм]
        """
        return self.P_line / self.CR