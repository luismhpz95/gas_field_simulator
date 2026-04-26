class LinearInterpolator:
    """
    Линейный интерполятор.
    
    Параметры
    ----------
    xs : list
        Узловые точки (отсортированы по возрастанию).
    ys : list
        Значения функции в узловых точках.
    """
    
    def __init__(self, xs, ys):
        if len(xs) != len(ys):
            raise ValueError("Ошибка: xs и ys должны иметь одинаковую длину")
        
        for i in range(1, len(xs)):
            if xs[i] <= xs[i-1]:
                raise ValueError("Ошибка: xs должен быть отсортирован по возрастанию")
        
        self.xs = xs
        self.ys = ys
    
    def predict(self, xp):
        """
        Вычислить интерполированное значение yp для заданного xp.
        
        Параметры
        ----------
        xp : float
            Точка, в которой нужно найти значение.
        
        Возвращает
        ----------
        float
            Интерполированное значение yp.
        """
        if xp < self.xs[0] or xp > self.xs[-1]:
            raise ValueError(f"xp={xp} вне диапазона [{self.xs[0]}, {self.xs[-1]}]")
        
        for i in range(len(self.xs) - 1):
            if self.xs[i] <= xp <= self.xs[i+1]:
                yp = self.ys[i] + (self.ys[i+1] - self.ys[i]) / (self.xs[i+1] - self.xs[i]) * (xp - self.xs[i])
                return yp
        
        raise ValueError("Не удалось найти интервал для интерполяции")
