class LinearInterpolator:
    """
    Линейный интерполятор.
    
    Параметры
    ----------
    x : list
        Узловые точки (отсортированы по возрастанию).
    y : list
        Значения функции в узловых точках.
    """
    
    def __init__(self, x, y):

        self.x = x
        self.y = y

        if len(self.x) != len(self.y):
            raise ValueError("Ошибка: x и y должны иметь одинаковую длину")
    
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

        if xp < self.x[0]:
            return self.y[0]
        if xp > self.x[-1]:
            return self.y[-1]

        for i in range(len(self.x) - 1):
            if self.x[i] <= xp <= self.x[i+1]:
                yp = (self.y[i] + (self.y[i+1] - self.y[i]) /
                    (self.x[i+1] - self.x[i]) * (xp - self.x[i]))
                return yp

        raise ValueError("Не удалось найти интервал для интерполяции")

