import numexpr as ne
import numpy as np

# Формула трапеций

def trapezoidal(func, a, b, n):
    if not isinstance(func, str):
        raise ValueError("{} must be string".format(func))

    h = (b - a) / n
    sum = 0
    x = a
    for i in range(1, n-1):
        sum += ne.evaluate(func)
        x += h

    return h * (ne.evaluate(func, {'x': a})/2 + ne.evaluate(func, {'x': b})/2 + sum)


# Формула Симпсона

def simpson(func, a, b, n):
    if not isinstance(func, str):
        raise ValueError("{} must be string".format(func))

    h = (b - a) / n
    k = 0.0
    x = a + h
    for i in range(1, int(n/2 + 1)):
        k += 4 * ne.evaluate(func)
        x += 2 * h

    x = a + 2 * h

    for i in range(1, int(n/2)):
        k += 2 * ne.evaluate(func)
        x += 2 * h

    return (h / 3) * (ne.evaluate(func, {'x':a}) + ne.evaluate(func, {'x': b}) + k)


def riemann_sum(func, a, b, num_of_partitions):
    # на вход принимает:
    # func - подынтегральную функцию
    # a - ниж предел интегрирования
    # b - верхний предел интегрирования
    # num_of_partitions - число разбиений

    if not isinstance(func, str):
        raise ValueError("{} must be string".format(func))

    step = (b - a) / num_of_partitions # определяем шаг

    x = a
    func_values = []

    for i in range(num_of_partitions):
        x += step/2
        func_values.append(ne.evaluate(func))
        x = x + step/2

    return step * sum(func_values)


# Монте-Карло: интеграл от ф-ии одной переменной

def monte_carlo(func, a, b, trials):
    sum = 0

    for i in range(trials):
        r = np.random.uniform(a, b)
        # you also can use x = 2r+1 here but generate uniform distribution from 0 to 1 on line above
        x = r
        sum += ne.evaluate(func)
    return (b-a) * sum / trials

