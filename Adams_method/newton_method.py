from math import *

# из трех вариантов метода пряямоугольника, я буду использовать средний, потому что
# при прочих равных условиях он дает приближение точнее, чем левый и правый

def riemann_sum(func, a, b, num_of_partitions):
    # на вход принимает:
    # func - подынтегральную функцию
    # a - ниж предел интегрирования
    # b - верхний предел интегрирования
    # num_of_partitions - число разбиений

    step = (b - a) / num_of_partitions # определяем шаг

    x = a
    func_values = []

    for i in range(num_of_partitions):
        func_values.append(func(x + step/2))
        x = x + step

    return step * sum(func_values)


def solve_nonlinear_equation(a, x0, func, b, eps):
    # принимает на вход:
    # a - нижний предел интегрирования
    # x0 - начальное прилижение (наша догадка)
    # func - подыинтегральная функци
    # b - верхний предел интегрирования
    # eps - допустимая погрешность

    x = x0
    prev = x0

    for i in range(1_000): # я думаю лимит в 100 итераций будет достаточен
        g = riemann_sum(func, a, x, 1_000) - b

        # если в какой-то точке функция обнуляется, то поднимаю исключение
        if func(x) == 0:
            raise Exception("it is impossible to converge\n")
        x = x - g / func(x)

        if abs(prev - x) < eps:
            break
        prev = x

    return x

