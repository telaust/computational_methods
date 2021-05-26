import matplotlib.pyplot as plt # для графика
import numpy as np # в данном коде тоже для графика
from math import *
import numexpr as ne


# Формула трапеций

def trapezoidal(func, a, b, n):

    h = (b - a) / n
    sum = 0
    x = a
    for i in range(1, n-1):
        sum += ne.evaluate(func)
        x += h

    return h * (ne.evaluate(func, {'x': a})/2 + ne.evaluate(func, {'x': b})/2 + sum)


# Формула Симсона

def simpson(func, a, b, n):

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


# из трех вариантов метода пряямоугольника, я буду использовать средний, потому что
# при прочих равных условиях он дает приближение точнее, чем левый и правый

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


# чтение из файла - 3 примера
with open("input.txt", 'r') as infile:
    eq = str(infile.readline())

    line0 = infile.readline() # pass the line

    a_file1 = int(next(infile))
    b_file1 = int(next(infile))
    x0_file1 = float(next(infile))
    eps1 = float(next(infile))

    line1 = infile.readline() # pass the line

    a_file2 = int(next(infile))
    b_file2 = int(next(infile))
    x0_file2 = float(next(infile))
    eps2 = float(next(infile))

    line2 = infile.readline() # pass the line

    a_file3 = int(next(infile))
    b_file3 = int(next(infile))
    x0_file3 = float(next(infile))
    eps3 = float(next(infile))


# f = lambda x: eval(eq)

space = np.linspace(-10, 10, 100)
# f2 = np.vectorize(f)

# чертим график
# plt.plot(space, f2(space), color='#fcba03', label=eq)
# plt.legend()
# plt.grid(True)
# plt.suptitle("Graph")
# plt.show()


# функция решения нелинейного уравнения

def solve_nonlinear_equation(a, x0, func, b, eps=1e-7):
    # принимает на вход:
    # a - нижний предел интегрирования
    # x0 - начальное прилижение (наша догадка)
    # func - подыинтегральная функци
    # b - верхний предел интегрирования
    # eps - допустимая погрешность

    x = x0
    prev = 0

    # while abs(prev - x) > eps:
    #     prev = x
    #
    #     # g = riemann_sum(func, a, x, 100_000) - b
    #
    #     g = simpson(func, a, x, 1000) - b
    #     # print(g)
    #
    #     # если в какой-то точке функция обнуляется, то приращаю х, чтобы избежать исключения
    #     if ne.evaluate(func) == 0:
    #         x += 1e-9
    #
    #     x = x - g / ne.evaluate(func)



    for i in range(1_00):

        g = riemann_sum(func, a, x, 10_000) - b

        # если в какой-то точке функция обнуляется, то приращаю х, чтобы избежать исключения
        if ne.evaluate(func) == 0:
            x += 1e-9

        x = x - g / ne.evaluate(func)

        if abs(prev - x) < eps:
            break

        prev = x

    return x


sol1 = solve_nonlinear_equation(a_file1, x0_file1, eq, b_file1, eps=eps1)
# sol2 = solve_nonlinear_equation(a_file2, x0_file2, f, b_file2, eps=eps2)
# sol3 = solve_nonlinear_equation(a_file3, x0_file3, f, b_file3, eps=eps3)

print(sol1)

# запись в файл и округление
with open('out.txt', 'w') as f:
    f.write(str(round(sol1, 4)) + "\n")

f.close()

