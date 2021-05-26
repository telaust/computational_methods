from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from eulers_method import *
from newton_method import solve_nonlinear_equation
from math import *
import numexpr as ne
from runge_kutta4 import runge_kutta4
import xlwt


""" Reading input file / Считывание с файла """

def passnlines(filename, n):
    for _ in range(n):
        filename.readline()


with open("input2.txt", "r") as file:

    try:
        passnlines(file, 4)

        f = str(file.readline()).strip("\n")

        passnlines(file, 3)

        y0 = float(next(file))

        passnlines(file, 3)

        interval = next(file)
        # t = float(next(file))
        # T = float(next(file))

        passnlines(file, 3)

        h = float(next(file))

    except Exception:
        exit("Incorrect input data\n")


interval = [float(i) for i in interval.split()]
t = interval[0]
T = interval[1]

# print data / просто вывод

print("f = {}, y0 = {}, t = {}, T = {}, h = {}".format(f, y0, t, T, h))

# """ Check args """



def round_float(num, limit=2):
    if limit <= 0:
        raise ValueError("limit = {} must be > 0".format(limit))

    return float(("%." + str(limit) + "f") % num)

# function value / возвращает значение функции (если x и y не заданы, то по умолчанию 0 и 0)

def value(f, x=0.0, y=0.0):
    res = 0
    try:
        # eval(f)
        res = ne.evaluate(f)
    except OverflowError:
        print("Some points don't belong to the domain of your function\n")
    if res == float("inf"):
        x += 1e-9
        y += 1e-9
        print("Кажется, где-то было деление на нуль. Нехорошо. Я добавил погрешность, чтобы этого избежать")
        res = ne.evaluate(f)

    return ne.evaluate(f)


# Это для разрешения неявной схемы, то есть, итераицонно приближаю у
# метод Адамса (простые итерации)
def iterations(f, step, x_n, y_1, x_1, y_2, x_2, n_iterations=2, eps=1e-4):

    # y = find_y0_euler(f, step, x_prev=x_1, y_prev=y_1)  # инициализация методом Эйлера


    # y = find_y0_adams2(f, step, y_prev1=y_1, y_prev2=y_2) # неявным методом Адамса 2ух шаговым
    y = adams2(f, step=step, x_prev1=x_1, y_prev1=y_1, x_prev2=x_2, y_prev2=y_2)

    y_prev = y
    for s in range(4):
        y = y + step * (1 / 12) * \
            (5 * value(f, x_n, y) + 8 * value(f, x_1, y_1) - value(f, x_2, y_2))
        if abs(y_prev - y) < 1e-6:
            break
        y_prev = y

    return y


# Для инициализации функции выше использую метод Адамса 1-го порядка, он же метод Эйлера

def find_y0_euler(f, step, x_prev, y_prev):
    return y_prev + step * value(f, x_prev, y_prev)


def find_y0_adams2(f, step, y_prev1, y_prev2):
    return adams_explicit2(f, y_prev1, y_prev2, t, T, step)

def adams2(f, step, x_prev1, y_prev1, x_prev2, y_prev2):

    return y_prev1 + (step / 2) * (3*value(f, x_prev1, y_prev1) - value(f, x_prev2, y_prev2))


# def find_y0_adams3(f, step, y_prev1)

# будем накапливать у и х
y_values = [y0]
x_values = [t]

def adams_method(f, y0, t, T, h, eps=1e-5, steps=2, return_values=True):
    x = t  # 0 . мы решаем задачу Коши, где начальное значение задается в нуле
    y = y0  # само начальное значение

    n_partitions = int((T - t) / h)  # кол-во сегментов
    print("T = ", T)
    print("t = ", t)
    print(n_partitions, " число сегментов")

    # Считаю методом Р-К 4-го порядка первые 2 значения

    for i in range(steps+1):

        """  Print """
        # print("x = ", x_values[-1])
        # print("y = ", y_values[-1])



        k1 = value(f, x, y)
        k2 = value(f, x + h / 2, y + (h * k1) / 2)
        k3 = value(f, x + h / 2, y + (h * k2) / 2)
        k4 = value(f, x + h, y + h * k3)

        delta_y = (1 / 6) * h * (k1 + 2 * k2 + 2 * k3 + k4)

        y += delta_y  # update y
        x += h  # update x

        y_values.append( round_float( y , 6) )
        x_values.append( round_float( x , 6))

    # Потом досчитываю 2шаговым методом Адамса, где использую для разрешения
    # неявной схемы простые итерации, где в качестве начального (нулевого значения) беру решение
    # явного метода Адамса 1-го порядка (он же метод Эйлера)

    # print("после РК4 x= ", x)
    # print("после РК4 y= ", y)


    for i in range(n_partitions - steps -1):  # тут было +1
        # Приближенно считаем y_n чтобы подставить его в f_n(x_n, y_n) в правую часть

        if isnan(y_values[-1]) or isnan(x_values[-1]):
            # exit("Кажется, вы ввели невозможное условие или функцию")
            exit("Кажется, вы ввели невозможное условие или функцию")

        appr_y = iterations(f,
                            step=h,
                            x_n=x,
                            y_1=y_values[-1],
                            x_1=x_values[-1],
                            y_2=y_values[-2],
                            x_2=x_values[-2])

        # print("итерация #", i, " у(прибл, до) = ", appr_y)

        delta_y = h * (1 / 12) * (5 * value(f, x=x, y=appr_y) +
                                  8 * value(f, x=x_values[-1], y=y_values[-1]) -
                                  value(f, x=x_values[-2], y=y_values[-2]))

        # print("dy = ", delta_y)

        if isnan(delta_y):
            # raise SystemExit("Кажется, вы ввели невозможное условие или функцию")
            exit("Кажется, вы ввели невозможное условие или функцию")  # Bad manner

        x += h  # update x
        y += delta_y  # update y

        # print("итерация #", i, " у(прибл, после) = ", y)

        y_values.append( y )
        x_values.append( round_float(x, 5) )

        # print("this y ", y_values[-1])
        # print("prev y ", y_values[-2])

        # if abs(y_values[-1] - y_values[-2]) < eps:
        #     print("Нужная точность достигнута\n")
        #     break

        if abs(y_values[-1] - y_values[-2]) > 1e+3:
            print("Расходится\n")
            break

    return (x, y) if not return_values else (x, y, x_values, y_values)


sol = adams_method(f, y0, t, T, h)


""" Plotting """

plt.figure(num=None, figsize=(7, 5))
plt.title("Решение ОДУ Коши\n")


plt.plot(x_values, y_values, label="Неявная схема Адамса, 3 порядок")
plt.legend(loc='best')
plt.grid(True)

# Решение Р-К4 для сравнения

# y, vals = runge_kutta4(f, y0, t, T, h, True)

# plt.plot(x_values, vals, label="Явная схема Р-К, 4 порядок")
# plt.legend()
# plt.grid(True)

plt.show()

""" Write solution to file """

with open("output.txt", "w") as file:
    file.write("Решение:\n\n")
    file.write("y(" + str(T) + ") = " + str(sol[1]) + "\n")


file.close()

###################################################################################################
###################################################################################################
###################################################################################################


analytical_solution = []
n_partitions = int((T-t)/h)

x = t

# str_fn_sol = '(x+1)*exp(-x)'
str_fn_sol = '1/(exp(x))'
# str_fn_sol = 'sqrt(25-2*x) + x - 2'


for i in range(n_partitions+1):
    analytical_solution.append(ne.evaluate(str_fn_sol).tolist())
    x += h




###################################################################################################
###################################################################################################
###################################################################################################


""" Write to Excel """


def write_to_excel(filename, x, y):
    wb = xlwt.Workbook()

    sheet = wb.add_sheet("sheet1")

    sheet.write(0, 0, "x")
    sheet.write(0, 1, "y")

    sheet.write(0, 2, "точный y")

    for i in range(len(x)):
        sheet.write(i+1, 0, str(x[i]))

    for i in range(len(y)):
        sheet.write(i+1, 1, str(y[i]))

    for i in range(len(analytical_solution)):
        sheet.write(i+1, 2, str(analytical_solution[i]))

    wb.save(filename)


write_to_excel("output_another.xls", x_values, y_values)



