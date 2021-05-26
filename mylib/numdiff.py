
""" Numerical differentiation """

import numexpr as ne
from matvec import *
import numpy as np
import string
import math


# Производные по опредлению
def derivativeX(func, X0, dx=1e-5, clear=True):

    if type(X0) is not list:
        raise Exception("{} must be list".format(X0))

    x = X0[0] + dx
    y = X0[1]

    if len(X0) == 3:
        z = X0[2]

    f1 = ne.evaluate(func)

    x -= dx

    f2 = ne.evaluate(func)

    return round_float((f1 - f2) / dx) if clear else (f1 - f2) / dx


# аналог для одномерной ф-ии
def derivativeX_(func, x, dx=1e-5, clear=True):

    f1 = ne.evaluate(func)
    x -= dx
    f2 = ne.evaluate(func)

    return round_float((f1 - f2) / dx) if clear else (f1 - f2) / dx


def derivativeY(func, X0, dy=1e-5, clear=True):

    if type(X0) is not list:
        raise Exception("{} must be list".format(X0))

    x = X0[0]
    y = X0[1] + dy

    if len(X0) == 3:
        z = X0[2]

    f1 = ne.evaluate(func)

    y -= dy

    f2 = ne.evaluate(func)


    return round_float((f1 - f2) / dy) if clear else (f1 - f2) / dy


def derivativeZ(func, X0, dz=1e-5, clear=True):

    if type(X0) is not list:
        raise Exception("{} must be list".format(X0))

    x = X0[0]
    y = X0[1]
    z = X0[2] + dz

    f1 = ne.evaluate(func)

    z -= dz

    f2 = ne.evaluate(func)


    return round_float((f1 - f2) / dz) if clear else (f1 - f2) / dz


# > 3
def derivative(func, x_list, index_of_derivative=0, dx=1e-5, clear=True):

    variables = string.ascii_lowercase[:len(x_list)]  # first n
    deriv_var = variables[index_of_derivative]

    # Fill local dict
    dict = {}

    keys = list(variables)
    values = x_list

    for i, key in enumerate(keys):
        if i == index_of_derivative:
            dict[key] = values[i] + dx  # даю приращение
        else:
            dict[key] = values[i]


    f1 = ne.evaluate(func, local_dict=dict)

    dict2 = {}
    for i, key in enumerate(keys):
        dict2[key] = values[i]

    f2 = ne.evaluate(func, dict2)


    return ((np.array(f1) - np.array(f2)) / dx) if not clear \
        else round_float((np.array(f1) - np.array(f2)) / dx)


# 2-dimensional jacobian

def jacobian2(eqs, X, clear=True):

    if type(eqs) is not list:
        raise Exception("{} must be list".format(eqs))

    n_eq = len(eqs)
    jac = zeromatrix(n_eq, n_eq) # квадратная матрица NxN

    jac[0][0] = derivativeX(eqs[0], X)
    jac[0][1] = derivativeY(eqs[0], X)

    jac[1][0] = derivativeX(eqs[1], X)
    jac[1][1] = derivativeY(eqs[1], X)

    return round_matrix(jac) if clear else jac


# 3-dimensional jacobian

def jacobian3(eqs, X, clear=True):

    if type(eqs) is not list:
        raise Exception("{} must be list".format(eqs))

    n_eq = len(eqs)
    jac = zeromatrix(n_eq, n_eq) # квадратная матрица NxN

    jac[0][0] = derivativeX(eqs[0], X)
    jac[0][1] = derivativeY(eqs[0], X)
    jac[0][2] = derivativeZ(eqs[0], X)

    jac[1][0] = derivativeX(eqs[1], X)
    jac[1][1] = derivativeY(eqs[1], X)
    jac[1][2] = derivativeZ(eqs[1], X)

    jac[2][0] = derivativeX(eqs[2], X)
    jac[2][1] = derivativeY(eqs[2], X)
    jac[2][2] = derivativeZ(eqs[2], X)

    return round_matrix(jac) if clear else jac


# N-dimensional jacobian

def jacobianN(eqs, X):

    if len(X) != len(eqs):
        raise ValueError("Dimension mismatch; {} must be equal to".format(len(X), len(eqs)))

    n_eq = len(eqs)
    jac = zeromatrix(n_eq, n_eq) # квадратная матрица NxN

    for i in range(n_eq):
        for j in range(n_eq):
            jac[i][j] = derivative(eqs[i], X, index_of_derivative=j)

    return jac



def valueof2(eqs, X0):
    if len(X0) != 2:
        raise ValueError("X0: expected 2, given {}".format(len(X0)))

    n_eq = len(eqs)
    res = zerovector(n_eq)

    x = X0[0]
    y = X0[1]

    for i in range(len(res)):
        res[i] = ne.evaluate(eqs[i])

    return np.array(res).tolist()

def value(f, x=0.0, y=0.0): return ne.evaluate(f)

def valueof3(eqs, X0):
    if len(X0) != 3:
        raise ValueError("X0: expected 3, given {}".format(len(X0)))

    n_eq = len(eqs)
    res = zerovector(n_eq)

    x = X0[0]
    y = X0[1]
    z = X0[2]

    for i in range(len(res)):
        res[i] = ne.evaluate(eqs[i])

    return np.array(res).tolist()


# Значение многомерной функции

def valueofN(eqs, X0):
    n_eq = len(eqs)
    res = zerovector(n_eq)

    variables = string.ascii_lowercase[:n_eq]

    # Fill local dict
    dict = {}

    keys = list(variables)
    values = X0

    for i, key in enumerate(keys):
        dict[key] = values[i]

    for i in range(len(res)):
        res[i] = ne.evaluate(eqs[i], local_dict=dict)

    return res


def how_much_vars(func_str):
    vars = 0

    for i in range(len(func_str)):
        if func_str[i] in string.ascii_lowercase:
            vars += 1

    return vars

def derivative_n(func, x, n):
    res = func

    for _ in range(n):
        res = ne.evaluate(func)

    return res

import matplotlib.pyplot as plt

def taylor(func, x, a, n=2):

    series = 0
    current_derivative = func
    for i in range(1, n):
        series += derivativeX_(current_derivative, x) * (x - a)**i / math.factorial(i)

        current_derivative = derivativeX_(current_derivative,  x)

    return series + ne.evaluate(func, {'x': x})  # прибавляем первый член последовательности


class DiffSolver(object):
    """

    this class defines methods for solving Cauchy problem for only two-dimensional functions, eg f(x, y)

    """
    def __init__(self, function, y0, a, b, h):
        """
        Agrs:
        :param function: function which will differentiate
        :param y0: initial y value
        :param a: left bound
        :param b: right bound
        :param h: step
        """

        assert isinstance(function, str)
        assert a < b
        assert h > 0

        self.function = function
        self.y0 = y0
        self.a = a
        self.b = b
        self.h = h
        self.n_partitions = int((self.b - self.a) / self.h)

    def euler(self, plot=False, return_values=False):

        x = self.a
        y = self.y0

        y_values = [self.y0]

        for i in range(self.n_partitions):
            f = value(self.function, x, y)

            x += self.h
            y += self.h * f

            y_values.append(y)

        if plot:
            plt.grid(True)
            plt.plot([(i + self.h) for i in range(self.n_partitions)], y_values[1:])
            plt.show()

        return y if not return_values else (y, y_values)

    def runge_kutta(self, order=4, plot=False, return_values=False):

        assert 0 < order < 5

        x = self.a
        y = self.y0

        y_values = [self.y0]

        if order == 1:
            print("it is Euler's method.. high inaccuracy")
            return self.euler()

        if order == 2:
            for i in range(self.n_partitions):
                k1 = value(self.function, x, y)
                k2 = value(self.function, x + self.h, y + self.h * k1)

                delta_y = (1 / 2) * self.h * (k1 + k2)

                x += self.h
                y += delta_y

                y_values.append(y)

            if plot:
                plt.title("Runge-Kutta 2nd order")

        if order == 3:

            for i in range(self.n_partitions):
                k1 = value(self.function, x, y)
                k2 = value(self.function, x + self.h/2, y + self.h * k1 / 2)
                k3 = value(self.function, x + self.h, y - self.h*k1 + 2*self.h*k2)

                delta_y = (1 / 6) * self.h * (k1 + 4*k2 + k3)

                x += self.h
                y += delta_y

                # x_values.append(x)
                y_values.append(y)
            if plot:
                plt.title("Runge-Kutta 3rd order")

        if order == 4:
            for i in range(self.n_partitions):
                k1 = value(self.function, x, y)
                k2 = value(self.function, x + self.h / 2, y + (self.h * k1) / 2)
                k3 = value(self.function, x + self.h / 2, y + (self.h * k2) / 2)
                k4 = value(self.function, x + self.h, y + self.h * k3)

                delta_y = (1 / 6) * self.h * (k1 + 2 * k2 + 2 * k3 + k4)

                x += self.h
                y = y + delta_y

                y_values.append(y)

            if plot:
                plt.title("Runge-Kutta 4th order")

        if plot:
            plt.grid(True)
            plt.plot([(i + self.h) for i in range(self.n_partitions)], y_values[1:])
            plt.show()

        return y if not return_values else (y, y_values)

    def iterations1(self, x_n, y_1, x_1, n_iterations=3):

        y = y_1 + self.h * value(self.function, x_1, y_1)

        for s in range(n_iterations):
            y += self.h * (1 / 2) * (value(self.function, x_n, y) + value(self.function, x_1, y_1))

        return y

    def iterations(self, x_n, y_1, x_1, y_2, x_2, n_iterations=3):

        y = y_1 + (self.h / 2) * (3 * value(self.function, x_1, y_1) - value(self.function, x_1, y_1))

        for s in range(n_iterations):
            y += self.h * (1 / 12) * \
                (5 * value(self.function, x_n, y) + 8 * value(self.function, x_1, y_1) - value(self.function, x_2, y_2))

        return y

    def iterations3(self, x_n, y_1, x_1, y_2, x_2, y_3, x_3, n_iterations=4):
        y = y_1 + (self.h / 12) * (23 * value(self.function, x_1, y_1) -
                                   16 * value(self.function, x_2, y_2) + 5*value(self.function, x_3, y_3))

        for s in range(n_iterations):
            y += (self.h / 24) * (55 * value(self.function, x_n, y) -
                                  59 * value(self.function, x_1, y_1) + 37 * value(self.function, x_2, y_2) - 9 * value(self.function, x_3, y_3))

        return y

    def iterations4(self, x_n, y_1, x_1, y_2, x_2, y_3, x_3, x_4, y_4, n_iterations=2) :
        y = y_1 + (self.h / 24) * (55 * value(self.function, x_1, y_1) - 59 * value(self.function, x_2, y_2) +
                                   37 * value(self.function, x_3, y_3) - 9 * value(self.function, x_4, y_4))

        for s in range(n_iterations):
            y += (self.h / 720) * (1901 * value(self.function, x_n, y) -
                                   2774 * value(self.function, x_1, y_1) +
                                   2616 * value(self.function, x_2, y_2) -
                                   1274 * value(self.function, x_3, y_3) +
                                   251 * value(self.function, x_4, y_4))

        return y

    # Explicit only method

    def adams_bashforth(self, order, plot=False, return_values=False):
        return self.adams(order=order, explicit=True, plot=plot, return_values=return_values)

    # Implicit only method

    def adams_moulton(self, order, plot=False, return_values=False):
        return self.adams(order=order, explicit=False, plot=plot, return_values=return_values)

    def adams(self, order, explicit=True, plot=False, return_values=False):

        x = self.a
        y = self.y0

        x_values, y_values = [], []

        if explicit:

            if order == 2:
                for i in range(order):
                    k1 = value(self.function, x, y)
                    k2 = value(self.function, x + self.h / 2, y + (self.h * k1) / 2)
                    k3 = value(self.function, x + self.h / 2, y + (self.h * k2) / 2)
                    k4 = value(self.function, x + self.h, y + self.h * k3)

                    delta_y = (1 / 6) * self.h * (k1 + 2 * k2 + 2 * k3 + k4)

                    y += delta_y
                    x += self.h

                    y_values.append(round_float(y, 6))
                    x_values.append(round_float(x, 6))

                for i in range(self.n_partitions - order):
                    k1 = (3 / 2) * value(self.function, x, y_values[-1])
                    k2 = (- 1 / 2) * value(self.function, x, y_values[-2])

                    delta_y = (k1 + k2) * self.h

                    y += delta_y
                    x += self.h

                    y_values.append(y)
                    x_values.append(x)

            if order == 3:
                for i in range(order):
                    k1 = value(self.function, x, y)
                    k2 = value(self.function, x + self.h / 2, y + (self.h * k1) / 2)
                    k3 = value(self.function, x + self.h / 2, y + (self.h * k2) / 2)
                    k4 = value(self.function, x + self.h, y + self.h * k3)

                    delta_y = (1 / 6) * self.h * (k1 + 2 * k2 + 2 * k3 + k4)

                    y += delta_y
                    x += self.h

                    y_values.append(round_float(y, 6))
                    x_values.append(round_float(x, 6))

                for i in range(self.n_partitions - order):
                    k1 = 23 * value(self.function, x, y_values[-1])
                    k2 = -16 * value(self.function, x, y_values[-2])
                    k3 = 5 * value(self.function, x, y_values[-3])

                    delta_y = (1 / 12) * self.h * (k1 + k2 + k3)

                    y += delta_y
                    x += self.h

                    y_values.append(y)
                    x_values.append(x)


            if order == 4:
                for i in range(order):
                    k1 = value(self.function, x, y)
                    k2 = value(self.function, x + self.h / 2, y + (self.h * k1) / 2)
                    k3 = value(self.function, x + self.h / 2, y + (self.h * k2) / 2)
                    k4 = value(self.function, x + self.h, y + self.h * k3)

                    delta_y = (1 / 6) * self.h * (k1 + 2 * k2 + 2 * k3 + k4)

                    y += delta_y
                    x += self.h

                    y_values.append(round_float(y, 6))
                    x_values.append(round_float(x, 6))

                for i in range(self.n_partitions - order):
                    k1 = 55 * value(self.function, x, y_values[-1])
                    k2 = -59 * value(self.function, x, y_values[-2])
                    k3 = 37 * value(self.function, x, y_values[-3])
                    k4 = -9* value(self.function, x, y_values[-4])

                    delta_y = (1 / 24) * self.h * (k1 + k2 + k3 + k4)

                    y += delta_y
                    x += self.h

                    y_values.append(y)
                    x_values.append(x)

                return y if not return_values else (y, y_values)

            if order == 5:

                # calculate first 5 values using Runge-Kutta method
                for i in range(order):
                    k1 = value(self.function, x, y)
                    k2 = value(self.function, x + self.h / 2, y + (self.h * k1) / 2)
                    k3 = value(self.function, x + self.h / 2, y + (self.h * k2) / 2)
                    k4 = value(self.function, x + self.h, y + self.h * k3)

                    delta_y = (1 / 6) * self.h * (k1 + 2 * k2 + 2 * k3 + k4)

                    y += delta_y
                    x += self.h

                    y_values.append(round_float(y, 6))
                    x_values.append(round_float(x, 6))

                for i in range(self.n_partitions - order):

                    k1 = 1901 * value(self.function, x,  y_values[-1])
                    k2 = -2774 * value(self.function, x, y_values[-2])
                    k3 = 2616 * value(self.function, x, y_values[-3])
                    k4 = -1274 * value(self.function, x, y_values[-4])
                    k5 = 251 * value(self.function, x, y_values[-5])

                    delta_y = (1 / 720) * self.h * (k1 + k2 + k3 + k4 + k5)

                    y += delta_y
                    x += self.h

                    y_values.append(y)
                    # x_values.append(x)

                if plot:
                    plt.title("Adams Method, order " + str(order))


        if not explicit:

            if order == 2:
                steps = order - 1

                # calculate first one value using Runge-Kutta method

                for i in range(steps + 1):

                    k1 = value(self.function, x, y)
                    k2 = value(self.function, x + self.h / 2, y + (self.h * k1) / 2)
                    k3 = value(self.function, x + self.h / 2, y + (self.h * k2) / 2)
                    k4 = value(self.function, x + self.h, y + self.h * k3)

                    delta_y = (1 / 6) * self.h * (k1 + 2 * k2 + 2 * k3 + k4)

                    y += delta_y
                    x += self.h

                    y_values.append(round_float(y, 6))
                    x_values.append(round_float(x, 6))

                for i in range(self.n_partitions - steps - 1):
                    appr_y = self.iterations1(
                                        x_n=x,

                                        y_1=y_values[-1],
                                        x_1=x_values[-1])

                    delta_y = self.h * (1 / 2) * (value(self.function, x, appr_y) +
                                                  value(self.function, x_values[-1], y_values[-1]))

                    x += self.h
                    y += delta_y

                    y_values.append(y)
                    x_values.append(x)


            if order == 3:
                steps = order - 1

                # Считаю методом Р-К 4-го порядка первые 2 значения
                for i in range(steps + 1):

                    k1 = value(self.function, x, y)
                    k2 = value(self.function, x + self.h / 2, y + (self.h * k1) / 2)
                    k3 = value(self.function, x + self.h / 2, y + (self.h * k2) / 2)
                    k4 = value(self.function, x + self.h, y + self.h * k3)

                    delta_y = (1 / 6) * self.h * (k1 + 2 * k2 + 2 * k3 + k4)

                    y += delta_y
                    x += self.h

                    y_values.append(round_float(y, 6))
                    x_values.append(round_float(x, 6))

                # Потом досчитываю 2шаговым методом Адамса, где использую для разрешения
                # неявной схемы простые итерации, где в качестве начального (нулевого значения) беру решение
                # явного метода Адамса 1-го порядка (он же метод Эйлера)

                for i in range(self.n_partitions - steps - 1):  # тут было +1
                    # Приближенно считаем y_n чтобы подставить его в f_n(x_n, y_n) в правую часть

                    if math.isnan(y_values[-1]) or math.isnan(x_values[-1]):
                        exit("Кажется, вы ввели невозможное условие или функцию")

                    appr_y = self.iterations(
                                        x_n=x,
                                        y_1=y_values[-1],
                                        x_1=x_values[-1],
                                        y_2=y_values[-2],
                                        x_2=x_values[-2])

                    delta_y = self.h * (1 / 12) * (5 * value(self.function, x, appr_y) +
                                              8 * value(self.function, x_values[-1], y_values[-1]) -
                                              value(self.function, x_values[-2], y_values[-2]))


                    if math.isnan(delta_y):
                        exit("Кажется, вы ввели невозможное условие или функцию")  # Bad manner

                    x += self.h  # update x
                    y += delta_y  # update y

                    y_values.append(y)
                    x_values.append(x)

                    if abs(y_values[-1] - y_values[-2]) > 1e+3:
                        print("Расходится\n")
                        break

            if order == 4:
                steps = order - 1

                # Calculate first 3 values using Runge-Kutta
                for i in range(steps + 1):
                    k1 = value(self.function, x, y)
                    k2 = value(self.function, x + self.h / 2, y + (self.h * k1) / 2)
                    k3 = value(self.function, x + self.h / 2, y + (self.h * k2) / 2)
                    k4 = value(self.function, x + self.h, y + self.h * k3)

                    delta_y = (1 / 6) * self.h * (k1 + 2 * k2 + 2 * k3 + k4)

                    y += delta_y
                    x += self.h

                    y_values.append(y)
                    x_values.append(x)

                for i in range(self.n_partitions - steps - 1):
                    appr_y = self.iterations3(
                                        x_n=x,

                                        y_1=y_values[-1],
                                        x_1=x_values[-1],

                                        y_2=y_values[-2],
                                        x_2=x_values[-2],

                                        y_3=y_values[-3],
                                        x_3=x_values[-3])

                    k1 = 9 * value(self.function, x, appr_y)
                    k2 = 19 * value(self.function, x_values[-1], y_values[-1])
                    k3 = -5 * value(self.function, x_values[-2], y_values[-2])
                    k4 =  value(self.function, x_values[-3], y_values[-3])


                    delta_y = (1 / 24) * self.h * (k1 + k2 + k3 + k4)

                    y += delta_y
                    x += self.h

                    y_values.append(y)
                    x_values.append(x)

            if order == 5:
                steps = order - 1

                # Calculate first 3 values using Runge-Kutta
                for i in range(steps + 1):
                    k1 = value(self.function, x, y)
                    k2 = value(self.function, x + self.h / 2, y + (self.h * k1) / 2)
                    k3 = value(self.function, x + self.h / 2, y + (self.h * k2) / 2)
                    k4 = value(self.function, x + self.h, y + self.h * k3)

                    delta_y = (1 / 6) * self.h * (k1 + 2 * k2 + 2 * k3 + k4)

                    y += delta_y
                    x += self.h

                    y_values.append(y)
                    x_values.append(x)

                for i in range(self.n_partitions - steps - 1):
                    appr_y = self.iterations4(
                        x_n=x,

                        y_1=y_values[-1],
                        x_1=x_values[-1],

                        y_2=y_values[-2],
                        x_2=x_values[-2],

                        y_3=y_values[-3],
                        x_3=x_values[-3],

                        y_4=y_values[-4],
                        x_4=x_values[-4])

                    k1 = 251 * value(self.function, x, appr_y)
                    k2 = 646 * value(self.function, x_values[-1], y_values[-1])
                    k3 = -264 * value(self.function, x_values[-2], y_values[-2])
                    k4 = 106 * value(self.function, x_values[-3], y_values[-3])
                    k5 = -19 * value(self.function, x_values[-4], y_values[-4])


                    delta_y = (1 / 720) * self.h * (k1 + k2 + k3 + k4 + k5)

                    y += delta_y
                    x += self.h

                    y_values.append(y)
                    x_values.append(x)




        if plot:
            plt.grid(True)
            plt.plot([(i + self.h) for i in range(self.n_partitions)],
                     y_values)
            plt.show()

        return y if not return_values else (y, y_values)