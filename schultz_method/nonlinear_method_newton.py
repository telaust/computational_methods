from math import *
import string
import random
# from main import *
# from matrix import *
import numexpr as ne
import numpy as np
from mylib.matvec import *

""" Read system of equations """

equations = []

A0 = []

with open("input.txt", "r") as infile:

    infile.readline()  # pass line
    infile.readline()  # pass line

    while True:
        line = infile.readline().strip()
        if line != '':
            equations.append(line)
        else:
            infile.readline()  # pass line
            infile.readline()  # pass line

            init = infile.readline()

            infile.readline()  # pass line
            infile.readline()  # pass line
            infile.readline()  # pass line

            # Обработка матрицы А0
            # while True:
            #
            #     matrix_line = infile.readline().strip()
            #
            #     if matrix_line == '':
            #         break
            #     else:
            #         elements = [float(i) for i in matrix_line.split()]
            #         A0.append(elements)

            break


infile.close()



# Начальное приближение

X0 = [float(x) for x in init.split()]

# Уравнения

equations = [x.strip() for x in equations]

n_eq = len(equations)

# Проверка входных данных

if len(X0) != len(equations):
    # raise ValueError("Кол-во уравнений должно быть равно кол-ву неизвестных")
    exit("Кол-во уравнений должно быть равно кол-ву неизвестных")



""" Сколько ур-ний? Если больше 3, то используем a,b,c,d..., если меньше 3, то x, y, z """


def matnoiseadd(matr, noise='gauss'):

    check_matrix(matr)

    if noise is 'gauss':
        noise = random.gauss(0, 1)

    for i in range(len(matr)):
        for j in range(len(matr[0])):
            matr[i][j] += noise

    return matr


"""Если уравнения записаны в явном виде, то есть справа 0"""

equations = [x[:-4] for x in equations]


def vecnorm(a, norm="cheb"):

    max = 0

    for i in range(len(a)):
        if abs(a[i]) > max:
            max = abs(a[i])

    return max


def absvec(vec):

    for i in range(len(vec)):
        vec[i] = abs(vec[i])

    return vec

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

    # print("f1 ", f1)
    # print("f2 ", f2)

    # f1 = np.array(f1).tolist()
    # f2 = np.array(f2).tolist()

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


    # f1 = np.array(f1).tolist()
    # f2 = np.array(f2).tolist()

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

    # print("dict1 ", dict)

    f1 = ne.evaluate(func, local_dict=dict)

    dict2 = {}
    for i, key in enumerate(keys):
        dict2[key] = values[i]

    f2 = ne.evaluate(func, dict2)

    # print("dict2 ", dict2)

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


# это будем приближать
X = X0


# инициализация матрицы А0
# inv_jac = [
#     [0, 0.08],
#     [1, -0.08]
# ]

from mylib.matvec import *

# inv_jac = matinv(jacobian2(equations, X0))
#
# Ak = inv_jac


# pure_print(A0, "A0")

# Флаг несходимости
flag = False

# Допустимая погрешность
# eps = 1e-5

# Счетчик итераций
n_iters = 0

X_prev = X


# Какой размерности используем якобиан?

if n_eq == 2:
    jacobian = jacobian2  # alias for function
    valueof = valueof2


elif n_eq == 3:
    jacobian = jacobian3
    valueof = valueof3

else:
    jacobian = jacobianN
    valueof = valueofN

inv_jac = matinv(jacobian(equations, X0))

Ak = inv_jac

for k in range(100):

    pure_print(X, ("Xk = ", k ))

    FXk = valueof(equations, X)

    pure_print(FXk, "значения в этой точке")
    pure_print(Ak, ("А k= ", k))

    # Ak = round_matrix(Ak)

    pure_print(matvecmul(Ak, FXk), " MAT VEC MULTIPLICATION")

    X = vecvecsub(X,  matvecmul(Ak, FXk) )

    # X = round_vector(X)

    # Невязка
    PHI = matmatsub(identity_matrix(n_eq, n_eq), matmatmul(jacobian(equations, X), Ak) )

    Ak = matmatadd(Ak, matmatmul(Ak, PHI))

    # Считаем итерации для вывода их в файл
    n_iters += 1

    pure_print(X, "X")
    pure_print(X_prev, "x prev")

    # Условие остановки
    if vecnorm(vecvecsub(X, X_prev)) < 1e-9:
        break

    if vecnorm(vecvecsub(X, X_prev)) > 1e+3:
        flag = True  # этот флаг исп при выводе "Не сходится"
        break

    X_prev = X
    print('\n')


# Небольшое округление результата

for i in range(len(X)):
    X[i] = round_float(X[i], limit=5)


# Запись в файл

with open("output.txt", "w") as outfile:

    outfile.write("Решение (с округлением до 5 знаков после запятой):\n")
    outfile.write(str(X))

    outfile.write("\n\nЗначения функций в этих точках:\n")
    outfile.write(str(np.array(valueof(equations, X)).tolist() ))

    outfile.write("\n\nКол-во итераций, затраченных на поиск:\n")
    outfile.write(str(n_iters))

    if flag:
        outfile.write("\n\nНе сходится или плохо сходится при данных параметрах!")


outfile.close()


pure_print(jacobian2(equations, X0))

# inverse_jac = schultz_method(jacobian(equations, X0=X), U0=U0, clear=True)
    # pure_print(inverse_jac, "inverse jacobian")
    #
    # value_of_x = valueof(equations, X0=X)
    # pure_print(value_of_x, "value of X")
    #
    # X = vecvecsub(X, matvecmul(inverse_jac, value_of_x))  # формула Ньютона
    # pure_print(X, "x")

# f = "x**2+1"
#
# print(ne.evaluate(f, local_dict={'x':1}))

# pure_print(jacobian3(equations, X0))


pure_print( matinv(jacobian(equations, X0), clear=False), "inverse matrix")
# pure_print(jacobian(equations, X0), "jac1")

