import numexpr as ne
from mylib.matvec import *
import numpy as np
import string

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