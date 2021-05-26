from itertools import starmap
from operator import mul
import random
from main import *

def round_float(num, limit=2):
    return float(("%." + str(limit) + "f") % num)


def round_matrix(matr, limit=2):
    for i in range(len(matr)):
        for j in range(len(matr[0])):
            matr[i][j] = round_float(matr[i][j])

    return matr


def matvecmul(A, B):

    res = zerovector(len(B))
    for i in range(len(B)):
        res[i] = ( sum(A[i][j]*B[j] for j in range(len(A))) )
    return res


# нулевая матрица
def zeromatrix(rows, cols):
    return [[0 for x in range(cols)] for y in range(rows)]


def identity_matrix(rows, cols):
    if rows != cols:
        raise Exception("Identity matrix must be squared\n")

    I = zeromatrix(rows, cols)
    for i in range(len(I)):
        I[i][i] = 1

    return I


def gen_random_matrix(n_rows, n_cols, uniform="random", round_limit=None):
    if n_rows <= 0 or n_cols <= 0:
        raise Exception("{} and {} must be > 0\n".format(n_rows, n_cols))
    res = zeromatrix(n_rows, n_cols)

    if uniform is "random":

        for i in range(n_rows):
            for j in range(n_cols):
                res[i][j] = random.random()

    if uniform is "gauss":

        for i in range(n_rows):
            for j in range(n_cols):
                res[i][j] = random.gauss()

    else:
        raise Exception("{} uniform is not found, sorry\n".format(uniform))

    return res if round_limit is None else round_matrix(res, limit=round_limit)


def matminor(m, i, j):
    # here must be arg checking
    return [row[:j] + row[j + 1:] for row in (m[:i] + m[i + 1:])]


def det(matr):
    if type(matr[0]) is not list:
        raise Exception("{} must be repr. as list of lists\n".format(matr))

    if len(matr) == 2:
        return matr[0][0] * matr[1][1] - matr[0][1] * matr[1][0]

    determinant = 0

    for i in range(len(matr)):
        determinant += ((-1) ** i) * matr[0][i] * det(matminor(matr, 0, i))

    return determinant


""" метод Данилевского для нахождения собственных чисел матрицы """
def danilevsky(matr):

    if len(matr) != len(matr[0]):
        raise Exception("{} must be square matrix\n".format(matr))

    n = len(matr)

    for k in range(n-1):
        piv = matr[n-1-k][n-2-k]

        for i in range(n):
            matr[i][n-2-k] /= piv

        for i in range(n):
            for ii in range(n):
                if ii != n-2-k:
                    matr[i][ii] -= matr[i][n-2-k] * matr[n-1-k][ii]

    ######
    # pivot = matr[n-1][n-2]
    #
    # for i in range(n):
    #     matr[i][n-2] /= pivot
    #
    # for i in range(n):
    #     for ii in range(n):
    #         if ii != n-2:
    #             matr[i][ii] -= matr[i][n-2]*matr[n-1][ii]


    #######

    # pivot = matr[n - 2][n - 3]
    #
    # for i in range(n):
    #     matr[i][n - 3] /= pivot
    #
    # for i in range(n):
    #     for ii in range(n):
    #         if ii != n - 3:
    #             matr[i][ii] -= matr[i][n - 3] * matr[n - 2][ii]


    return matr

from main import pure_print

a = [
    [1.6, 2.3, 1.2, 4.3],
    [2.3, 0.6, 1.5, 1.1],
    [1.2, 1.5, 3.8, 0.8],
    [0.9, 1.1, 4.2, 6.3]
]






