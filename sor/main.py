import cmath
from mylib.matvec import *

import numpy as np

def is_conj(z1, z2):
    return complex.conjugate(z1) == z2


def is_hermitian(matr):

    for i in range(len(matr)):
        for j in range(i):
            if not is_conj(matr[i][j], matr[j][i]):
                return False
    return True


def test_comlex():

    a1 = [
        [5, complex(1, 2)],
        [complex(1, -2), 7]
    ]
    a2 = [
        [5, complex(1, 5)],
        [complex(1, -2), 7]
    ]

    print(is_hermitian(a2))


def numvecmul(l, a):
    check_vector(a)

    res = a

    for i in range(len(a)):
        res[i] *= l

    return res



def vecnorm(a, norm="cheb"):

    max = 0

    for i in range(len(a)):
        if abs(a[i]) > max:
            max = abs(a[i])

    return max


def check_args(A, b):
    check_matrix(A)
    check_vector(b)

    if len(A) != len(b):
        exit("{} and {} must have the same dimensions".format(A, b))

    if not is_diagdominant(A):
        print("Matrix {} hasn't diagonal dominant\nIt couldn't converge\n".format(A))


# Метод релаксации

def sor(A, b, omega, initial_guess=None, eps=1e-5):

    check_args(A, b)

    if initial_guess is None:
        initial_guess = zerovector(len(b))

    x = initial_guess


    residual = vecnorm(vecvecsub(matvecmul(A, x), b))

    while residual > eps:

        for i in range(len(A)):
            sigma = 0

            for j in range(len(A[0])):

                if j != i:

                    sigma += A[i][j] * x[j]

            x[i] = (1 - omega) * x[i] + (omega / A[i][i]) * (b[i] - sigma)

        residual = vecnorm(vecvecsub(matvecmul(A, x), b))

        if residual > 1e+2:
            exit("It doesnt converge")

    return x


if __name__ == '__main__':

    A = [
        [10, 1, 1],
        [2, 10, 1],
        [12, 2, 10]
    ]

    b = [12, 13, 24]
    guess = [0, 0, 0]



    print(sor(A, b, omega=0.1, initial_guess=guess, eps=1e-3))



