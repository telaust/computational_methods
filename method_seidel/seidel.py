from mylib.matvec import *
import numpy as np


A = [
    [10, 1, 1],
    [2, 10, 1],
    [12, 2, 10]
]

b = [12, 13, 24]

x0 = [0, 0, 0]



A2 = [[1, -1],
      [2, 1]
]

b2 = [5, 7]


# матрица итерирования / матрица перехода

def matrix_B(a):
    check_matrix(a)
    res = zeromatrix(len(a), len(a[0]))

    for i in range(len(res)):

        for j in range(len(res[0])):
            if i == j:
                res[i][j] = 0
            else:
                res[i][j] = - a[i][j]/ a[i][i]

    return res


def vector_c(A, b):
    check_matrix(A)
    check_vector(b)

    res = zerovector(len(b))

    for i in range(len(res)):
        res[i] = b[i] / A[i][i]

    return res



def norm(a):
    abs_vec = absvec(a)
    max = 0
    for i in range(len(a)):
        if a[i] > max:
            max = a[i]

    return max


def vecnorm(a, norm="cheb"):

    max = 0

    for i in range(len(a)):
        if abs(a[i]) > max:
            max = abs(a[i])

    return max


def seidel(A, b, x0, max_iter=10, eps=1e-4):

    if not is_diagdominant(A):
        print("Matrix hasn't diagoanl dominant in this case")

    check_matrix(A)
    check_vector(b)

    x = x0

    C = vector_c(A, b)
    B = matrix_B(A)  # new b

    for k in range(max_iter):

        x_prev = x
        for i in range(len(A)):

            x[i] = vecvecmul(B[i], x) + C[i]

        if vecnorm(x, x_prev) < eps:
            break

        if vecnorm(x, x_prev) > 1e+3:
            print(">"*10  + "\tIt doesn't converge\n")
            exit()

    return x


print(seidel(A, b, [0, 0, 0]))