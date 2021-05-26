from itertools import starmap
from operator import mul
import math

def is_close_to_zero(x, eps=1e-6):
    return -eps < x < eps

# меняет местами строки
def change_rows(matrix, a, b):
    if a >= len(matrix) or a < 0 or b >= len(matrix) or b < 0:
        raise Exception("out of matrix bounds")

    temp = matrix[b]
    matrix[b] = matrix[a]
    matrix[a] = temp

    return matrix


# нулевая матрица
def zeromatrix(rows, cols):
    return [[0 for x in range(cols)] for y in range(rows)]

# нулевой вектор
def zerovector(size):
    if size <= 0:
        raise ValueError("{} must be > 0".format(size))

    res = []
    for _ in range(size):
        res.append(0)
    return res

# Решение СЛАУ методом Гаусса
def gauss(A, b):

    lenA = len(A)

    for i in range(lenA):
        if is_close_to_zero(A[i][i]):
            # print("matrix is seems to be singular")

            if i == 0:
                A = change_rows(A, i, lenA-1)
            if i == lenA-1:
                A = change_rows(A, i, 0)

    x = zerovector(lenA)

    for k in range(lenA-1):
        for i in range(k+1, lenA):
            if A[k][k] == 0:
                A[k][k] += 1e-5

            temp = A[i][k] / A[k][k]
            b[i] = b[i] - temp * b[k]

            for j in range(k+1, lenA):
                A[i][j] = A[i][j] - temp * A[k][j]

    if A[lenA-1][lenA-1] == 0:
        A[lenA - 1][lenA - 1] += 1e-5

    x[lenA-1] = b[lenA-1] / A[lenA-1][lenA-1]

    for k in reversed(range(0, lenA-1)):

        s = 0
        for p in range(k+1, lenA):
            s = s + A[k][p] * x[p]

        x[k] = (b[k] - s) / A[k][k]

    return x


# Транспонирование матрицы
def transpose(matrix):

    rows = len(matrix)
    cols = len(matrix[0])

    res = zeromatrix(cols, rows)

    for i in range(cols):
        for j in range(rows):
            res[i][j] = matrix[j][i]

    return res

# Умножение матрицы на веткор
def matmatmul(A, B):

    if len(A) == len(B):
        res = zeromatrix(len(A), len(A))
    else:
        res = zeromatrix(len(A), len(B[0]))

    for i in range(len(A)):
        for j in range(len(B[0])):
            sum = 0
            for k in range(len(A[0])):
               sum += A[i][k] * B[k][j]

            res[i][j] = sum

    return res


# Умножение матрицы на вектор
def multiply(first, second):
    return [sum(starmap(mul, zip(first, col))) for col in zip(*second)]

# пример 2
# x_data = [0.75, 1.5, 2.25, 3, 3.75]
# y_data = [2.5, 1.2, 1.12, 2.25, 4.28]
m = 2

# пример 1
# x_data = [0.75, 1.5, 2.25, 3, 3.75]
# y_data = [math.sin(2*x) for x in x_data ]
m = 2

# пример 3
x_data = [i for i in range(1, 21)]
y_data = [math.log(x) for x in x_data]
m = 3

n = len(x_data)

phi = [[0] * (m+1) for i in range(n)]

for i in range(len(x_data)):
    for j in range(m+1):
        phi[i][j] = x_data[i] ** j


transposed_phi = transpose(phi)

FtF = matmatmul(transposed_phi, phi)


Fty = multiply(y_data, phi)

params = gauss(FtF, Fty)

print(params)