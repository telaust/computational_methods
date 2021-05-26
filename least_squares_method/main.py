from itertools import starmap
from operator import mul
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
from gaussian_elimination import gauss



def pure_print(matr, name="Matrix"):
    print('-'*5, name, '-'*5, '\n')
    for row in range(len(matr)):
        print('\t', end='')
        for col in range(len(matr[0])):

            print( "{}  ".format(matr[row][col]), end='')

        print("\n")
    print('-'*(12+len(name))) # 10 + 2 spaces around a name


def vecvecmul(a, b): # vector(transposed) by vector multiplication
    sum = 0

    if len(a) != len(b):
        raise Exception("not the same dimensions")

    for i in range(len(a)):
        sum += a[i] * b[i]

    return sum


# нулевая матрица
def zeromatrix(rows, cols):
    return [[0 for x in range(cols)] for y in range(rows)]


# 2D matrices only
def dim(A):
    return [len(A), len(A[0])]


# транспонирование матрицы - работает правильно
def transpose(matrix):
    rows = len(matrix)
    cols = len(matrix[0])


    res = zeromatrix(cols, rows)

    for i in range(cols):
        for j in range(rows):
            res[i][j] = matrix[j][i]

    return res

# matrix multiplication
def matmul(A, B):
    # dimension checking
    if dim(A) == dim(B):
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



def matvectormul(A, b):
    return [sum(starmap(mul, zip(b, col))) for col in zip(*A)]


# Эта функция сделана так, как описывается метод наим квадратов в Турчаке
def least_squares(x_data, y_data, m):
    n = len(x_data)

    phi = [[0] * (m+1) for i in range(n)]

    for i in range(len(x_data)):
        for j in range(m+1):
            phi[i][j] = x_data[i] ** j

    transposed_phi = transpose(phi)

    FtF = matmul(transposed_phi, phi)

    Fty = matvectormul(phi, y_data)

    params = gauss(FtF, Fty)

    return params


# read data
with open("input.txt", "r") as infile:
    lines = [line.split() for line in infile]

infile.close()

matr = zeromatrix(len(lines), len(lines[0]))

for i in range(len(lines)):
    for j in range(len(lines[0])):
        matr[i][j] = float(lines[i][j])


x_vec = []

for i in range(len(matr)):
    x_vec.append(matr[i][0])

y_vec = []

for i in range(len(matr)):
    y_vec.append(matr[i][1])

mycoefs = least_squares(x_vec, y_vec, 2)
print(mycoefs)
plt.scatter(y_vec, x_vec)

def plot_smooth_graph(x, y, k, color='black', legend='Legend'):
    x_new = np.linspace(min(x), max(y), len(y))
    spline = make_interp_spline(x_new, y, k=k)
    y_smooth = spline(x_new)

    plt.plot(x_new, y_smooth, color=color, label=legend)
    plt.legend(loc="upper right")
    plt.show()





# return string representation of polynomial function
def make_polynomial(coefs):
    degree = len(coefs)

    print(degree)

    str_func = ""
    for i in range(degree):
        if coefs[i] > 0:
            str_coef = '+' + str(coefs[i])
        elif coefs[i] < 0:
            str_coef = str(coefs[i])
        else:
            str_coef = '0'

        str_func += str_coef + "*x**" + str(i)

    return str_func


def plot_string(str_repr, x_range):
    x = np.array(x_range)
    y = eval(str_repr)
    plt.plot(x, y)
    plt.show()


print(make_polynomial(mycoefs))
plot_string(make_polynomial(mycoefs), x_vec)
