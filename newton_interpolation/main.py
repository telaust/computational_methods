from mylib.matvec import *
from math import factorial

x = [0, 0.2, 0.4, 0.6, 0.8, 1]

y = [1.272, 4.465, 5.644, 5.809, 3.961, 2.101]


def newton_interp(x, y, x_pivot):

    if len(x) != len(y):
        raise ValueError("{} and {} haven't the same dims".format(x, y))

    # теперь нужно еще len(x)-1 векторов, где будем хранить разд разности
    # пусть это будет матрица..

    differences = zeromatrix(len(x)-1, len(x)-1)

    for i in range(len(differences)):
        differences[i][0] = y[i+1] - y[i]

    for j in range(1, len(differences[0])):
        for i in range(0, len(differences)-j):
            differences[i][j] = differences[i+1][j-1] - differences[i][j-1]

    t = (x_pivot - x[0]) / (x[1] - x[0])

    sum = 0

    for i in range(1, len(x)):
        piece = 1

        for j in range(i):
            piece *= (t-j)

        piece /= factorial(i)

        sum += differences[0][i-1] * piece

    return sum + y[0]



print(newton_interp(x, y, 0.1))


