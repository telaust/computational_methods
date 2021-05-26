from matvec import *
from numint import monte_carlo
from methods import delta2aitken
from numdiff import *



# print(derivative_n("(sin(6*x))**2", x=1, n=2))



def gauss_interp(x_array, y_array, x0, h):

    if len(x_array) != len(y_array):
        exit("not the same dimensions\n")

    differences = zeromatrix(len(x_array), len(y_array) - 2)

    for i in range(len(differences)):
        differences[i][0] = x_array[i]

    for i in range(len(differences)):
        differences[i][1] = y_array[i]

    # for i in range(len(differences)-1):
    #     for j in range(2, len(differences[0])):
    #         differences[i][j] = -differences[i][j-1] + differences[i+1][j-1]

    for j in range(2, len(differences[0])):
        for i in range(0, len(differences)-j+1):
            differences[i][j] = differences[i+1][j-1] - differences[i][j-1]

    differences = round_matrix(differences, limit=4)

    poly = str(y_array[int((len(y_array)-1)/2)])

    for k in range(3):
        poly += '+' + str(differences[3-k][2+k]) if differences[3-k][2+k] > 0 else str(differences[3-k][2+k])
        numerator = ...





    pure_print(differences)

    return poly

from methods import least_squares
x = [0.75, 1.5, 2.25, 3, 3.75]
y = [2.5, 1.2, 1.12, 2.25, 4.28]


# pure_print(gauss_interp(x, y, x0=0.35, h=0.05))
# print(gauss_interp(x, y, x0=0.35, h=0.05))
print(least_squares(x, y, m=2))

