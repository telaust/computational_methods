from numint import monte_carlo
from numdiff import DiffSolver
import math
f = '- (x*y) / (x + 1)'
y0 = 1
a = 0
b = 3
h = 0.01

# method = DiffSolver(f, y0, a, b, h)

# print(method.adams(order=3, explicit=False))
# print("метод Р-К 4 порядка ", method.runge_kutta(order=4, plot=False))
# print(method.runge_kutta(order=3, plot=True))
# print(method.runge_kutta(order=2, plot=False))
# print(method.runge_kutta(order=1, plot=False))

# print("Адамс, явный, 5 порядок ", method.adams(order=5, explicit=True, return_values=False, plot=False))
# print("Адамс, явный, 2 порядок ", method.adams(order=2, explicit=True, return_values=False, plot=False))
#
# print('\n')
#
# print("Адамс, неявный, 3 порядок ", method.adams(order=3, explicit=False, return_values=False, plot=False))
# print("Адамс, неявный, 2 порядок ", method.adams(order=2, explicit=False, return_values=False, plot=False))
# print("Адамс, неявный, 4 порядок ", method.adams(order=4, explicit=False, return_values=False, plot=False))
# print("Адамс, неявный, 5 порядок ", method.adams(order=5, explicit=False, return_values=False, plot=False))

# from numint import *
# res = riemann_sum("sin(x)", 0, 2, 120)
# print(res)

from matvec import *


a = [
    [1, 3],
    [2, 4],
    [5, 6]
]

b = [1, 2]


print(matvecmul(transpose(a), b))


from methods import least_squares
x_data = [0.75, 1.5, 2.25, 3, 3.75]
y_data = [2.5, 1.2, 1.12, 2.25, 4.28]
m=2


n = len(x_data)

phi = [[0] * (m+1) for i in range(n)]

for i in range(len(x_data)):
    for j in range(m+1):
        phi[i][j] = x_data[i] ** j


transposed_phi = transpose(phi)

FtF = matmatmul(transposed_phi, phi)

A = [
    [1.0,  0.75,  0.5625],
    [1.0,  1.5,  2.25],
    [1.0,  2.25 , 5.0625 ] ,
    [1,  3,  9]  ,
    [1.0,  3.75,  14.0625]]


pure_print(phi, "Фт")
pure_print(y_data, "y data")

Fty = matvecmul( transpose(A), y_data)

params = gauss(FtF, Fty)

# print(least_squares(x, y, m=2))