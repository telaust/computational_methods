from math import *
from scipy.misc import *

f = lambda x: x**2 + 5*x
d = derivative(f, 1.0, dx=1e-6)

print(d)


n = 3 # сколько переменных и функций

# значит, надо создать якобиан размерности n*n
#
# def jacobian(funclist, x0):
#     if len(funclist) != n:
#         raise Exception("not the same dimensions\n")
#
#     matrix = []
#
#     for j in range(len(funclist)):
#         row = []
#
#         for i in range(n):
#             row.append(funclist[j])
#
#




a = [
    [1, 2],
    [3, 4]
]

def det2(matrix2x2):
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]

print(det2(a))