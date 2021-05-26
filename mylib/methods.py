""" Здесь будут методы для решения различных задач из линейной алгебры """
from utils import *
from matvec import *
import numexpr as ne

# gauss elimination main function; returns vector
def gauss(A, b, clear=False):

    if not is_square(A):
        raise ValueError("{} must be square matrix".format(A))

    if len(A) != len(b):
        raise ValueError("Dimension mismatch")


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

    return x if not clear else round_vector(x)


# Эта функция сделана так,
# как описывается метод наим квадратов в Турчаке
def least_squares(x_data, y_data, m):
    n = len(x_data)

    phi = [[0] * (m+1) for i in range(n)]


    for i in range(len(x_data)):
        for j in range(m+1):
            phi[i][j] = x_data[i] ** j

    pure_print(phi)

    transposed_phi = transpose(phi)

    FtF = matmatmul(transposed_phi, phi)
    pure_print(phi)
    print('\n')
    pure_print(y_data)

    Fty = matvecmul( phi , y_data)

    params = gauss(FtF, Fty)

    return params


def value_of_polynom(coefs, x_value):
    n = len(coefs)
    first = []
    for i in range(n):
        first.append(x_value ** i)

    return vecvecmul(first, coefs)


def local_smooth(x, y, pivot, k, m):  # здесь pivot - это индекс опорной точки

    """

    k - number of data points left and right of pivot point (k/2 on the left, k/2 on the right)
    m - degree of approximate polynom

    How to use:
    in your code write something like that:

    res = [] # smoothed values will be here

    for i in range(len(x_noise)):  # loop over given data
        res.append(local_smooth(list(x_noise), list(y_noise), pivot=i, k=k_file, m=m_file))

    """

    if pivot == 0:
        return y[0]

    if (pivot - k/2) <= 0:
        return y[pivot]

    if (pivot + k/2) >= len(y)-1:
        return y[pivot]

    if pivot >= len(x) or pivot < 0:
        raise Exception("Out of data bounds")

    ranged_y = y[pivot - int(k/2) - 1 : pivot  + int(k/2)]
    ranged_x = x[pivot - int(k/2) - 1 : pivot  + int(k/2)]

    params = least_squares(ranged_x, ranged_y, m)

    y_hat = value_of_polynom(params, x[pivot])

    return y_hat


# TODO: доделать

def delta2aitken(func, x0, q, eps):

    x1 = ne.evaluate(func, {'x': x0})
    x2 = ne.evaluate(func, {'x': x1})

    for i in range(100):

        p = (x2 - 2*x1 + x0)
        if p == 0:
            p += 1e-7

        delta2 = (x0 * x2 - x1**2) / p

        x3 = ne.evaluate(func, {'x': delta2})

        if abs(x3 - delta2) < eps:
            break

        x0 = delta2
        x1 = x3
        x2 = ne.evaluate(func, {'x': x1})

    return x3


# Method for matrix inverse

def schultz_method(A, U0, eps=1e-5, m=1, n_iterations=20, clear=False):

    if not is_square(A):
        raise Exception("Only square matrix can be inverse")

    rows = len(A)
    cols = len(A[0])

    U = U0

    for i in range(n_iterations):
        k = i
        PHI = matmatsub(identity_matrix(rows, cols), matmatmul(A, U))
        if matrix_norm(PHI) <= eps:
            break  # we find inverse matrix
        U = matmatmul(U, matmatadd(identity_matrix(rows, cols), PHI))

        if matrix_norm(PHI) >= 1:
            break

    return U if not clear else round_matrix(U)