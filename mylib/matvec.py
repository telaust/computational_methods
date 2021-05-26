from itertools import starmap
from operator import mul
import random

""" Argument checking """


def check_vector(a):
    if type(a) is not list or type(a[0]) is list:
        raise ValueError("Type mismatch. {} must be repr. as list".format(a))


def check_matrix(A):
    if type(A) is not list or type(A[0]) is not list:
        raise ValueError("Type mismatch. {} must be repr. as list of lists".format(A))


""" Vector operations """

def vecnorm(a, norm="cheb"):

    max = 0

    for i in range(len(a)):
        if abs(a[i]) > max:
            max = abs(a[i])

    return max


def absvec(vec):
    check_vector(vec)

    res = zerovector(len(vec))

    for i in range(len(vec)):
        res[i] = abs(vec[i])

    return res



def round_vector(vec, limit=2):
    check_vector(vec)

    for i in range(len(vec)):
        vec[i] = round_float(vec[i], limit=limit)

    return vec


def vecvecdiv(a, b, eps=1e-5):
    check_vector(a)  # ASSERT ??
    check_vector(b)

    if len(a) != len(b):
        raise ValueError("Shape mismatch")

    for i in range(len(a)):
        if b[i] == 0:
            b[i] += eps
        a[i] = a[i] / b[i]

    return a


def vecvecmul(a, b):  # vector(transposed) by vector multiplication
    check_vector(a)
    check_vector(b)

    sum = 0

    if len(a) != len(b):
        raise Exception("not the same dimensions")

    for i in range(len(a)):
        sum += a[i] * b[i]

    return sum


def vecvecadd(a, b):
    check_vector(a)
    check_vector(b)

    res = zerovector(len(a))

    for i in range(len(res)):
        res[i] = a[i] + b[i]

    return res


def vecvecsub(a, b):
    check_vector(a)
    check_vector(b)

    res = zerovector(len(a))

    for i in range(len(res)):
        res[i] = a[i] - b[i]

    return res


def reverse_vector(vec):
    check_vector(vec)

    res = []
    for i in range(len(vec)):
        res.append(vec[-i-1])

    return res


def onesvector(size):
    if size <= 0:
        raise ValueError("size {} must be > 0".format(size))

    return [1 for x in range(size)]


def zerovector(size):
    if size <= 0:
        raise ValueError("{} must be > 0".format(size))

    res = []
    for _ in range(size):
        res.append(0)
    return res


def constvector(size, const=1):
    if size <= 0:
        raise ValueError("size {} must be > 0\n".format(size))

    return [const for x in range(size)]


def onehotvec(size, index):
    if size <= 0:
        raise ValueError("size {} must be > 0\n".format(size))

    res = zerovector(size)
    for i in range(size):
        if i == index:
            res[i] = 1
            break

    return res


# Sum of vector elements except one exluded

def sum_except(X, excluded):
    s = 0
    for i in range(len(X)):
        if i == excluded:
            continue
        s += abs(X[i])

    return s

""" Matrix operations """


# нулевая матрица
def zeromatrix(rows, cols):
    return [[0 for x in range(cols)] for y in range(rows)]


# проверка на диагональность матрицы
def isdiag(D):
    check_matrix(D)

    for i in range(len(D)):
        if D[i][i] == 0:
            return False
    return True


def matmatadd(A, B):
    check_vector(A)
    check_vector(B)

    res = zeromatrix(len(A), len(A[0]))

    for i in range(len(res)):
        for j in range(len(res[0])):
            res[i][j] = A[i][j] + B[i][j]

    return res


def matmatsub(A, B):
    check_matrix(A)
    check_vector(B)

    res = zeromatrix(len(A), len(A[0]))

    for i in range(len(res)):
        for j in range(len(res[0])):
            res[i][j] = A[i][j] - B[i][j]

    return res


def transpose(matrix):
    check_matrix(matrix)

    rows = len(matrix)
    cols = len(matrix[0])

    res = zeromatrix(cols, rows)

    for i in range(cols):
        for j in range(rows):
            res[i][j] = matrix[j][i]

    return res


# 2D matrices only
def dim(A):
    check_matrix(A)
    return [len(A), len(A[0])]


# matrix multiplication
def matmatmul(A, B):

    check_matrix(A)
    check_matrix(B)

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



def identity_matrix(rows, cols):
    if rows != cols:
        raise Exception("Identity matrix must be squared\n")

    I = zeromatrix(rows, cols)
    for i in range(len(I)):
        I[i][i] = 1

    return I


def onesmatrix(n_rows, n_cols):

    if n_rows <= 0 or n_cols <= 0:
        raise ValueError("size ({},{}) must be > 0".format(n_rows, n_cols))

    return [[1 for i in range(n_cols)] for j in range(n_rows)]


def constmatrix(n_rows, n_cols, const=1):

    if n_rows <= 0 or n_cols <= 0:
        raise ValueError("size ({},{}) must be > 0".format(n_rows, n_cols))

    return [[const for i in range(n_cols)] for j in range(n_rows)]


# Make diagonal elements of matrix zeros

def diagzeros(matr):

    check_matrix(matr)

    for i in range(len(matr)):
        for j in range(len(matr[0])):
            if i == j:
                matr[i][j] = 0

    return matr


def round_matrix(matr, limit=2):
    check_matrix(matr)

    for i in range(len(matr)):
        for j in range(len(matr[0])):
            matr[i][j] = round_float(matr[i][j], limit=limit)

    return matr


def matminor(m, i, j):
    check_matrix(m)

    assert (i >= 0 and j >= 0)

    return [row[:j] + row[j + 1:] for row in (m[:i] + m[i + 1:])]


def det(matr):
    check_matrix(matr)

    if len(matr) == 1:
        return matr[0][0]

    if len(matr) == 2:
        return matr[0][0] * matr[1][1] - matr[0][1] * matr[1][0]

    determinant = 0

    for i in range(len(matr)):
        determinant += ((-1) ** i) * matr[0][i] * det(matminor(matr, 0, i))

    return determinant


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


# Division matrix by scalar

def matnumdiv(matr, x, clear=False):
    check_matrix(matr)

    if x == 0:
        raise ValueError("Division be zero\n")

    for i in range(len(matr)):
        for j in range(len(matr[0])):
            print("matrix i j is >>> ", matr[i][j])
            matr[i][j] = matr[i][j] / x

    return round_matrix(matr, limit=3) if clear else matr


# Присоединенная матрица

def adjugate_matrix(matr):
    check_matrix(matr)

    n_rows = len(matr)
    n_cols = len(matr[0])

    res = zeromatrix(n_rows, n_cols)

    for i in range(n_rows):
        for j in range(n_cols):
            if (i+j) % 2 == 0:
                res[i][j] = det(matminor(matr, i, j))
            else:
                res[i][j] = (-1) * det(matminor(matr, i, j))

    return res


# Inverse matrix

def matinv(matr, clear=True):
    check_matrix(matr)

    return transpose(matnumdiv(adjugate_matrix(matr), det(matr), clear=clear))


# trace

def trace(matr):
    check_matrix(matr)

    return sum([matr[i][i] for i in range(len(matr))])

# Диагональное преобладание

def is_diagdominant(A):
    check_matrix(A)

    for i in range(len(A)):
        if abs(A[i][i]) <= sum_except(A[i], i):
            return False

    return True


def is_equalmat(A, B):
    check_matrix(A)
    check_matrix(B)

    for i in range(len(A)):
        for j in range(len(B)):
            if A[i][j] != B[i][j]:
                return False

    return True


def is_square(A):
    check_matrix(A)
    return len(A) == len(A[0])


# is equal diagonal elements of matrix?

def is_eqdiagmat(A):
    check_matrix(A)

    for i in range(len(A)-1):
        if A[i][i] != A[i+1][i+1]:
            return False
    return True


# is matrix symmetrical ?

def is_symmat(A):
    check_matrix(A)
    return is_equalmat(A, transpose(A)) and is_eqdiagmat(A)


# is lower triangular matrix ?

def is_lowermat(A):
    check_matrix(A)

    for i in range(len(A)):
        for j in range(len(A[0])):
            if j >= i and A[i][j] != 0:
                return False
            if j < i and A[i][j] == 0:
                return False
    return True

# is upper triangular matrix ?

def is_uppermat(A):
    check_matrix(A)

    for i in range(len(A)):
        for j in range(len(A[0])):
            if j < i and A[i][j] != 0:
                return False
            if j >= i and A[i][j] == 0:
                return False
    return True


""" Vector-Matrix operations """


def matvecmul(A, B):
    check_matrix(A)
    check_vector(B)

    res = zerovector(len(B))

    for i in range(len(B)):
        res[i] = sum(A[i][j] * B[j] for j in range(len(A)))
    return res


# возвращает вектор, состоящий из диагональных элементов данной матрицы

def diag(X):
    check_matrix(X)

    res = []

    for i in range(len(X)):
        res.append(X[i][i])

    return res


def pure_print(matr, name="Matrix", clear=False):

    if type(matr[0]) is not list:  # if a 1D vector
        print('-' * 5, name, '-' * 5)
        print(str(matr).replace('[', '(').replace(']', ')'))

    else: # if a 2D matrix
        print('-' * 5, name, '-' * 5, '\n')
        for row in range(len(matr)):
            print('\t', end='')
            for col in range(len(matr[0])):

                print( "{}  ".format(
                    (matr[row][col]) if not clear
                    else (round_float(matr[row][col])) ), end='')

            print("\n")
        print('-'*(12+len(name)))  # 10 + 2 spaces around a name


def round_float(num, limit=2):
    if limit <= 0:
        raise ValueError("limit = {} must be > 0".format(limit))

    return float(("%." + str(limit) + "f") % num)


def matnoiseadd(matr, noise='gauss'):

    check_matrix(matr)

    if noise is 'gauss':
        noise = random.gauss(0, 1)

    for i in range(len(matr)):
        for j in range(len(matr[0])):
            matr[i][j] += noise

    return matr


def matrix_norm(matr, p="cheb", clear=True):

    check_matrix(matr)

    n_rows = len(matr)
    n_cols = len(matr[0])

    res_norm = None

    if p is not None and p.isdigit(): # эта норма захватывает Фробениуса(евклидова) при p=2
        p = float(p) # is p was a string
        s = 0

        for i in range(n_rows):
            for j in range(n_cols):
                s += pow(abs(matr[i][j]), p)

        res_norm = pow(s, 1/p)


    elif p is "maxnorm" or p is "max":
        res_norm = matr[0][0]

        for i in range(n_rows):
            for j in range(n_cols):
                if matr[i][j] > res_norm:
                    res_norm = matr[i][j]


    elif p is "one":

        res_norm = 0

        for i in range(n_cols):
            s = 0
            for j in range(n_rows):
                s += abs(matr[i][j])
            if s > res_norm:
                res_norm = s

    elif p is "cheb":

        res_norm = 0

        for i in range(n_rows):
            s = 0
            for j in range(n_cols):
                s += abs(matr[i][j])
            if s > res_norm:
                res_norm = s


    else:
        raise Exception("I cannot realize this norm, sorry\n")

    return round_float(res_norm) if clear else res_norm
