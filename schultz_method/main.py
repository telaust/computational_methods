from itertools import starmap
from operator import mul


def matvecmul(A, b):
    return [sum(starmap(mul, zip(b, col))) for col in zip(*A)]


# нулевая матрица
def zeromatrix(rows, cols):
    return [[0 for x in range(cols)] for y in range(rows)]


def zerovector(size):
    if size <= 0:
        raise Exception("size {} must be > 0\n".format(size))

    return [0 for x in range(size)]


def onesvector(size):
    if size <= 0:
        raise Exception("size {} must be > 0\n".format(size))

    return [1 for x in range(size)]


def constvector(size, const=1):
    if size <= 0:
        raise Exception("size {} must be > 0\n".format(size))

    return [const for x in range(size)]


def onehotvec(size, index):
    if size <= 0:
        raise Exception("size {} must be > 0\n".format(size))

    res = zerovector(size)
    for i in range(size):
        if i == index:
            res[i] = 1
            break

    return res


def identity_matrix(rows, cols):
    if rows != cols:
        raise Exception("Identity matrix must be squared\n")

    I = zeromatrix(rows, cols)
    for i in range(len(I)):
        I[i][i] = 1

    return I


# 2D matrices only
def dim(A):
    return [len(A), len(A[0])]

def vecvecadd(a, b):
    if type(a) is not list or type(b) is not list:
        raise Exception("{} and {} must be represented as lists".format(a, b))
    if len(a) != len(b):
        raise Exception("{} and {} must have the same size".format(a, b))

    res = zerovector(len(a))

    for i in range(len(res)):
        res[i] = a[i] + b[i]

    return res


def vecvecsub(a, b):
    if type(a) is not list or type(b) is not list:
        raise Exception("{} and {} must be represented as lists".format(a, b))

    if len(a) != len(b):
        raise Exception("{} and {} must have the same size".format(a, b))

    res = zerovector(len(a))

    for i in range(len(res)):
        res[i] = a[i] - b[i]

    return res



# matrix multiplication
def matmatmul(A, B):
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


def matmatadd(A, B):
    if dim(A) != dim(B):
        raise Exception("different size of {} and {}".format(A, B))

    res = zeromatrix(len(A), len(A[0]))

    for i in range(len(res)):
        for j in range(len(res[0])):
            res[i][j] = A[i][j] + B[i][j]

    return res

def matmatsub(A, B):
    if dim(A) != dim(B):
        raise Exception("different size of {} and {}".format(A, B))

    res = zeromatrix(len(A), len(A[0]))

    for i in range(len(res)):
        for j in range(len(res[0])):
            res[i][j] = A[i][j] - B[i][j]

    return res


def pure_print(matr, name="Matrix", clear=False):

    if type(matr[0]) is not list: # if a 1D vector
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
    return float(("%." + str(limit) + "f") % num)


def round_matrix(matr, limit=2):
    for i in range(len(matr)):
        for j in range(len(matr[0])):
            matr[i][j] = round_float(matr[i][j])

    return matr


def matrix_norm(matr, p="cheb", clear=True):

    if type(matr) is not list or type(matr[0]) is not list:
        raise Exception("{} must be matrix (repr. as list of lists)\n".format(matr))

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


def vecnorm(a, norm="cheb"):

    max = 0

    for i in range(len(a)):
        if abs(a[i]) > max:
            max = abs(a[i])

    return max


A = [
    [1, 2],
    [3, 4]
]


U0 = [
    [0.6, -0.5],
    [0.1, 0.6]
]


def is_square(matr):
    if len(matr) != len(matr[0]):
        return False

    return True


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
