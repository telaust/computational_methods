from mylib.matvec import *

A = [
    [100, 2, 3],
    [4, 500, 6],
    [7, 8, 100]
]



def is_square(A):
    check_matrix(A)
    return len(A) == len(A[0])


def summ(a, b, start, end):

    check_vector(a)
    check_vector(b)

    if end <= start:
        return 0

    res = zerovector(end - start)

    for i in range(start, end):
        res[i] = a[i]*b[i]

    return sum(res)


def lu_decomp(a, clear=False):

    if not is_square(a):
        exit("LU decomposition can be applied to square matrix only")


    check_matrix(a)

    n_rows = len(a)
    n_cols = len(a[0])

    res = zeromatrix(n_rows, n_cols)

    for i in range(n_rows):

        for j in range(n_cols):

            if i == 1:
                res[i][j] = a[i][j] # first row of result matrix

            if i <= j:
                res[i][j] = a[i][j] - summ(res[i], res[j], start=0, end=i)

            if i > j:
                res[i][j] = (1 / res[j][j]) * (a[i][j] - summ(res[i], res[j], start=0, end=j))

    return res if not clear else round_matrix(res)


pure_print(lu_decomp(A), clear=False)



