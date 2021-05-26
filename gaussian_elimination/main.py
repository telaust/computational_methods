# Gauss elimination method for solving linear equations Ax = b

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


def zero_vector(size):
    return [0 for i in range(size)]


# gauss elimination main function; returns vector
def gauss(A, b):
    lenA = len(A)

    for i in range(lenA):
        if is_close_to_zero(A[i][i]):
            print("matrix is seems to be singular")

            if i == 0:
                A = change_rows(A, i, lenA-1)
            if i == lenA-1:
                A = change_rows(A, i, 0)

    x = zero_vector(lenA)

    for k in range(lenA-1):
        for i in range(k+1, lenA):
            temp = A[i][k] / A[k][k]
            b[i] = b[i] - temp * b[k]

            for j in range(k+1, lenA):
                A[i][j] = A[i][j] - temp * A[k][j]

    x[lenA-1] = b[lenA-1] / A[lenA-1][lenA-1]

    for k in reversed(range(0, lenA-1)):

        s = 0
        for p in range(k+1, lenA):
            s = s + A[k][p] * x[p]

        x[k] = (b[k] - s) / A[k][k]

    return x


A = [
    [2, 1, 4, 5],
    [7, 9, 3, 6],
    [9, 3, 1, 1],
    [8, 8, 3, 6]
]

b = [11, 13, 43, 23]

res = gauss(A, b)

print(res)