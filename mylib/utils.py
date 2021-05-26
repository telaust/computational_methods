
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