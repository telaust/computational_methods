from mylib.matvec import *

def vecnorm(a, norm="cheb"):

    max = 0

    for i in range(len(a)):
        if abs(a[i]) > max:
            max = abs(a[i])

    return max


def check_args(A, b):
    check_matrix(A)
    check_vector(b)

    if len(A) != len(b):
        exit("{} and {} must have the same dimensions".format(A, b))

    if not is_diagdominant(A):
        print("Matrix {} hasn't diagonal dominant\nIt couldn't converge\n".format(A))



def numvecmul(l, a):
    check_vector(a)

    res = a

    for i in range(len(a)):
        res[i] *= l

    return res


def richardson(A, b, x0=None, omega=1.0, eps=1e-5):
    check_args(A, b)

    if x0 is None:
        x0 = zerovector(len(b))

    x = x0

    residual = vecnorm(vecvecsub(matvecmul(A, x), b))

    while residual > eps:

        addend = numvecmul(omega, vecvecsub(b, matvecmul(A, x)))

        x = vecvecadd(x, addend)

        residual = vecnorm(vecvecsub(matvecmul(A, x), b))

        if residual > 1e+2:
            exit("It doesnt converge. Change omega")


    return x


if __name__ == '__main__':
    A = [
        [10, 1, 1],
        [2, 10, 1],
        [12, 2, 10]
    ]

    b = [12, 13, 24]
    guess = [0, 0, 0]

    print(richardson(A, b, omega=1))