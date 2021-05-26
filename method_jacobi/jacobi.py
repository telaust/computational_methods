from mylib.matvec import *
import matplotlib.pyplot as plt


mat = []


def passnlines(filename, n):
    for _ in range(n):
        filename.readline()


with open("input.txt", "r") as infile:
    passnlines(infile, 2)

    while True:
        line = infile.readline().strip()

        if line != '':
            mat.append(line)
        else:
            passnlines(infile, 2)
            b = infile.readline()
            passnlines(infile, 3)
            eps = float(next(infile))
            break

infile.close()

b = [float(i) for i in b.split()]

A = [[float(i) for i in mat[i].split()] for i in range(len(mat))]

pure_print(A)
pure_print(b)

def mysum(X, exclude):
    s = 0
    for i in range(len(X)):
        if i == exclude:
            continue
        s += abs(X[i])

    return s


# деление вектора на диагональную матрицу !
def divondiag(X, D):
    if len(X) != len(D):
        raise Exception("разные размерности")

    for i in range(len(D)):
        if D[i] == 0:
            raise Exception("диагональная матрица не диагональна")

        X[i] = X[i] / D[i]

    return X

def converge(A):

    for i in range(len(A)):
        if abs(A[i][i]) <= mysum(A[i], i):
            return False

    return True

def norm(vec):

    vec = absvec(vec)
    max = 0

    for i in range(len(vec)):
        if vec[i] > max:
            max = vec[i]

    return max

x_history = []

max_iter = 0


def jacobi(A, b, eps=eps, n_iter=0, init='zeros'):

    if not converge(A):
        print("Matrix hasn't diagonal dominant, может не сойтись\n")
        exit()

    D = diag(A)
    R = diagzeros(A)

    if init is 'zeros':
        x = zerovector(len(A))
    elif init is 'ones':
        x = onesvector(len(A))
    else:
        raise ValueError("{} doesn't support\n".format(init))


    for i in range(n_iter):

        x_history.append(x)
        x = divondiag(vecvecsub(b, matvecmul(R, x)), D)

        if norm(vecvecsub(x, x_history[-1])) < eps:
            break

    return x



sol = jacobi(A, b, n_iter=25, init='zeros')
print("solution is ", sol)


with open("output.txt", "w") as file:

    file.write("Решение:\n")
    file.write(str(sol) )


file.close()

plt.title("Сходимость начальных точек к оптимальным, метод Якоби")
plt.ylabel("Решения")
plt.xlabel("Кол-во итераций")
plt.plot(range(len(x_history)), x_history)
plt.grid(True)
plt.show()
