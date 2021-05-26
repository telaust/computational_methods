import numpy as np
import matplotlib.pyplot as plt
import time


start_time = time.time()



with open("input.txt", 'r') as f:

    filestr = f.read()
    tokens = filestr.split() # strings
    tokens = list(map(int, tokens)) # tokens (integers)
    print(tokens) # ints


    N = tokens[0]
    a = []
    coeflist = tokens[1:-N]
    print(coeflist)

    for i in range(N):
        a.append(coeflist[0: N])
        coeflist = coeflist[N:]

    b = tokens[-N:]


f.close()

def mysum(X, exclude):
    s = 0
    for i in range(len(X)):
        if i == exclude:
            continue
        s += abs(X[i])

    return s

def isconverge(A):

    for i in range(len(A)):
        if abs(A[i][i]) <= mysum(A[i], i):
            return False

    return True

def myzeros(n): # создает нулевой вектор
    res = []
    for i in range(n):
        res.append(0)

    return res

# возвращает вектор, состоящий из диагональных элементов данной матрицы WORK
def mydiag(X):
    res = []

    for i in range(len(X)):
        res.append(X[i][i])

    return res


def pureprint(matrix): # pure print
    for i in range(len(matrix)):
        print(matrix[i], "\n")


# произведение матрицы на вектор, в результате вектор
def mydot(A, B):

    res = myzeros(len(B))
    for i in range(len(B)):
        res[i] = ( sum(A[i][j]*B[j] for j in range(len(A))) )
    return res



def zeromatrix(n): # создает нулевую матрицу
    res = []
    for i in range(n):
        res.append([0] * n)
    return res



def mydif(A, B): # разность двух матриц
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        raise Exception("разные размерности")

    res = zeromatrix(len(A))

    for i in range(len(A)):
        for j in range(len(B)):
            res[i][j] = A[i][j] - B[i][j]

    return res

# проверка на диагональность матрицы
def isdiag(D):
    for i in range(len(D)):
        if D[i][i] == 0:
            return False
    return True

# деление вектора на диагональную матрицу !
def divondiag(X, D):
    if len(X) != len(D):
        raise Exception("разные размерности")

    for i in range(len(D)):
        if D[i] == 0:
            raise Exception("диагональная матрица не диагональна")

        X[i] = X[i] / D[i]

    return X




def mydiagMatrix(X): # создает диагональную матрицу из данной

    res = zeromatrix(len(X))

    for i in range(len(X)):
        for j in range(len(X[0])):
            if i != j:
                res[i][j] = 0

    return res

def dif2vectors(a, b): # difference of 2 vectors
    res = myzeros(len(a))
    for i in range(len(a)):
        res[i] = a[i] - b[i]
    return res


    # return [map(lambda x, y: x-y, ii, jj) for ii, jj in zip(a,b)]


# сделать диагональные элементы матрицы нулями
def diagzeros(A):

    if type(A) is not list or type(A[0]) is not list:
        raise Exception("{} must be repr. as list of lists\n".format(A))

    for i in range(len(A)):
        for j in range(len(A[0])):
            if i == j:
                A[i][j] = 0

    return A


def normalvector(size): # return list
    return list(np.random.normal(size=size))

def absvec(vec):
    if type(vec) is not list:
        raise Exception("{} must be list".format(vec))

    res = myzeros(len(vec))

    for i in range(len(vec)):
        res[i] = abs(vec[i])

    return res

def norm(vec):

    vec = absvec(vec)
    max = 0

    for i in range(len(vec)):
        if vec[i] > max:
            max = vec[i]

    return max


h = []
ind = []
def jacobi(A, b,eps=1e-6, N=100, init=None):

    if not isconverge(a):
        raise Exception("не сходится")

    # создать вектор, состоящий из только диагональных элементов
    D = mydiag(A)

    # R - это матрицы А, но с нулями вместо диагональных элементов
    R = diagzeros(A)

    if init == 'zeros' or init is None:
        x = myzeros(len(A))
    else:
        x = init

    x_prev = x

    for i in range(N):
        h.append(x)
        ind.append(i)
        x = divondiag(dif2vectors(b, mydot(R, x)), D)

        if norm(dif2vectors(x, x_prev)) < eps:
            break

        x_prev = x

    # чертим график сходимости скаляров вектора
    plt.title("Сходимость начальных значений к правильным, метод Якоби")
    plt.plot(ind, h)
    plt.xlabel("номер итарации")
    plt.ylabel("значения вектора")
    plt.grid(True)
    plt.show()


    return x

sol = jacobi(a, b, N=25, init=[1, 1])
print("solution is ", sol)



with open("out.txt", 'w') as outfile:
    outfile.write(str(sol))

outfile.close()


print("\nTime you spent on that task is {:.2f} sec".format(time.time() - start_time))


