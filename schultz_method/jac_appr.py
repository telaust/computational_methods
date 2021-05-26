import numexpr as ne
from mylib.matvec import *
from calculus import *

""" Read system of equations """

equations = []

A0 = []

with open("input.txt", "r") as infile:

    infile.readline()  # pass line
    infile.readline()  # pass line

    while True:
        line = infile.readline().strip()
        if line != '':
            equations.append(line)
        else:
            infile.readline()  # pass line
            infile.readline()  # pass line

            init = infile.readline()

            infile.readline()  # pass line
            infile.readline()  # pass line
            infile.readline()  # pass line

            # Обработка матрицы А0
            # while True:
            #
            #     matrix_line = infile.readline().strip()
            #
            #     if matrix_line == '':
            #         break
            #     else:
            #         elements = [float(i) for i in matrix_line.split()]
            #         A0.append(elements)

            break


infile.close()

# Начальное приближение

X0 = [float(x) for x in init.split()]

# Уравнения

equations = [x.strip() for x in equations]

n_eq = len(equations)

# Проверка входных данных

if len(X0) != len(equations):
    # raise ValueError("Кол-во уравнений должно быть равно кол-ву неизвестных")
    exit("Кол-во уравнений должно быть равно кол-ву неизвестных")

"""Если уравнения записаны в явном виде, то есть справа 0"""

equations = [x[:-4] for x in equations]

# это будем приближать
X = X0

# Флаг несходимости
flag = False

# Допустимая погрешность
# eps = 1e-5

# Счетчик итераций
n_iters = 0

X_prev = X

# Какой размерности используем якобиан?

if n_eq == 2:
    jacobian = jacobian2  # alias for function
    valueof = valueof2
elif n_eq == 3:
    jacobian = jacobian3
    valueof = valueof3
else:
    jacobian = jacobianN
    valueof = valueofN

Ak = matinv(jacobian(equations, X0))

for k in range(100):

    FXk = valueof(equations, X)

    X = vecvecsub(X,  matvecmul(Ak, FXk) )

    # Невязка
    PHI = matmatsub(identity_matrix(n_eq, n_eq), matmatmul(jacobian(equations, X), Ak) )

    Ak = matmatadd(Ak, matmatmul(Ak, PHI))

    # Считаем итерации для вывода их в файл
    n_iters += 1

    # Условие остановки
    if vecnorm(vecvecsub(X, X_prev)) < 1e-9:
        break

    if vecnorm(vecvecsub(X, X_prev)) > 1e+3:
        flag = True  # этот флаг исп при выводе "Не сходится"
        break

    X_prev = X
    print('\n')


# Небольшое округление результата

for i in range(len(X)):
    X[i] = round_float(X[i], limit=5)

# Запись в файл

with open("output.txt", "w") as outfile:

    outfile.write("Решение (с округлением до 5 знаков после запятой):\n")
    outfile.write(str(X))

    outfile.write("\n\nЗначения функций в этих точках:\n")
    outfile.write(str(np.array(valueof(equations, X)).tolist() ))

    outfile.write("\n\nКол-во итераций, затраченных на поиск:\n")
    outfile.write(str(n_iters))

    if flag:
        outfile.write("\n\nНе сходится или плохо сходится при данных параметрах!")


outfile.close()
