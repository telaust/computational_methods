from math import *
import numexpr as ne



""" Reading input file / Считывание с файла """

def passnlines(filename, n):
    for _ in range(n):
        filename.readline()


with open("input.txt", "r") as file:

    try:
        passnlines(file, 4)

        f = str(file.readline()).strip("\n")

        passnlines(file, 3)

        y0 = float(next(file))

        passnlines(file, 3)

        interval = next(file)
        # t = float(next(file))
        # T = float(next(file))

        passnlines(file, 3)

        h = float(next(file))

    except Exception:
        exit("Incorrect input data\n")


interval = [float(i) for i in interval.split()]
t = interval[0]
T = interval[1]


def check_arguments(y0, t, T, h):

    # if y0 != t:
    #     raise Exception("y_0 ({}) doesn't match t ({})".format(y0, t))

    if t >= T:
        raise Exception("t = {} must be less than T = {}".format(t, T))

    if h >= T - t:
        raise Exception("step {} is too large or not to allow here".format(h))


def euler_method(f, y0, t, T, h, return_values=False):

    check_arguments(y0, t, T, h)

    n_partitions = int((T - t) / h)

    x = 0  # начальный икс всегда равен нулю в задаче Коши
    y = y0

    y_values = [y0]

    for i in range(n_partitions):

        fi = ne.evaluate(f)

        hf = h * fi

        y += hf
        x += h

        y_values.append(y)

    return y if not return_values else (y, y_values)


# Метод Адамса невяный 2-шаговый
def adams_explicit2(f, y_prev1, y_prev2, t, T, h, return_values=False):

    # ДОБАВИТЬ ПРОВЕРКУ АРГУМЕНТОВ

    n_partitions = int((T-t)/h)

    x = 0
    y_values = [y_prev2, y_prev1]

    y = 0
    for i in range(n_partitions):
        y = y_values[-1]  # define x
        f_prev1 = ne.evaluate(f)

        y = y_values[-2]  # define y
        f_prev2 = ne.evaluate(f)

        hf = h * ((3/2)*f_prev1 - (1/2)*f_prev2)

        y += hf
        x += h

        y_values.append(y)

    return y if not return_values else (y, y_values)



# eu_sol = euler_method(f, y0, t, T, h, True)
#
# print(eu)