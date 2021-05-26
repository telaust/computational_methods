from math import *
import matplotlib.pyplot as plt

# метод Эйлера
# на вход подается производная
# начальное условие Коши
# отрезок

with open("input.txt", "r") as file:
    f = str(file.readline()).strip("\n")

    a = file.readline() # pass line

    y0 = float(next(file))

    b = file.readline() # pass line

    # t = float(next(file)) # start of interval
    # T = float(next(file)) # end of interval

    interval = next(file)
    t = float(interval[0])
    T = float(interval[2])


    c = file.readline() #

    h = float(next(file))


print("f = {}, y0 = {}, t = {}, T = {}, h = {}".format(f, y0, t, T, h))


def check_arguments(y0, t, T, h):

    # if y0 != t:
    #     raise Exception("y_0 ({}) doesn't match t ({})".format(y0, t))

    if t >= T:
        raise Exception("t = {} should be less than T = {}".format(t, T))

    if h >= T - t:
        raise Exception("step {} is too large or not to allow here", format(h))


def euler_method(f, y0, t, T, h, return_values=False):

    check_arguments(y0, t, T, h)

    n_partitions = int((T - t) / h)

    x = 0
    y = y0

    y_values = [y0]

    for i in range(n_partitions):

        fi = eval(f)

        hf = h * fi

        y += hf
        x += h

        y_values.append(y)

    return y if not return_values else (y, y_values)


y, values = euler_method(f, y0, t, T, h, True)

print(y)


plt.plot(values)
plt.show()