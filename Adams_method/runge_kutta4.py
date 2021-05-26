from math import *
import matplotlib.pyplot as plt
import numpy as np
import numexpr as ne

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

print("logging : f = {}, y0 = {}, t = {}, T = {}, h = {}".format(f, y0, t, T, h))


def value(f, x=0.0, y=0.0): return ne.evaluate(f)


def runge_kutta4(f, y0, t, T, h, return_values=False):
    x = 0
    y = y0

    n_partitions = int((T - t) / h)

    y_values = [y0]

    for i in range(n_partitions):
        k1 = value(f, x, y)
        k2 = value(f, x + h / 2, y + (h * k1) / 2)
        k3 = value(f, x + h / 2, y + (h * k2) / 2)
        k4 = value(f, x + h, y + h * k3)

        delta_y = (1 / 6) * h * (k1 + 2 * k2 + 2 * k3 + k4)

        x += h
        y = y + delta_y

        print("i=",i, "\tdy=", delta_y)

        y_values.append(y)

    return y if not return_values else (y, y_values)


y, vals = runge_kutta4(f, y0, t, T, h, True)

print(vals)

# print(vals)
#
# plt.plot(vals, label="Method Runge-Kutta")
# plt.legend()
# plt.grid(True)
# plt.suptitle("Methods")
# plt.xlim([0, 10])
# plt.ylim([0, 1])

# plt.show()