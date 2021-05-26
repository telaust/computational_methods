import numexpr as ne
from matvec import pure_print
from numdiff import derivativeX, derivativeY
f = 'x + 2*y'
g = 'x**2 + y**2 - 5'


def brawn(f, g, X0, eps=1e-7):

    p, q = [1e+5], [1e+5]
    xk, yk = 1, 3
    x, y = [], []

    while max(p[-1], q[-1]) > eps:

        if derivativeX(f, [xk, yk]) == 0:
            xk_temp = xk - ne.evaluate(f, {'x': xk, 'y': yk}) / ( derivativeX(f, [xk, yk]) + 1e-7)

        xk_temp = xk - ne.evaluate(f, {'x': xk, 'y': yk}) / derivativeX(f, [xk, yk])
        qk = ne.evaluate(g, {'x': xk_temp, 'y': yk} ) * derivativeX(f, [xk, yk]) / (derivativeX(f, [xk, yk])*derivativeY(g, [xk, yk]) - derivativeY(f, [xk, yk])*derivativeX(f, [xk_temp, yk]))

        pk = (ne.evaluate(f, {'x': xk, 'y': yk}) - qk * derivativeY(f, [xk, yk])) / derivativeX(f, [xk, yk])

        p.append(pk)
        q.append(qk)

        xk = xk - pk
        yk = yk - qk

        x.append(xk)
        y.append(yk)


    return [x, y]


pure_print(brawn(f, g, [0, 0]))