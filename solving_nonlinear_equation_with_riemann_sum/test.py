import numpy as np
import matplotlib.pyplot as plt

lexems = {

    "cos(x)": lambda x: np.cos(x),
    "sin(x)": lambda x: np.sin(x),
    "tan(x)": lambda x: np.tan(x),
    "arcsin(x)": lambda x: np.arcsin(x),
    "arccos(x)": lambda x: np.arccos(x),
    "sinh(x)": lambda x: np.sinh(x),
    "cosh(x)": lambda x: np.cosh(x), # hyperbolic cos
    "tanh(x)": lambda x: np.tanh(x), # hyperbolic sin
    "arcsinh(x)": lambda x: np.arcsinh(x), # hyperbolic arcsin
    "arccosh(x)": lambda x: np.arccosh(x), # hyperbolic arccos
    "exp(x)": lambda x: np.exp(x),
    "log(x)": lambda x: np.log(x),
    "log10(x)": lambda x: np.log10(x),
    "log2(x)": lambda x: np.log2(x),
    # ...
    "sqrt(x)": lambda x: np.sqrt(x), # square root
    "cbrt(x)": lambda x: np.cbrt(x), # cube root

}


def mul_lambdas(lambda1, lambda2):
    return lambda x: lambda1(x) * lambda2(x)


def sum_lambdas(lambda1, lambda2):
    return lambda x: lambda1(x) + lambda2(x)


def sub_lambdas(lambda1, lambda2):
    return lambda x: lambda1(x) - lambda2(x)


def div_lambdas(lambda1, lambda2):
    # проверка на 0 не нужна, тк функции
    return lambda x: lambda1(x) / lambda2(x)


signs = ['+', '-', "\\", "*"]

str = "cos(x)+sin(x)*sin(x)"

print(str)

space = np.linspace(1, 100, 1000)


def alternativa(str):

    list_str = list(str)
    n = len(list_str)
    operators = []

    for i in range(n-2): #????
        if list_str[i] in signs:
            operators.append(list_str[i])
            del list_str[i]

    print("resulting list ", list_str)
    print("ops ", operators)

    operands = []

    p = len(list_str)

    temp = []

    for i in range(p):
        temp.append(list_str[i])

        if list_str[i] == ')':
            operands.append(''.join(temp))
            temp.clear()


    op_len = len(operators)


    print("op len is " , len(operators))

    # сначала цикл по самым приоритетным операторам
    # в конце цикла они должны удалиться из списка операторов, если были там
    for sign in range(len(operators)):



        # if len(operators) == 0:
        #     break

        if sign == "*" and len(operators) > 1:
            res_lambda = mul_lambdas(lexems[operands[sign]],
                                     lexems[operands[sign + 1]])
            del operands[sign]
            del operands[sign + 1]

            del operators[sign]

        elif operators[sign] == "*" and len(operators) == 1:
            res_lambda = mul_lambdas(res_lambda,
                                     lexems[operands[sign]])
            del operands[sign]

            del operators[sign]


        elif sign == "//" and len(operators) > 1:
            res_lambda = div_lambdas(lexems[operands[sign]],
                                     lexems[operands[sign + 1]])
            del operands[sign]
            del operands[sign + 1]

            del operators[sign]

        elif operators[sign] == "//" and len(operators) == 1:
            res_lambda = div_lambdas(res_lambda, lexems[operands[sign]])
            del operands[sign]

            del operators[sign]


    # теперь цикл по самым неприоритетным операторам, если такие есть
    for sign in range(len(operators)):
        if len(operators) == 0:
            return res_lambda

        if sign == "+" and len(operators) > 1:
            res_lambda = sum_lambdas(lexems[operands[sign]],
                                     lexems[operands[sign+1]])

            del operands[sign]
            del operands[sign+1]

            del operators[sign]

        elif operators[sign] == "+" and len(operators) == 1:
            res_lambda = div_lambdas(res_lambda, lexems[operands[sign]])
            del operands[sign]

            del operators[sign]
    return res_lambda


alternativa(str)
