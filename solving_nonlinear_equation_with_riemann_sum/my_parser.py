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

    print("length if str, is ", n)
    for i in range(n-3): #????
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


    for i in range(len(operators)):
        if len(operators) == 0:
            break

        print("operators ", len(operators))
        print("iter ", i)

        if i > len(operators)-1 or i < 0:
            break

        if operators[i] == "*" and len(operands) > 1:

            result = mul_lambdas(lexems[operands[i]],
                                 lexems[operands[i + 1]])
            del operands[i + 1]
            del operands[i]


            del operators[i]

        elif operators[i] == "*" and len(operands) == 1:

            result = mul_lambdas(result,
                                 lexems[operands[i]])
            del operands[i]

            del operators[i]


        elif operators[i] == "//" and len(operators) > 1:
            result = div_lambdas(lexems[operands[i]],
                                 lexems[operands[i + 1]])
            del operands[i + 1]
            del operands[i]


            del operators[i]

        elif operators[i] == "//" and len(operators) == 1:
            result = div_lambdas(result, lexems[operands[i]])
            del operands[i]

            del operators[i]


    # теперь цикл по самым неприоритетным операторам, если такие есть
    for i in range(len(operators)):
        if len(operators) == 0:
            return result

        if operators[i] == "+" and len(operators) > 1:
            result = sum_lambdas(lexems[operands[i]],
                                 lexems[operands[i+1]])

            del operands[i + 1]
            del operands[i]


            del operators[i]

        elif operators[i] == "+" and len(operands) == 1:
            result = sum_lambdas(result, lexems[operands[i]])
            del operands[i]

            del operators[i]

        # elif operators[i] == "-" and len(operands) > 1:


    return result


# for j in range(len(operators)):
#
#     if len(operators) == 0:
#         return res_lambda
#
#     print("iteration ", j, " len operators ", len(operators))
#     if operators[j] == '*' and len(operators) > 1:
#         # print(lexems[operands[j]])
#         res_lambda = mul_lambdas(lexems[operands[j]], lexems[operands[j+1]])
#         del operands[j]
#         # del operands[j+1]
#
#
#
#         del operators[j]
#
#
#
#     elif operators[j] == "*" and len(operators) == 1:
#         res_lambda = mul_lambdas(res_lambda, lexems[operands[j]])
#         del operands[j]
#
#         del operators[j]

# к этому моменту списки операторов и операндов должны быть удалены
# print(res_lambda)




# l = alternativa("cos(x)*cos(x)+sin(x)*sin(x)")
# plt.plot(space, l(space), color="g")
# plt.plot(space, np.cos(space)*np.cos(space)+np.sin(space)*np.sin(space), color="red")
#



# a = lambda x: np.sin(x)
# b = lambda x: np.cos(x)
# c = lambda x: a(x) * b(x)
#
#
# print(mul_lambdas(a, b) )
#
#
# plt.plot(space, mul_lambdas(a, b)(space), color="pink")
plt.show()
#
# def process_str(str):
#     list_str = list(str)
#     print(list_str)
#
#     a, b = [], []
#     cur_operator = []
#
#
#     # res_lambda
#
#     operand_stack = []
#
#     for i in range(len(list_str)):
#
#
#
#         if list_str[i] in signs:
#             pass
#
#         if list_str[i] == "*":
#             print("a >> ", a)
#
#             rr = ''.join(a)
#             func = lexems[rr]
#
#             if not operand_stack:
#                 operand_stack.append(func)
#             else:
#                 res_lambda = mul_lambdas(operand_stack[0], func)
#                 operand_stack.clear()
#
#
#             print(rr)
#             a.clear()
#             i += 2



        # if list_str[i] == "+" or list_str[i] == "*": # if + is met
        #     print(a)
        #
        #     res1 = ''.join(a) # concat chars into a single list
        #     print("concat into list: ", res1)
        #     lexem_lambda = lexems[res1]
            # cur_operator.append(list_str[i])


        # if list_str[i] == "*":
        #
        #     cur_operator.append(list_str[i])
        #     res = ''.join(a) # concat chars into a single list
        #     print(res)
        #     f_res = lexems[res] # take a value (lambda function) from corresponding key
        #     print(f_res)
        #     plt.plot(space, [f_res(x) for x in space])
        #     b = a
        #     print(b)
        #     a.clear() # clear a list
        #     i += 1

        # a.append(list_str[i])



# process_str(str)




# фунция, проверяющая на равенство количества скобок ()
def is_equally_parenthesis(str):
    parOpened = 0
    parClosed = 0

    for i in range(len(str)):
        if str[i] == '(':
            parOpened += 1
        elif str[i] == ')':
            parClosed += 1

    return parClosed == parOpened


lin = lambda x: x


str = "x*22.12"

num = 0

def is_number(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

for i in range(len(str)):
    if str[i] == '*' or str[i] == '+':
        if is_number(str[:i]):
            num = float(str[:i])
        elif is_number(str[i+1:]):
            num = float(str[i+1:])


def lin_func(a, b):
    return lambda x: x * a + b

def power_func(pow):
    return lambda x: x ** pow

def log():
    return lambda x: np.log(x)

# f1 = lin_func(num)

# эта ф-ия разносит по двум листам низкоприоритетные операции (+ или -)
#  и высокоприоритетные операнды

def distribute(str):
    operands = []
    operations = []


    for i in range(len(str)): # len(str) is immutable here
        if str[i] == "+" or str[i] == "-":
            if i == 0 or i == len(str)-1:
                raise Exception("illegal format")
            operations.append(str[i])

    operands = str.split("+")
    return operands, operations


# цель: проходя через все вложенности, получить окончательную лямбду
# и все манипуляции, которые делаются, делаются с лямбдой

# def eval_expr(str):
#     len_str = len(str)
#
#     for i in range(len_str):

str = ""
print()

# print(f)
# plt.plot(space, f1(space), color="g")
# plt.plot(space, 22.12*space, color="r")
# plt.show()