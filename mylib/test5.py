import numexpr as ne

with open("input.txt", "r") as file:

    f = str(file.readline()).strip("\n")

    y0 = float(next(file))

    interval = next(file)

    h = float(next(file))


interval = [float(i) for i in interval.split()]
a = interval[0]
b = interval[1]



def value(f, x=0.0, y=0.0):
    return ne.evaluate(f)

def euler(function, y0, a, b, h):

    x = a
    y = y0

    y_values = [y0]

    n_partitions = int((b - a) / h)

    print(n_partitions)

    for i in range(n_partitions):
        f = value(function, x, y)

        x += h
        y += h * f

        y_values.append(y)


    return (y, y_values)

# пример (чтобы поменять пример, надо поменять название входного файла)
result = euler(f, y0, a, b, h)

print(len(result[1]))


with open("output.txt", "w") as file:
    file.write(str(result[1]))


file.close()

import xlwt
def write_to_excel(filename, y):
    wb = xlwt.Workbook()

    sheet = wb.add_sheet("sheet1")

    # sheet.write(0, 0, "x")
    sheet.write(0, 1, "y")
    sheet.write(0, 2, "точный y")
    #
    # for i in range(len(x)):
    #     sheet.write(i+1, 0, str(x[i]))

    for i in range(len(y)):
        sheet.write(i+1, 1, str(y[i]))



    wb.save(filename)


write_to_excel("output_another.xls", result[1])
