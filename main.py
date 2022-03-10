from typing import Callable
import numpy as np
from matplotlib import pyplot as plt
from scipy import sparse
from scipy.sparse import linalg
from scipy import interpolate
import math
from time import time


def m_x_y_(x: float, y: float) -> float:
    """
    Функция пористости
    :param x: координата x
    :param y: координата y
    :return: значение пористости
    """
    return 0.2


def p_i(i: int, j: int) -> int:
    """
    индекс p[i][j] в матрице неизвестных
    :param i: индекс i
    :param j: индекс j
    :return: индекс в матрице неизвестных
    """
    return j + i * (N + 1)


def q0_i(i: int, j: int) -> int:
    """
    индекс q0[i][j] в матрице неизвестных
    :param i: индекс i
    :param j: индекс j
    :return: индекс в матрице неизвестных
    """
    return j + i * N + q_p


def q1_i(i: int, j: int) -> int:
    """
    индекс q1[i][j] в матрице неизвестных
    :param i: индекс i
    :param j: индекс j
    :return: индекс в матрице неизвестных
    """
    return j + i * M + q_p + q_q0


def cell_i(i: int, j: int) -> int:
    """
    индекс ячейки[i][j] в матрице неизвестных (N*M)
    :param i: индекс i
    :param j: индекс j
    :return: индекс в матрице неизвестных
    """
    return j + i * N


def gu_const_val(side: str, type_: str, value: float) -> None:
    """
    Задание граничного условия в виде константного значения функции
    :param side: сторона области (right, top, left, bottom)
    :param type_: тип ГУ: давление "p" или расход "q"
    :param value: значение
    """
    global gu
    for i in range(len(gu[side])):
        gu[side][i][0] = type_
        gu[side][i][1] = value


def gu_from_to_value(side: str, type_: str, value_from: float, value_to: float) -> None:
    """
    Задание граничного условия в виде линейной функции от первого значения до второго
    Значения задаются слева-направо либо снизу-вверх
    :param side: сторона области (right, top, left, bottom)
    :param type_: тип ГУ: давление "p" или расход "q"
    :param value_from: первое значение
    :param value_to: второе значение
    """
    global gu
    values = np.linspace(value_from, value_to, len(gu[side]))
    for i in range(len(gu[side])):
        gu[side][i][0] = type_
        gu[side][i][1] = values[i]


def gu_from_function(side: str, type_: str, func: Callable, xlim: list[float, float] = None) -> None:
    """
    Задание граничного условия в виде функции
    Значения задаются слева-направо либо снизу-вверх
    :param side: сторона области (right, top, left, bottom)
    :param type_: тип ГУ: давление "p" или расход "q"
    :param func: функция
    :param xlim: пределы изменения значения функции, которые будут интерполированы на узлы
    """
    if not xlim:
        if side == 'left' or side == 'right':
            xlim = [1, M - 1]
        elif side == 'top' or side == 'bottom':
            xlim = [1, N - 1]
    values = np.linspace(xlim[0], xlim[1], len(gu[side]))
    for i in range(values.size):
        gu[side][i][0] = type_
        gu[side][i][1] = func(values[i])


def add_poiseuille_equation_horizontal(i: int, j: int) -> None:
    """
    Добавление в матрицу уравнения Пуазейля для горизонтального канала
    :param i: индекс i горизонтального канала
    :param j: индекс j горизонтального канала
    """
    global data, row_ind, col_ind, eq_ind
    h = 0.5 * (b[i - 1][j - 1] + b[i][j - 1])
    coef = h ** 3 / 12
    data.extend([coef, -coef, -1])
    row_ind.extend([eq_ind, eq_ind, eq_ind])
    col_ind.extend([p_i(i + 1, j), p_i(i + 1, j + 1), q0_i(i, j)])
    eq_ind += 1


def add_poiseuille_equation_vertical(i: int, j: int) -> None:
    """
    Добавление в матрицу уравнения Пуазейля для вертикального канала
    :param i: индекс i вертикального канала
    :param j: индекс j вертикального канала
    """
    global data, row_ind, col_ind, eq_ind
    h = 0.5 * (b[j][i] + b[j][i + 1])
    coef = h ** 3 / 12
    data.extend([coef, -coef, -1])
    row_ind.extend([eq_ind, eq_ind, eq_ind])
    col_ind.extend([p_i(j, i + 1), p_i(j + 1, i + 1), q1_i(i, j)])
    eq_ind += 1


def add_balance_equation(i: int, j: int) -> None:
    """
    Добавление в матрицу уравнения баланса для внутреннего узла
    Суммарный расход через узел равен нулю
    :param i: индекс i узла
    :param j: индекс j узла
    """
    global data, row_ind, col_ind, eq_ind
    data.extend([1.0, 1.0, -1.0, -1.0])
    row_ind.extend([eq_ind, eq_ind, eq_ind, eq_ind])
    col_ind.extend([q0_i(i - 1, j), q1_i(j - 1, i), q0_i(i - 1, j - 1), q1_i(j - 1, i - 1)])
    eq_ind += 1


def add_value_p(i: int, j: int, value) -> None:
    """
    Добавление в систему уравнений конкретного значения для давления в узле
    :param i: индекс i узла
    :param j: индекс j узла
    :param value: значение давления
    """
    global data, row_ind, col_ind, eq_ind, rhs
    data.append(1)
    row_ind.append(eq_ind)
    col_ind.append(p_i(i, j))
    rhs[eq_ind] = value
    eq_ind += 1


def add_value_q_horizontal(i: int, j: int, value) -> None:
    """
    Добавление в систему уравнений конкретного значения для расхода в горизонтальном канале
    :param i: индекс i горизонтального канала
    :param j: индекс j горизонтального канала
    :param value: значение расхода
    """
    global data, row_ind, col_ind, eq_ind, rhs
    data.append(1)
    row_ind.append(eq_ind)
    col_ind.append(q0_i(i, j))
    rhs[eq_ind] = value
    eq_ind += 1


def add_value_q_vertical(i: int, j: int, value) -> None:
    """
    Добавление в систему уравнений конкретного значения для расхода в вертикальном канале
    :param i: индекс i вертикального канала
    :param j: индекс j вертикального канала
    :param value: значение расхода
    """
    global data, row_ind, col_ind, eq_ind, rhs
    data.append(1)
    row_ind.append(eq_ind)
    col_ind.append(q1_i(i, j))
    rhs[eq_ind] = value
    eq_ind += 1


def add_corner_p_equations() -> None:
    """
    Добавление в систему уравнений метода нахождения четырех угловых давлений
    Даления находятся как среднее арифметическое двух соседних давлений
    """
    global data, row_ind, col_ind, eq_ind
    data.extend([-2, 1, 1])
    row_ind.extend([eq_ind, eq_ind, eq_ind])
    col_ind.extend([p_i(0, 0), p_i(0, 1), p_i(1, 0)])
    eq_ind += 1
    data.extend([-2, 1, 1])
    row_ind.extend([eq_ind, eq_ind, eq_ind])
    col_ind.extend([p_i(0, N), p_i(1, N), p_i(0, N - 1)])
    eq_ind += 1
    data.extend([-2, 1, 1])
    row_ind.extend([eq_ind, eq_ind, eq_ind])
    col_ind.extend([p_i(M, N), p_i(M - 1, N), p_i(M, N - 1)])
    eq_ind += 1
    data.extend([-2, 1, 1])
    row_ind.extend([eq_ind, eq_ind, eq_ind])
    col_ind.extend([p_i(M, 0), p_i(M, 1), p_i(M - 1, 0)])
    eq_ind += 1


def extrap_func(x: list[float], y: list[float]) -> Callable:
    if len(x) == 1:
        return lambda xx: y[0]
    elif len(x) == 2:
        return interpolate.interp1d(x, y, kind='linear', fill_value='extrapolate')
    elif len(x) == 3:
        return interpolate.interp1d(x, y, kind='quadratic', fill_value='extrapolate')
    else:
        return interpolate.interp1d(x, y, kind='cubic', fill_value='extrapolate')


# Просто пример функции
def function(x):
    return x ** 2


global_start = time()
start = time()
# Количество ячеек по горизонтали и вертикали соответственно
N, M = 3, 4
# Количество ячеек
Q = N * M
# Количество уравнений
eq_Q = 3 * M * N + 1
# Количество неизвестных давлений
q_p = (N + 1) * (M + 1)
# Количество неизвестных Q0 (→)
q_q0 = (M - 1) * N
# Количество неизвестных Q1 (↑)
q_q1 = (N - 1) * M
# координаты центров ячеек
x = np.array([[j + 0.5 for j in range(N)] for i in range(M)])
y = np.array([[i + 0.5 for j in range(N)] for i in range(M)])
# Пористость в ячейках, ширины пор
m = np.array([[m_x_y_(x[i][j], y[i][j]) for j in range(N)] for i in range(M)])
b = np.array([[1 - math.sqrt(1 - m[i][j]) for j in range(N)] for i in range(M)])
# Значения, индексы по строке и столбцу для спарс матрицы
data, row_ind, col_ind = list(), list(), list()
rhs = np.zeros(eq_Q)
# Типы граничных условий на гранях (→↑←↓)
gu = {'right': [['p', 0] for _ in range(M - 1)],
      'top': [['p', 0] for _ in range(N - 1)],
      'left': [['p', 0] for _ in range(M - 1)],
      'bottom': [['p', 0] for _ in range(N - 1)]}
print(f'Система {N}x{M} ячеек')
print('=' * 40)
print(f'Инициализация: {time() - start} секунд')
start = time()
# Задание граничных условий
gu_const_val('top', 'q', 0)
gu_const_val('bottom', 'q', 0)
gu_const_val('left', 'p', 10)
gu_const_val('right', 'p', 0)
print(f'Задание граничных условий: {time() - start} секунд')

start = time()
# Формирование матрицы
eq_ind = 0

# Уравнения Пуазейля для горизонтальных каналов + ГУ
# Уравнения Пуазейля для вертикальных каналов + ГУ
# Уравнения баланса
# 4 угловых давления

for i in range(M - 1):
    if gu['left'][i][0] == 'q':
        add_value_q_horizontal(i, 0, gu['left'][i][1])
    else:
        add_value_p(i + 1, 0, gu['left'][i][1])

    for j in range(N):
        add_poiseuille_equation_horizontal(i, j)

    if gu['right'][i][0] == 'q':
        add_value_q_horizontal(i, N - 1, gu['right'][i][1])
    else:
        add_value_p(i + 1, N, gu['right'][i][1])

for i in range(N - 1):
    if gu['bottom'][i][0] == 'q':
        add_value_q_vertical(i, 0, gu['bottom'][i][1])
    else:
        add_value_p(0, i + 1, gu['bottom'][i][1])

    for j in range(M):
        add_poiseuille_equation_vertical(i, j)

    if gu['top'][i][0] == 'q':
        add_value_q_vertical(i, M - 1, gu['top'][i][1])
    else:
        add_value_p(M, i + 1, gu['top'][i][1])

for i in range(1, M):
    for j in range(1, N):
        add_balance_equation(i, j)

add_corner_p_equations()

sA = sparse.csc_matrix((tuple(data), (tuple(row_ind), tuple(col_ind))), shape=(eq_Q, eq_Q))
print(f'Формирование матрицы: {time() - start} секунд')

start = time()
# Решение системы
p_q0_q1 = linalg.spsolve(sA, rhs)
print(f'Решение системы: {time() - start} секунд')
print('=' * 40)

start = time()
# Расчет функции тока на смещенной сетке / Составление матрицы
row_ind = list()
col_ind = list()
data = list()
rhs = np.zeros(2 * M * N - M - N + 1)
eq_ind = 0
# горизонтальные расходы
for i in range(M):
    for j in range(N - 1):
        row_ind.extend([eq_ind, eq_ind])
        col_ind.extend([cell_i(i, j + 1), cell_i(i, j)])
        data.extend([1.0, -1.0])
        rhs[eq_ind] = p_q0_q1[q1_i(j, i)]
        eq_ind += 1
# вертикальные расходы
for j in range(N):
    for i in range(M - 1):
        row_ind.extend([eq_ind, eq_ind])
        col_ind.extend([cell_i(i + 1, j), cell_i(i, j)])
        data.extend([1.0, -1.0])
        rhs[eq_ind] = p_q0_q1[q0_i(i, j)]
        eq_ind += 1
# граничное условие
row_ind.append(eq_ind)
col_ind.append(cell_i(0, 0))
data.append(1.0)
sA = sparse.csc_matrix((tuple(data), (tuple(row_ind), tuple(col_ind))), shape=(rhs.size, N * M))
print(f'Формирование системы для функции тока: {time() - start} секунд')
start = time()
# Решение системы
ppsi = sparse.linalg.lsqr(sA, rhs)[0]
print(f'Решение системы: {time() - start} секунд')

start = time()
# Перерасчет psi на обычную сетку

# внутренняя сетка
psi = np.zeros((M + 1, N + 1))
for i in range(1, M):
    for j in range(1, N):
        psi[i][j] = 0.25 * (
                ppsi[cell_i(i - 1, j - 1)] + ppsi[cell_i(i, j - 1)] + ppsi[cell_i(i, j)] + ppsi[cell_i(i - 1, j)])
# экстраполяция
# вертикаль
for j in range(1, N):
    func = extrap_func(np.linspace(1, M - 1, M - 1), psi[:, j][1:-1])
    psi[0][j] = func(0)
    psi[M][j] = func(M)
# горизонталь
for i in range(1, M):
    func = extrap_func(np.linspace(1, N - 1, N - 1), psi[i][1:-1])
    psi[i][0] = func(0)
    psi[i][N] = func(N)
# углы
func_right = extrap_func(np.linspace(1, M - 1, M - 1), psi[:, N][1:-1])
func_left = extrap_func(np.linspace(1, M - 1, M - 1), psi[:, 0][1:-1])
func_top = extrap_func(np.linspace(1, N- 1, N - 1), psi[M][1:-1])
func_bottom = extrap_func(np.linspace(1, N - 1, N - 1), psi[0][1:-1])
psi[0][0] = (func_bottom(0) + func_left(0))*0.5
psi[0][N] = (func_bottom(N) + func_right(0))*0.5
psi[M][N] = (func_top(N) + func_right(M))*0.5
psi[M][0] = (func_top(0) + func_left(M))*0.5
print(f'Перерасчет psi на обычную сетку: {time() - start} секунд')
print('=' * 40)

start = time()
# Выделение величин
p = np.zeros((M + 1, N + 1))
for i in range(M + 1):
    for j in range(N + 1):
        p[i][j] = p_q0_q1[p_i(i, j)]
print(f'Выделение функции давления и функции тока из неизвестных величин: {time() - start} секунд')
print('=' * 40)


print(f'Время работы программы: {time() - global_start} секунд')
