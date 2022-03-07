import numpy as np
from matplotlib import pyplot as plt
from scipy import sparse
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


start = time()
# Количество ячеек по горизонтали и вертикали соответственно
N, M = 3, 2
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
# Пористость в ячейках, длины твердых включений, ширины пор
m = np.array([[m_x_y_(x[i][j], y[i][j]) for j in range(N)] for i in range(M)])
l = np.array([[math.sqrt(1 - m[i][j]) for j in range(N)] for i in range(M)])
b = np.array([[1 - l[i][j] for j in range(N)] for i in range(M)])
# Значения, индексы по строке и столбцу, форма для спарс матрицы
data, row_ind, col_ind = tuple(), tuple(), tuple()
shape = (eq_Q, eq_Q)

print(f'Инициализация: {time() - start} секунд')
