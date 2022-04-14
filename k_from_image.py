from typing import Callable

import numpy as np
from PIL import Image
from scipy import interpolate
import math
from scipy.misc import derivative
from time import time


def m_from_image(image_name: str, min_m: float = 0.0, max_m: float = 1.0) -> None:
    """
    Функция для генерации пористости по картинке.
    Пористость интерполируется от min_m (для самого темного оттенка) до max_m(для самого светлого)
    :param image_name: название файла картинки
    :param min_m: минимальное значение пористости
    :param max_m: масимальное значение пористости
    """

    def interpolate(color):
        return min_m + (max_m - min_m) * (color - minimum_color) / (maximum_color - minimum_color)

    global m, N, M
    m = np.zeros((M, N))
    image = Image.open(image_name).convert('L')
    image = image.resize((N, M))
    pixels = image.load()
    colors = np.array([[pixels[x, y] for x in range(N)] for y in range(M - 1, -1, -1)])
    maximum_color = np.max(colors)
    minimum_color = np.min(colors)
    m = np.array([[interpolate(colors[i][j]) for j in range(N)] for i in range(M)])


def extrap_func(x: list[float], y: list[float]) -> Callable:
    if len(x) == 1:
        return lambda xx: y[0]
    elif len(x) == 2:
        return interpolate.interp1d(x, y, kind='linear', fill_value='extrapolate')
    elif len(x) == 3:
        return interpolate.interp1d(x, y, kind='quadratic', fill_value='extrapolate')
    else:
        return interpolate.interp1d(x, y, kind='cubic', fill_value='extrapolate')


def k0_i(i: int, j: int) -> int:
    return j + i * N


def k1_i(i: int, j: int) -> int:
    return j + i * M


N = 200
M = 100
print(f'Система {N}x{M} ячеек')
print('=' * 40)
start = time()
m_from_image('field.png', 0.2, 0.5)
b = np.array([[1 - math.sqrt(1 - m[i][j]) for j in range(N)] for i in range(M)])
# Количество k0 (→)
q_k0 = (M - 1) * N
# Количество k1 (↑)
q_k1 = (N - 1) * M
k_0 = np.array([(b[i][j] * 0.5 + b[i + 1][j] * 0.5) ** 3 / 12 for i in range(M - 1) for j in range(N)])
k_1 = np.array([(b[i][j] * 0.5 + b[i][j + 1] * 0.5) ** 3 / 12 for j in range(N - 1) for i in range(M)])
print(f'Инициализация: {time() - start} секунд')

k1 = np.zeros((M, N))
k2 = np.zeros((M, N))

dk_dx = np.zeros((M, N))
dk_dy = np.zeros((M, N))

start = time()
# горизонтальная аппроксимация
for i in range(M):
    x = np.linspace(1, N - 1, N - 1)
    y = np.array([k_1[k1_i(j, i)] for j in range(N - 1)])
    f = extrap_func(x, y)
    x1 = np.linspace(0.5, N - 0.5, N)
    y1 = f(x1)
    for j in range(N):
        k1[i][j] = y1[j]
        dk_dx[i][j] = derivative(f, x1[j], 0.0001)
# вертикальная аппроксимация
for j in range(N):
    x = np.linspace(1, M - 1, M - 1)
    y = np.array([k_0[k0_i(i, j)] for i in range(M - 1)])
    f = extrap_func(x, y)
    x1 = np.linspace(0.5, M - 0.5, M)
    y1 = f(x1)
    for i in range(M):
        k2[i][j] = y1[i]
        dk_dy[i][j] = derivative(f, x1[i], 0.0001)
# общая аппроксимация
k = 0.5 * (k1 + k2)
print(f'Аппроксимация: {time() - start} секунд')

start = time()
with open(f'k_from_image/k{N}x{M}.dat', 'w') as file:
    for i in range(M):
        for j in range(N):
            file.write('%15.6f' % (k[i][j]) + '\t' * 2)
        file.write('\n')
print(f'Выгрузка k: {time() - start} секунд')

start = time()
with open(f'k_from_image/dk_dx{N}x{M}.dat', 'w') as file:
    for i in range(M):
        for j in range(N):
            file.write('%15.6f' % (dk_dx[i][j]) + '\t' * 2)
        file.write('\n')
print(f'Выгрузка dk_dx: {time() - start} секунд')

start = time()
with open(f'k_from_image/dk_dy{N}x{M}.dat', 'w') as file:
    for i in range(M):
        for j in range(N):
            file.write('%15.6f' % (dk_dy[i][j]) + '\t' * 2)
        file.write('\n')
print(f'Выгрузка dk_dy: {time() - start} секунд')