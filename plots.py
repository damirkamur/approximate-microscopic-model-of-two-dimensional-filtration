from matplotlib import pyplot as plt
import numpy as np

N, M = 400, 200
X, Y = np.meshgrid(np.linspace(0, N + 1, N + 1), np.linspace(0, M + 1, M + 1))
X1, Y1 = np.meshgrid(np.linspace(0.5, N - 0.5, N), np.linspace(0.5, M - 0.5, M))
psi = np.zeros((M + 1, N + 1))
p = np.zeros((M + 1, N + 1))
ppsi = np.zeros((M, N))
m = np.zeros((M, N))
with open(f'results/psi{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        psi[i, :] = list(map(float, line.split()))
with open(f'results/p{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        p[i, :] = list(map(float, line.split()))

with open(f'results/m{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        m[i, :] = list(map(float, line.split()))

try:
    p_an = np.zeros((M + 1, N + 1))
    with open(f'results/analit_p{N}x{M}.dat', 'r') as file:
        for i, line in enumerate(file):
            p_an[i, :] = list(map(float, line.split()))
except Exception as _ex:
    p_an = None

plt.figure(1)
plt.axes().set_aspect('equal')
plt.contourf(X, Y, psi, np.linspace(np.min(psi), np.max(psi), 100))
plt.colorbar()
plt.title('Функция тока')
plt.show()

plt.figure(2)
plt.axes().set_aspect('equal')
plt.contourf(X, Y, p, np.linspace(np.min(p), np.max(p), 100))
plt.colorbar()
plt.title('Функция давления')
plt.show()

plt.figure(3)
plt.axes().set_aspect('equal')
plt.contourf(X1, Y1, m, np.linspace(np.min(m), np.max(m), 100))
plt.colorbar()
plt.title('пористость')
plt.show()

if not p_an is None:
    plt.figure(3)
    plt.contourf(p_an, np.linspace(np.min(p_an), np.max(p_an), 30))
    plt.colorbar()
    plt.title('Функция давления аналитическая')
    plt.show()

    dif_p = np.zeros((M + 1, N + 1))
    maximum_p = np.max(p_an)
    for i in range(M + 1):
        for j in range(N + 1):
            dif_p[i][j] = abs((p_an[i][j] - p[i][j]) / maximum_p)

    plt.figure(4)
    plt.contourf(dif_p)
    plt.colorbar()
    plt.title('Относительная погрешность давления')
    plt.show()

    print(f'Максимальная пгрешность: {np.max(dif_p)}')
