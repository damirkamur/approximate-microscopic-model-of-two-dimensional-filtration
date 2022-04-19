from matplotlib import pyplot as plt
import numpy as np


def psi_i(i: int, j: int) -> int:
    return j + i * (N1 + 1)


def psi_i_i(k: int) -> tuple[int, int]:
    return k // (N1 + 1), k % (N1 + 1)


N, M = 1000, 500
X, Y = np.meshgrid(np.linspace(0, N + 1, N + 1), np.linspace(0, M + 1, M + 1))

N1, M1 = 500, 250
X1, Y1 = np.meshgrid(np.linspace(0, N + 1, N1 + 1), np.linspace(0, M + 1, M1 + 1))

Xm, Ym = np.meshgrid(np.linspace(0.5, N - 0.5, N), np.linspace(0.5, M - 0.5, M))
m = np.zeros((M, N))
with open(f'results/m{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        m[i, :] = list(map(float, line.split()))

psi = np.zeros((M + 1, N + 1))
with open(f'results/psi{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        psi[i, :] = list(map(float, line.split()))

v_psi = np.zeros((M1 + 1, N1 + 1))
with open('res140_70.dat', 'r') as file:
    line = file.readline()  # TITLE
    line = file.readline()  # VARIABLES
    psi_ind = 2
    line = file.readline()  # ZONE
    for k, line in enumerate(file):
        i, j = psi_i_i(k)
        v_psi[i][j] = list(map(float, line.split()))[psi_ind]

min_psi = np.min(psi)
max_psi = np.max(psi)
min_v_psi = np.min(v_psi)
max_v_psi = np.max(v_psi)

v_psi -= min_v_psi
v_psi *= (max_psi - min_psi) / (max_v_psi - min_v_psi)
v_psi += min_psi

plt.axes().set_aspect('equal')
CS0 = plt.contourf(Xm, Ym, m, np.linspace(np.min(m), np.max(m), 100))

CS1 = plt.contour(X, Y, psi, np.linspace(np.min(psi), np.max(psi), 15), colors='red')
plt.clabel(CS1, fontsize=7)
plt.colorbar()
CS2 = plt.contour(X1, Y1, v_psi, np.linspace(np.min(v_psi), np.max(v_psi), 15), colors='blue')
plt.clabel(CS2, fontsize=7)
plt.colorbar()
plt.show()
