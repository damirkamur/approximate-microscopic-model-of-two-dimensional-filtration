from matplotlib import pyplot as plt
import numpy as np

M = 100
N = 2 * M
M2 = 10
N2 = 2 * M2
X, Y = np.meshgrid(np.linspace(0, N2, N + 1), np.linspace(0, M2, M + 1))
X1, Y1 = np.meshgrid(np.linspace(0, N2, 2 * N2 + 1), np.linspace(0, M2, 2 * M2 + 1))
X2, Y2 = np.meshgrid(np.linspace(0.5, N2 - 0.5, N), np.linspace(0.5, M2 - 0.5, M))
psi = np.zeros((M + 1, N + 1))
a_psi = np.zeros((2 * M2 + 1, 2 * N2 + 1))
m = np.zeros((M, N))
with open(f'results3/psi{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        psi[i, :] = list(map(float, line.split()))

with open(f'results3/m{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        m[i, :] = list(map(float, line.split()))

with open(f'result2.dat', 'r') as file:
    file.readline()
    file.readline()
    file.readline()
    for j in range(2 * N2 + 1):
        for i in range(2 * M2 + 1):
            a_psi[i][j] = list(map(float, file.readline().split()))[2]

print(f'Пределы a_psi: {np.min(a_psi)} – {np.max(a_psi)}')
min_psi = np.min(psi)
max_psi = np.max(psi)
min_v_psi = np.min(a_psi)
max_v_psi = np.max(a_psi)

a_psi -= min_v_psi
a_psi *= (max_psi - min_psi) / (max_v_psi - min_v_psi)
a_psi += min_psi

plt.figure(1)
plt.axes().set_aspect('equal')
plt.contour(X, Y, psi, np.linspace(np.min(psi), np.max(psi), 30))
plt.colorbar()
plt.title('Функция тока')
plt.show()
#
# plt.figure(2)
# plt.axes().set_aspect('equal')
# plt.contour(X1, Y1, a_psi, np.linspace(np.min(a_psi), np.max(a_psi), 20))
# plt.colorbar()
# plt.title('Функция тока (детальная)')
# plt.show()

plt.figure(3)
plt.axes().set_aspect('equal')
plt.contourf(X2, Y2, m, np.linspace(np.min(m), np.max(m), 100))
CS1 = plt.contour(X, Y, psi, np.linspace(np.min(psi), np.max(psi), 10), colors='red')
plt.clabel(CS1, fontsize=7)
plt.colorbar()
CS2 = plt.contour(X1, Y1, a_psi, np.linspace(np.min(a_psi), np.max(a_psi), 10), colors='blue')
plt.clabel(CS2, fontsize=7)
plt.colorbar()
plt.title('2 функции тока')
plt.show()
#
# plt.figure(4)
# plt.axes().set_aspect('equal')
# plt.contourf(X2, Y2, m, np.linspace(np.min(m), np.max(m), 6))
# plt.colorbar()
# plt.title('пористость')
# plt.show()
#
print(f'Пределы psi: {np.min(psi)} – {np.max(psi)}')
print(f'Пределы a_psi: {np.min(a_psi)} – {np.max(a_psi)}')
