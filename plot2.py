from matplotlib import pyplot as plt
import numpy as np

N = 480
M = N * 2 // 3
X, Y = np.meshgrid(np.linspace(0, N + 1, N + 1), np.linspace(0, M + 1, M + 1))
X1, Y1 = np.meshgrid(np.linspace(0.5, N - 0.5, N), np.linspace(0.5, M - 0.5, M))
psi = np.zeros((M + 1, N + 1))
a_psi = np.zeros((M + 1, N + 1))
m = np.zeros((M, N))
with open(f'results2/psi{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        psi[i, :] = list(map(float, line.split()))
with open(f'results2/analit_psi{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        a_psi[i, :] = list(map(float, line.split()))

with open(f'results2/m{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        m[i, :] = list(map(float, line.split()))

# plt.figure(1)
# plt.axes().set_aspect('equal')
# plt.contourf(X, Y, psi, np.linspace(np.min(psi), np.max(psi), 30))
# plt.colorbar()
# plt.title('Функция тока')
# plt.show()
#
# plt.figure(2)
# plt.axes().set_aspect('equal')
# plt.contour(X, Y, a_psi, np.linspace(np.min(a_psi), np.max(a_psi), 20))
# plt.colorbar()
# plt.title('Функция тока (аналитическая)')
# plt.show()

plt.figure(3)
plt.axes().set_aspect('equal')
plt.contourf(X1, Y1, m, np.linspace(np.min(m), np.max(m), 100))
CS1 = plt.contour(X, Y, psi, np.linspace(np.min(psi), np.max(psi), 20), colors='red')
plt.clabel(CS1, fontsize=7)
plt.colorbar()
CS2 = plt.contour(X, Y, a_psi, np.linspace(np.min(a_psi), np.max(a_psi), 20), colors='blue')
plt.clabel(CS2, fontsize=7)
plt.colorbar()
plt.title('2 функции тока')
plt.show()

# plt.figure(4)
# plt.axes().set_aspect('equal')
# plt.contourf(X1, Y1, m, np.linspace(np.min(m), np.max(m), 100))
# plt.colorbar()
# plt.title('пористость')
# plt.show()

print(f'Пределы psi: {np.min(psi)} – {np.max(psi)}')
print(f'Пределы a_psi: {np.min(a_psi)} – {np.max(a_psi)}')
