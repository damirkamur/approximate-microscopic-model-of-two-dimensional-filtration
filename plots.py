from matplotlib import pyplot as plt
import numpy as np

N, M = 3, 4
psi = np.zeros((M + 1, N + 1))
p = np.zeros((M + 1, N + 1))
with open(f'results/psi{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        psi[i, :] = list(map(float, line.split()))
with open(f'results/p{N}x{M}.dat', 'r') as file:
    for i, line in enumerate(file):
        p[i, :] = list(map(float, line.split()))

plt.figure(1)
plt.contourf(psi)
plt.colorbar()
plt.title('Функция тока')
plt.show()

plt.figure(2)
plt.contourf(p)
plt.colorbar()
plt.title('Функция давления')
plt.show()
