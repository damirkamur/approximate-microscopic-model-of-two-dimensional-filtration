from PIL import Image


def interpolate(color):
    return minM + color * (maxM - minM) / 255


image = Image.open('field.png').convert('L')
N, M = 400, 200
maxM = 0.2
minM = 0.05

image = image.resize((N, M))
data = image.load()

with open('m.txt', 'w') as file:
    for j in range(M - 1, -1, -1):
        for i in range(N):
            file.write('%10.5f' % interpolate(data[i, j]))
        file.write('\n')
