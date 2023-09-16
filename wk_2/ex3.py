import numpy as np
from collections import namedtuple

N = 10
MatrixArray = np.zeros(shape = (N, N))

for i in range(N):
    MatrixArray[i][i] = 2
    if (i != 0):
        MatrixArray[i][i - 1] = -1
    if (i != (N - 1)):
        MatrixArray[i][i + 1] = -1

Rights = np.zeros(N)
Rights[0] = 100
Rights[N - 1] = 1

x = np.linalg.solve(MatrixArray, Rights)

print(x)
