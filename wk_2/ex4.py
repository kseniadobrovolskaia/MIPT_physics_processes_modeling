import numpy as np
import  matplotlib.pyplot as plt

A = 12
B = -13

def Func(x):
    return A * x + B

mu, sigma = 0, 25
Gaus = np.random.normal(mu, sigma, 100)
One = np.random.uniform(-25, 25, 100)
Stand = np.random.standard_normal(100)

X = np.linspace(1, 10, 100)
Y = [Func(x) for x in X]
Yone = np.zeros(100)
YoneL = np.zeros(100)
Ygaus = np.zeros(100)
YgausL = np.zeros(100)
Ystand = np.zeros(100)
YstandL = np.zeros(100)

fig, axs = plt.subplots(1, 3)
fig.suptitle('Least squares fit')

for i in range(100):
    F = Func(X[i])
    Yone[i] = F + One[i]
    Ygaus[i] = F + Gaus[i]
    Ystand[i] = F + Stand[i]

A, B = np.polyfit(X, Yone, 1)
YoneL = [A * x + B for x in X]

A, B = np.polyfit(X, Ygaus, 1)
YgausL = [A * x + B for x in X]

A, B = np.polyfit(X, Ystand, 1)
YstandL = [A * x + B for x in X]

axs.flat[0].set(xlabel='Равномерное')
axs[0].scatter(X, Yone)
axs[0].plot(X, Y, "g")
axs[0].plot(X, YoneL, "r")

axs.flat[1].set(xlabel='Гауссовское')
axs[1].scatter(X, Ygaus)
axs[1].plot(X, Y, "g")
axs[1].plot(X, YgausL, "r")

axs.flat[2].set(xlabel='Стандартное')
axs[2].scatter(X, Ystand)
axs[2].plot(X, Y, "g")
axs[2].plot(X, YstandL, "r")

for ax in axs.flat:
    ax.label_outer()

axs.flat[0].grid()
axs.flat[1].grid()
axs.flat[2].grid()
plt.show()
