import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

w, e = np.loadtxt('/Users/edoardo/Desktop/simulazione_prova/fit.txt', unpack=True)


def F(w, a, b):
    return a*w/np.tanh(b*w)


init = (0.01, 0.01)

popt, pcov = curve_fit(F, e, w, init)

print(popt)

x = np.linspace(0.01, 2, 10000)
y = F(x, 0.5, 2)

fig, ax = plt.subplots()

ax.scatter(w, e, c='blue', s=2, label='computed')
ax.scatter(x, y, c='green', s=0.01, label='theoretical')
ax.set_title('Energy of the quantum harmonic oscillator')
ax.set_ylabel('Energy [a.u.]')
ax.set_xlabel('$\omega_{ext}$ [a.u]')
ax.legend()
plt.show()

