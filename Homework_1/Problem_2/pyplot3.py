import numpy as np
import matplotlib.pyplot as plt

# Load the Fortran data:
data = np.loadtxt('./fort.17')

# Separate x and y values
x = data[:, 0]
y0 = data[:, 1]
plt.plot(x, y0, linestyle='-', label = "n=0")
plt.legend()
plt.xlabel('x')
plt.ylabel(r'$\psi(x)$')
plt.title(r'Wavefunctions for 1D Harmonic Oscillator')
plt.grid(True)
plt.show()
