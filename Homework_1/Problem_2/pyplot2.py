import numpy as np
import matplotlib.pyplot as plt

# Load the Fortran data:
data = np.loadtxt('./fort.18')

# Separate x and y values
x = data[:, 0]
y1 = data[:, 1]
#print(x, y)
# Create the plot
plt.plot(x, y1, linestyle='-', label = "n=1")
plt.legend()
plt.xlabel('x')
plt.ylabel(r'$\psi(x)$')
plt.title(r'Wavefunctions for 1D Harmonic Oscillator')
plt.grid(True)
plt.show()
