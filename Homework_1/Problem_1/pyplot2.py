import numpy as np
import matplotlib.pyplot as plt

# Load the Fortran data:
data = np.loadtxt('./fort.7')

# Separate x and y values
x = data[:, 0]
y1 = data[:, 1]
y2 = data[:, 2]
y3 = data[:, 3]
y4 = data[:, 4]
y5 = data[:, 5]
#print(x, y)
# Create the plot
plt.plot(x, y1, linestyle='-', label = "n=1")
plt.plot(x, y2, linestyle='-', label = "n=2")
plt.plot(x, y3, linestyle='-', label = "n=3")
plt.plot(x, y4, linestyle='-', label = "n=4")
plt.plot(x, y5, linestyle='-', label = "n=5")
plt.legend()
plt.xlabel('x')
plt.ylabel(r'$\psi(x)$')
plt.title(r'Wavefunctions for 1D Harmonic Oscillator')
plt.grid(True)
plt.show()
