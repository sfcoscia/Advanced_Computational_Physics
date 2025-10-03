import numpy as np
import matplotlib.pyplot as plt

# Load the Fortran data:
data = np.loadtxt('./fort.1')
data2 = np.loadtxt('./fort.2')

# Separate x and y values
x = data[:, 0]
y0 = data[:, 2]
y1 = data[:,3]

x1 = data2[:, 0]
y2 = data2[:, 2]
y3 = data2[:,3]

#t = data2[:,0]
#current = data2[:, 1]
#print(x, y)
# Create the plot
plt.figure()
plt.plot(x, y0, linestyle='-', label = "n=0")
plt.legend()
plt.xlabel('x')
plt.ylabel(r'y')
plt.title(r'y vs x')
plt.grid(True)
plt.plot(x, y1, linestyle='-', label = "n=0")
plt.plot(x1, y2, linestyle='-', label = "n=0")
plt.plot(x1, y3, linestyle='-', label = "n=0")
plt.legend()
plt.grid(True)
plt.savefig('y-v-x.png')
plt.show()
