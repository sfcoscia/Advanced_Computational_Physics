import numpy as np
import matplotlib.pyplot as plt

# Load the Fortran data:
data = np.loadtxt('./fort.1')
data2 = np.loadtxt('./fort.3')
data3 = np.loadtxt('./fort.4')


# Separate x and y values
x0 = data[:, 0]
y0 = data[:, 1]

x1 = data2[:,0]
y1 = data2[:, 1]

x2 = data3[:,0]
y2 = data3[:,1]
#print(x, y)
# Create the plot
plt.figure()
plt.plot(x0, y0, linestyle='-', label = "n=0")
plt.legend()
plt.xlabel('x')
plt.ylabel(r'y')
plt.title(r'y vs x')
#plt.show()

plt.plot(x1, y1, linestyle='-', label = "n=0")
plt.plot(x2, y2, linestyle='-', label = "n=0")
plt.savefig('y-v-x.png')
plt.show()

