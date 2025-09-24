import numpy as np
import matplotlib.pyplot as plt

# Load the Fortran data:
data = np.loadtxt('./fort.1')
data2 = np.loadtxt('./fort.2')

# Separate x and y values
x = data[:, 0]
y0 = data[:, 1]

t = data2[:,0]
current = data2[:, 1]
#print(x, y)
# Create the plot
plt.figure()
plt.plot(x, y0, linestyle='-', label = "n=0")
plt.legend()
plt.xlabel('x')
plt.ylabel(r'y')
plt.title(r'y vs x')
plt.savefig('y-v-x.png')
plt.grid(True)
#plt.show()

plt.figure()
plt.plot(t, current, linestyle='-', label = "n=0")
plt.legend()
plt.xlabel('t')
plt.ylabel(r'W')
plt.title(r'W vs t')
plt.savefig('W-v-t.png')
plt.grid(True)
#plt.show()

