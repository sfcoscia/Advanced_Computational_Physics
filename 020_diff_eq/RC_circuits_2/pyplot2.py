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
plt.xlabel('t')
plt.ylabel(r'Q')
plt.title(r'Q vs t')
plt.savefig('Q-v-t.png')
plt.grid(True)
#plt.show()

plt.figure()
plt.plot(t, current, linestyle='-', label = "n=0")
plt.legend()
plt.xlabel('t')
plt.ylabel(r'I')
plt.title(r'I vs t')
plt.savefig('I-v-t.png')
plt.grid(True)
#plt.show()

