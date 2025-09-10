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
plt.plot(x, y1, marker='o', linestyle='-')
plt.plot(x, y2, marker='o', linestyle='-')
plt.plot(x, y3, marker='o', linestyle='-')
plt.plot(x, y4, marker='o', linestyle='-')
plt.plot(x, y5, marker='o', linestyle='-')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Fortran Data Plot')
plt.grid(True)
plt.show()
