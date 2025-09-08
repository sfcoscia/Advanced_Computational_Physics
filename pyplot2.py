import numpy as np
import matplotlib.pyplot as plt

# Read data from the text file
# For formatted data, np.loadtxt is convenient
data = np.loadtxt('./200/fort.19')

# Separate x and y values
x = data[:, 0]
y = data[:, 2]
print(x, y)
# Create the plot
plt.plot(x, y, marker='o', linestyle='-')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Fortran Data Plot')
plt.grid(True)
plt.show()
