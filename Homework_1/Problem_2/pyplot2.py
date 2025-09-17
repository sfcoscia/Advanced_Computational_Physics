import numpy as np
import matplotlib.pyplot as plt

# Function to plot the figures
def plotfig(data,l,D):

    plt.figure()
    x = data[:, 0]
    y0 = data[:, 1]
    y1 = data[:, 2]
    y2 = data[:, 3]
    y3 = data[:, 4]
    y4 = data[:, 5]
    y5 = data[:, 6]

    plt.plot(x, y0, linestyle='-', label = "n=0")
    plt.plot(x, y1, linestyle='-', label = "n=1")
    plt.plot(x, y2, linestyle='-', label = "n=2")
    plt.plot(x, y3, linestyle='-', label = "n=3")
    plt.plot(x, y4, linestyle='-', label = "n=4")
    plt.plot(x, y5, linestyle='-', label = "n=5")
    plt.legend()
    plt.xlabel('r')
    plt.ylabel(r'$\psi_{n\ell}(v,r)$')
    plt.title(f'Radial Wavefunctions for the {D}D Harmonic Oscillator with ' + r'$\ell$ = ' + f'{l}')
    plt.grid(True)
    plt.savefig(f'./plots/Wf_{l}l_{D}-D.png')

# Load the Fortran data:
data = np.loadtxt('./fort.2')
print(data)
plotfig(data,0,2)
data = np.loadtxt('./fort.3')
print(data)
plotfig(data,0,3)
data = np.loadtxt('./fort.12')
print(data)
plotfig(data,1,2)
data = np.loadtxt('./fort.13')
print(data)
plotfig(data,1,3)
