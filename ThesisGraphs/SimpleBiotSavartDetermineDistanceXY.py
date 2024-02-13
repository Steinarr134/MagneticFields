import numpy as np
from matplotlib import pyplot as plt

"""
This program was written by Steinarr Hrafn to generate a plot for his Master Thesis

The output is a plot of the strength of the x and y components of the B field
as a function of sag. 

In this simulation the approximation used are that the currents in the three phases
are perfectly balanced and the wires are infinitely long and straight.

sag in this context means lowering the wires from the starting position of Y0

"""

f = 50  # [Hz] frequency of transmission
x_p = 9  # [m] Phase spacing
Y0 = 20  # [m] Starting position of the wires
r = np.array([0, 0])  # [m] sensing location, in (x, y)
mu = 1.25663753e-6  # [?] magnetic permeability of air
I_rms = 400  # [A] RMS current in wires
I_0 = np.sqrt(2)*I_rms  # [A] amplitude of current
N = 401  # number of points in time
cycles = 2  # length of simulation in cycles of the frequency f
t_space = np.linspace(0, cycles/f, N)  # [s] timespace

N_sag = 100  # asldfkj
sagspace = np.linspace(0.1, 12, N_sag)  # [m] all the values of sag that are to be considered

# h will be calculated by Y0-sag
hspace = Y0 - sagspace

# Arrays to hold the results
Bx = np.zeros(N_sag)
By = np.zeros(N_sag)


def rms(signal):
    return np.sqrt(np.mean(np.square(signal)))


for i, h in enumerate(hspace):
    phase_locations = np.array([[-x_p, h], [0, h], [x_p, h]])  # (x, y), in meters
    # phase_locations = np.array([[0, 9], [0, 14], [0, 19]])  # (x, y), in meters  -  vertical arrangement

    I = np.zeros((N, 3))  # current in each phase as a function of time
    _B = np.zeros((N, 2))  # magnetic flux density at r in (x, y) components as a function of time

    k = np.array([0, 0, 1])  # (x, y, z), basis vector for z-direction
    for p in range(3):  # loop through phases, p is the phase number
        # calculate r' in (x, y, z) coordinates
        r_prime = np.hstack((r - phase_locations[p, :], 0))  # in (x, y, z)
        lenght_r_prime = np.sqrt(np.sum(np.square(r_prime)))

        # calculate the current through this phase over the timespace
        I_p = I_0 * np.sin(2 * np.pi * f * t_space + (p - 1) * (2 * np.pi / 3))

        # Save the current of this phase to I
        I[:, p] = I_p

        # change the current array into a column vector
        I_p = np.atleast_2d(I_p).T

        # calculate the affect this current has on the B field (equation ?? in thesis)  TODO fix equation number
        # The current is flowing parallell to the z-axis thus the cross product
        # is taken between k and r'
        B_p = (mu/(2*np.pi)) * (I_p / lenght_r_prime ** 2) * np.cross(k, r_prime)

        # add the x and y components to the B array
        _B += B_p[:, :2]
    # Add the rms values of the x and y components to Bx and By
    Bx[i] = rms(_B[:, 0])
    By[i] = rms(_B[:, 1])

# The Bx and By arrays now hold the x and y components of the B field over sagspace
# The rest is just about plotting the results

# plt.rcParams['text.usetex'] = True
# set the figure size to be pretty
goldenratio = (5**.5-1)/2
width = 6
figsize = (width, width*goldenratio)
f, ax = plt.subplots(1, 1, figsize=figsize)
f.suptitle("Simulation of 3 horizontal phases \nusing Biot-Savart infinite line simplification")

ax.plot(sagspace, Bx*1e6, 'tab:orange')
ax.plot(sagspace, By*1e6, 'm')
ax.set_xlabel("sag [m]")
ax.set_ylabel(r" $\mathbf{B}$ [ $\mu T$ ]")
ax.legend(("x", "y"), loc='upper left')

# f.savefig("Graphs\\SimpleBiotSavartThreePhaseHorizontal.pdf", format='pdf', bbox_inches='tight')
plt.show()