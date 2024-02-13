import numpy as np
from matplotlib import pyplot as plt

"""
This program was written by Steinarr Hrafn to generate a plot for his Master Thesis

The output is a plot of the current in three phases of an otl and the x and y components
of the magnetic field (B field) at a sensing location 12 meters below the otl wires.
"""

f = 50  # [Hz] frequency of transmission
x_p = 9  # [m] Phase spacing
h = 12  # [m] Height of wires
phase_locations = np.array([[-x_p, h], [0, h], [x_p, h]])  # (x, y), in meters
# phase_locations = np.array([[0, 9], [0, 14], [0, 19]])  # (x, y), in meters  -  vertical arrangement
r = np.array([0, 0])  # sensing location, (x, y), in meters
mu = 1.25663753e-6  # magnetic permeability of air
I_rms = 400  # [A] RMS current in wires
I_0 = np.sqrt(2)*I_rms  # [A] amplitude of current

N = 401  # number of points in time
cycles = 2  # length of simulation in cycles of the frequency f
t_space = np.linspace(0, cycles/f, N)  # timespace

I = np.zeros((N, 3))  # current in each phase as a function of time
B = np.zeros((N, 2))  # magnetic flux density at r in (x, y) components as a function of time

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
    B += B_p[:, :2]

# The B array now holds the x and y components of the B field over tspace
# And I array holds the currents in the three phases over tspace
# The rest is just about plotting the results

# plt.rcParams['text.usetex'] = True
# set the figure size to be pretty
goldenratio = (5**.5-1)/2
width = 6
figsize = (width, width*goldenratio)
f, (ax1, ax2) = plt.subplots(2, 1, sharex='col', figsize=figsize)
f.suptitle("Simulation of 3 horizontal phases \nusing Biot-Savart infinite line simplification")
ax1.plot(t_space*1000, I)
ax1.set_ylabel(r"Current, $I$ [ $A$ ]")
ax1.legend(("left", "center", "right"))
# ax2.plot(t_space*1000, B*1e6)
ax2.plot(t_space*1000, B[:, 0]*1e6, 'tab:orange')
ax2.plot(t_space*1000, B[:, 1]*1e6, 'm')
ax2.set_xlabel("Time [ms]")
ax2.set_ylabel(r" $\mathbf{B}$ [ $\mu T$ ]")
ax2.legend(("x", "y"), loc='upper right')

print(t_space.shape, I.shape)
print(np.hstack((np.transpose([t_space]), I))[:55,:])
# f.savefig("Graphs\\SimpleBiotSavartThreePhaseHorizontal.pdf", format='pdf', bbox_inches='tight')
plt.show()