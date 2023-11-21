import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker

"""
This program was written by Steinarr Hrafn to generate a plot for his Master Thesis

The output is a vector and contour plot (vectour) showing the magnetic field (B field) around
a single current carrying wire
"""


f = 50  # frequency in Hz
wire_location = np.array([0, 0])
mu = 1.25663753e-6  # magnetic permeability of air
I = 100  # current in wire, in amps

N = 100  # number of points being simulated
size = (1, 1)  # the size of the simulation in meters

x_space = np.linspace(wire_location[0]-size[0]/2, wire_location[0]+size[0]/2, N)
y_space = np.linspace(wire_location[1]-size[1]/2, wire_location[1]+size[1]/2, N)
X, Y = np.meshgrid(x_space, y_space)

B = np.zeros((N, N, 2))  # magnetic flux density
B_strenght = np.zeros((N, N))
B_direction = np.zeros((N, N, 2))
k = np.array([0, 0, 1])  # (x, y, z)
for i, x in enumerate(x_space):
    for j, y in enumerate(y_space):
        r_prime = np.array([x, y]) - wire_location
        lenght_r_prime = np.sqrt(np.sum(np.square(r_prime)))
        _B = (mu/(2*np.pi))*(I/lenght_r_prime**2)*np.cross(k, r_prime)
        B[i, j, :] = _B[:2]
        B_strenght[i, j] = np.sqrt(np.sum(np.square(B[i, j, :])))
        B_direction[i, j, :] = np.cross(k, r_prime/lenght_r_prime)[:2]
#
# for p in range(3):  # loop through phases
#     r_prime = np.hstack((r - phase_locations[p, :], 0))  # in (x, y, z)
#     lenght_r_prime = np.sqrt(np.sum(np.square(r_prime)))
#     for i, t in enumerate(t_space):  # loop through timespace
#         # calculate the current through this phase at this time
#         I[i, p] = I_0 * np.sin(2*np.pi*f*t + (p-1)*(2*np.pi/3))
#         # calculate the affect this current has on the B field (equation ?? in thesis)  TODO fix equation number
#         _B = (mu/(2*np.pi))*(I[i, p]/lenght_r_prime**2)*np.cross(k, r_prime)
#
#         # add the x and y components to the B array
#         B[i, :] += _B[:2]


# plt.rcParams['text.usetex'] = True
goldenratio = (5**.5-1)/2
width = 6
figsize = (width, width)
f = plt.figure(figsize=figsize)
plt.contour(X, Y, B_strenght, locator=ticker.LogLocator(subs=(1, 2, 3, 4, 5, 6, 7, 8, 9)))
a = 4
# logBdirection = B_direction*np.dstack((B_strenght, B_strenght))
print(B_direction[::a, ::a, 0])
plt.quiver(X[::a, ::a], Y[::a, ::a], B_direction[::a, ::a, 1], B_direction[::a, ::a, 0])
f.suptitle("Magnetic field lines in vicinity of current, $I=100 A$ \n flowing into the page")
plt.xlabel("x location [m]")
plt.ylabel("y location [m]")
f.savefig("Graphs\\SimpleBiotSavartSingleWireVectour.pdf", format='pdf', bbox_inches='tight')
plt.show()