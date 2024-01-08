import numpy as np
from matplotlib import pyplot as plt

"""
This program was written by Steinarr Hrafn to generate a plot for his Master Thesis

The output is a plot of the current in three phases of an otl and the x and y components
of the magnetic field (B field) at a sensing location 12 meters below the otl wires.
"""

f = 50  # frequency in Hz
phase_locations = np.array([[-9, 12], [0, 12], [9, 12]])  # (x, y), in meters
# r = np.array([0, 0])  # sensing location, (x, y), in meters
mu = 1.25663753e-6  # magnetic permeability of air
I_rms = 400  # RMS current in wires, in amps
I_0 = np.sqrt(2)*I_rms

N_t = 201  # number of points in time
cycles = 1  # length of simulation in cycles of the frequency f
t_space = np.linspace(0, cycles/f, N_t)  # timespace

N = 201
hspace = np.linspace(0, 8, N)
B_rms_vs_h = np.zeros((N, 2))  # magnetic flux density at r in (x, y) components

k = np.array([0, 0, 1])  # (x, y, z)
for i, h in enumerate(hspace):  # loop through the height, (y coordinate)
    r = np.array([0, h])

    B = np.zeros((N_t, 2))
    for p in range(3):  # loop through phases
        r_prime = np.hstack((r - phase_locations[p, :], 0))  # in (x, y, z)
        lenght_r_prime = np.sqrt(np.sum(np.square(r_prime)))

        for j, t in enumerate(t_space):  # loop through timespace
            # calculate the current through this phase at this time
            I = I_0 * np.sin(2*np.pi*f*t + (p-1)*(2*np.pi/3))
            # calculate the affect this current has on the B field (equation ?? in thesis)  TODO fix equation number
            _B = (mu/(2*np.pi))*(I/lenght_r_prime**2)*np.cross(k, r_prime)

            # add the x and y components to the B array
            B[j, :] += _B[:2]
    B_rms_vs_h[i, :] = np.sqrt(np.mean(np.square(B), axis=0))
# plt.rcParams['text.usetex'] = True
goldenratio = (5**.5-1)/2
width = 6
figsize = (width, width*goldenratio)
# f, (ax1, ax2) = plt.subplots(2, 1, sharex='col', figsize=figsize)
plt.figure(figsize=figsize)
plt.suptitle(f"Simulation of 3 horizontal phases, "
             r"$I_{rms} = " + f"{I_rms}" + r" A$" + "\n"
             f"using Biot-Savart infinite line simplification")
plt.plot(hspace, B_rms_vs_h*1e6)
plt.xlabel("Height on y-axis [m]")
plt.ylabel(r" $\mathbf{B}$ [ $\mu T$ ]")
plt.legend(("x", "y"), loc='upper left')

# also plot analytical solution
# plt.plot(hspace, 80*(1/(12-hspace) - (12-hspace)/((12-hspace)**2 + 81)))
# plt.plot(hspace, 720*np.sqrt(3)/((12 - hspace)**2+81))

# print(t_space.shape, I.shape)
# print(np.hstack((np.transpose([t_space]), I))[:55,:])
# plt.savefig("Graphs\\Brms_over_yaxis_simpleBiotSavart.pdf", format='pdf', bbox_inches='tight')
plt.show()