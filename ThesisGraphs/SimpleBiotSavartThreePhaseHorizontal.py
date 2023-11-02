import numpy as np
from matplotlib import pyplot as plt

f = 50  # frequency in Hz
phase_locations = np.array([[-9, 12], [0, 12], [9, 12]])  # (x, y), in meters
r = np.array([0, 0])  # sensing location, (x, y), in meters
mu = 1.25663753e-6  # magnetic permeability of air
I_rms = 400  # RMS current in wires, in amps
I_0 = np.sqrt(2)*I_rms

N = 250  # number of points in time
cycles = 2  # length of simulation in cycles of the frequency f
t_space = np.linspace(0, cycles/f, N)  # timespace

I = np.zeros((N, 3))  # current in each phase
B = np.zeros((N, 2))  # magnetic flux density at r in (x, y) components

k = np.array([0, 0, 1])  # (x, y, z)
for p in range(3):  # loop through phases
    r_prime = np.hstack((r - phase_locations[p, :], 0))  # in (x, y, z)
    lenght_r_prime = np.sqrt(np.sum(np.square(r_prime)))
    for i, t in enumerate(t_space):  # loop through timespace
        # calculate the current through this phase at this time
        I[i, p] = I_0 * np.sin(2*np.pi*f*t + (p-1)*(2*np.pi/3))
        # calculate the affect this current has on the B field (equation ?? in thesis)  TODO fix equation number
        _B = (mu/(2*np.pi))*(I[i, p]/lenght_r_prime**2)*np.cross(k, r_prime)

        # add the x and y components to the B array
        B[i, :] += _B[:2]


# plt.rcParams['text.usetex'] = True
f, (ax1, ax2) = plt.subplots(2, 1, sharex='col')
f.suptitle("Simulation of 3 horizontal phases \nusing Biot-Savart infinite line simplification")
ax1.plot(t_space*1000, I)
ax1.set_ylabel(r"Current, $I$ [ $A$ ]")
ax1.legend(("left", "center", "right"))
ax2.plot(t_space*1000, B*1e6)
ax2.set_xlabel("Time [ms]")
ax2.set_ylabel(r" $\mathbf{B}$ [ $\mu T$ ]")
ax2.legend(("x-component", "y-component"), loc='upper right')

plt.savefig("Graphs\\SimpleBiotSavartThreePhaseHorizontal.png")
plt.show()