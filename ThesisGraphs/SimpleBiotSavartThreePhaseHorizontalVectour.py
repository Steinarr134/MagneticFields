import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker

"""
This program was written by Steinarr Hrafn to generate a plot for his Master Thesis

The output is a vector and contour plot (vectour) showing the magnetic field (B field) around
a three phases arranged horizontally
"""


f = 50  # frequency in Hz
phase_locations = np.array([[-9, 12], [0, 12], [9, 12]])  # (x, y), in meters
mu = 1.25663753e-6  # magnetic permeability of air
I_rms = 400  # RMS current in wires, in amps
I_0 = np.sqrt(2)*I_rms

N = 100  # number of points being simulated

space = ((-15, -2), (15, 18))
x_space = np.linspace(space[0][0], space[1][0], N)
y_space = np.linspace(space[0][1], space[1][1], N)
X, Y = np.meshgrid(x_space, y_space)


width = 6
figsize = (6, 10)

timepoints = np.array([0, 5])/1000
for n, t in enumerate(timepoints):
    fig, axs = plt.subplots(2,1, figsize=figsize, sharex=True)
    axs = np.ravel(axs)
    Ip = [I_0 * np.sin(2*np.pi*f*t + (p-1)*(2*np.pi/3)) for p in range(3)]  # current in each phase

    B = np.zeros((N, N, 2))  # magnetic flux density
    B_strength = np.zeros((N, N))
    B_direction = np.zeros((N, N, 2))
    B_phase = np.zeros((N, N, 2, 3))
    B_phase_strength = np.zeros((N, N, 3))
    B_phase_direction = np.zeros((N, N, 2, 3))
    k = np.array([0, 0, 1])  # (x, y, z)
    for j, x in enumerate(x_space):
        for i, y in enumerate(y_space):
            for p in range(3):
                r_prime = np.array([x, y]) - phase_locations[p, :]
                lenght_r_prime = np.sqrt(np.sum(np.square(r_prime)))
                _B = (mu/(2*np.pi))*(Ip[p]/lenght_r_prime**2)*np.cross(k, r_prime)
                B[i, j, :] += _B[:2]
                B_phase[i, j, :, p] = _B[:2]
                B_phase_strength[i, j, p] = np.sqrt(np.sum(np.square(_B[:2])))
                if B_phase_strength[i, j, p] == 0:
                    B_phase_strength[i, j, p] = 1e-12
                B_phase_direction[i, j, :, p] = B_phase[i, j, :, p]/B_phase_strength[i, j, p]
            B_strength[i, j] = np.sqrt(np.sum(np.square(B[i, j, :])))
            B_direction[i, j, :] = B[i, j, :]/B_strength[i, j]

    a = 6
    colors = "brg"
    for p in range(3):
        axs[0].contour(X, Y, B_phase_strength[:, :, p], locator=ticker.LogLocator(subs=(1, 2, 3, 4, 5, 6, 7, 8, 9)),
                       colors=colors[p], linewidths=1)
        axs[0].quiver(X[::a, ::a], Y[::a, ::a], B_phase_direction[::a, ::a, 0, p], B_phase_direction[::a, ::a, 1, p], scale=35, color=colors[p])
    axs[1].contour(X, Y, B_strength, locator=ticker.LogLocator(subs=(1, 2, 3, 4, 5, 6, 7, 8, 9)))
    axs[1].quiver(X[::a, ::a], Y[::a, ::a], B_direction[::a, ::a, 0], B_direction[::a, ::a, 1], scale=35, )
    fig.suptitle(r"Vector and contour plot of $\mathbf{B}$ in"
                 +"\nvicinity of 3 phase OTL, at " + f"$t={int(t*1000)}$ ms")
    axs[1].set_xlabel("x location[m]")
    axs[0].set_ylabel("y location[m]")
    axs[1].set_ylabel("y location[m]")
    fig.savefig(f"Graphs\\SimpleBiotSavartThreePhaseHorizontalVectour{int(t*1000)}.pdf", format='pdf', bbox_inches='tight')
plt.show()