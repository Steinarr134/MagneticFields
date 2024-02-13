import numpy as np
import pickle
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

with open("D250_Yn20_z.pkl", 'rb') as f:
    Is, sags, Y0, D, Bs1, Bs2 = pickle.load(f)

Bs1 = Bs1*1e6
Bs2 = Bs2*1e6

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X, Y = np.meshgrid(sags, Is)

# Plot the surface.
surf = ax.plot_surface(X, Y, Bs1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# surf = ax.plot_surface(X, Y, Bs2, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# surf = ax.plot_surface(X, Y, Bs1-Bs2, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(np.min(Bs1), np.max(Bs1))
# ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
# ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.xlabel("Sag [m]")
plt.ylabel("Peak Current [A]")
ax.set_zlabel("Peak magnetic strength [uT]")
# plt.show()

"""contour plot of the two surfaces"""
plt.figure()
levels = np.arange(np.min(Bs2), np.max(Bs1), 0.5)
plt.contour(X, Y, Bs1, levels=levels)
# plt.contour(X, Y, Bs2, levels=levels)
plt.contour(X, Y, Bs1-Bs2, levels=np.arange(0, 1, 0.1), colors="red")
plt.show()

# it might work, the contour lines of Bs1-Bs2 and Bs1 are not perpendicular and each pair of lines appears
# to intersect only once so one might be able to induce the location in (sag, I) coordinates based on the
# intersection point of B1 and B1-B2 contour lines
