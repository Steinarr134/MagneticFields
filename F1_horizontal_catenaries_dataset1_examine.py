import numpy as np
import pickle
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

with open("dataset1.pkl", 'rb') as f:
    Is, sags, Y0, D, Data = pickle.load(f)

Data = Data*1e6

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X, Y = np.meshgrid(sags, Is)

# Plot the surface.
surf = ax.plot_surface(X, Y, Data, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(np.min(Data), np.max(Data))
# ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.xlabel("Sag [m]")
plt.ylabel("Peak Current [A]")
ax.set_zlabel("Peak magnetic strength [uT]")


print(Data)
plt.show()