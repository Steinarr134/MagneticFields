import numpy  as np
from matplotlib import pyplot as plt
from matplotlib import cm
plt.style.use("dark_background")

e0 = 0.05
e1 = 0.02
h = 12
d_P = 8
d = 0.5


def func(x):
    # return 2.09385992/np.sinh(x - 1.0043624) +30.68953421
    # return 4.66487133e+03/np.cosh(x -8.58502502e-01) -4.65016916e+03
    return 4.66487133e+03/np.cosh(x -8.58502502e-01) -4.65016916e+03


th0 = np.linspace(-np.pi, np.pi, 250)
th1 = np.linspace(-np.pi, np.pi, 250)


fig, (ax11, ax12) = plt.subplots(ncols=2, subplot_kw={"projection": "3d"})

# Make data.
TH0, TH1 = np.meshgrid(th0, th1)
Cy = 2*e0*np.cos(TH0)/np.sqrt(3) - e1*np.sin(TH1)
Dy = 1 + e1*np.cos(TH1) + 2*e0*np.sin(TH0)/np.sqrt(3)
Zy = np.sqrt(np.square(Cy) + np.square(Dy))

# Plot the surface.
surf = ax11.plot_surface(TH0, TH1, Zy, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

ax11.view_init(elev=12, azim=-150, roll=0)
a = (1 + 3*h**2/(d_P**2))
Cx = 1+ e1*np.cos(TH1) - a*e0*np.cos(TH0)
Dx = e1*np.sin(TH1) + a*e0*np.sin(TH0)
Zx = np.sqrt(np.square(Cx) + np.square(Dx))

# Plot the surface.
surf = ax12.plot_surface(TH0, TH1, Zx, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)

ax12.view_init(elev=12, azim=-150, roll=0)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)


fig2, (ax21, ax22) = plt.subplots(ncols=2, subplot_kw={"projection": "3d"})
a = (1 + 3*(h+d)**2/(d_P**2))
Cx2 = 1+ e1*np.cos(TH1) - a*e0*np.cos(TH0)
Dx2 = e1*np.sin(TH1) + a*e0*np.sin(TH0)
Zx2 = np.sqrt(np.square(Cx2) + np.square(Dx2))

surf = ax21.plot_surface(TH0, TH1, Zy/Zx, cmap=cm.coolwarm,
                         linewidth=0, antialiased=False)
surf = ax22.plot_surface(TH0, TH1, 10*Zx2/Zx-8.7, cmap=cm.coolwarm,
                         linewidth=0, antialiased=False)


fig3, ax3 = plt.subplots(subplot_kw={"projection": "3d"})

r1 = h
r2 = h+d
A = (Zx2/(r2*(r2**2 + d_P**2)))/(Zx/(r1*(r1**2 + d_P**2)))

def mid(mat):
    return np.min(mat) + 0.5*(np.max(mat) - np.min(mat))

# for a in np.linspace(1, 20, 200):
a = 1
stuff = a*A
stuff = stuff - mid(stuff) + 1
stuff = 0.5*(Zy/Zx + (stuff))
stuff = stuff - mid(stuff) + 1
    # print(a, np.max(stuff) - np.min(stuff))

surf = ax3.plot_surface(TH0, TH1, A, cmap=cm.coolwarm,
                         linewidth=0, antialiased=False)



for ax in [ax11, ax12, ax21, ax22, ax3]:
    ax.set_xlabel("Zero phase")
    ax.set_ylabel("Negative phase")


plt.show()