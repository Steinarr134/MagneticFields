import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

mu = 4*np.pi*1e-7
Imax = 500
d = 9
f = 50
h = 12
N_t = 100
N_h = 100

phase_xs = [-d, 0, d]
hspace = np.linspace(5, 20, N_h)
tspace = np.linspace(0, 1/f, N_t, endpoint=False)



Bx = np.zeros((len(hspace)))
By = np.zeros((len(hspace)))
# B45 = np.zeros((len(hspace), 2))


# for h in hspace:
phis = np.zeros(hspace.shape)
for j, h in enumerate(hspace):
    Ihat = np.array([0, 0, 1])
    B = np.zeros((tspace.shape[0], 3))
    for i, t in enumerate(tspace):
        for p, phase_x in enumerate(phase_xs):
            r = np.array([phase_x, h, 0])
            I = Imax*np.sin(2*np.pi*f*t + (p-1)*2*np.pi/3)*Ihat

            B[i, :] += (mu*np.cross(I, r))/(2*np.pi*np.sum(np.square(r)))
    Bx[j] = np.max(B[:, 0])
    By[j] = np.max(B[:, 1])


y = (By/Bx)*(d/np.sqrt(3))
plt.plot(hspace, y, '.')
a, b = np.polyfit(hspace, y, 1)
print(a,b)
plt.xlabel("h [m] - height of wires")
plt.ylabel("By*d/(Bx*sqrt(3))")
plt.text(np.min(hspace), np.max(y), 'y = ' + '{:.2e}'.format(b) + ' + {:.5f}'.format(a) + 'x', size=14)
plt.title("Confirm By/Bx = sqrt(3)*h/d")
plt.show()
