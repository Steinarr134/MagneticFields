import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

mu = 4*np.pi*1e-7
Imax = 500
d = 9
f = 50
h = 12
N_t = 100

phase_xs = [-d, 0, d]
hspace = np.arange(5, 20, 1)
tspace = np.arange(0, 1/f, (1/f)/N_t)


def get_cos_params(x, samples):
    N = len(samples)
    # x = np.linspace(0, 2*np.pi, N, endpoint=False)
    template = np.exp(1j * x)
    corr = 2 / N * template@samples
    R = np.abs(corr)
    phi = np.log(corr).imag
    return R, phi/(2*np.pi)

def get_sin_params(tspace, samples):
    def sin(t, A, p, w):
        return A*np.sin(w*t + p)

    p0 = [(np.max(samples) - np.min(samples))/2, np.pi, np.pi*2*f]
    eh, ehannad = curve_fit(sin, tspace, samples, p0)
    print(eh, ehannad)
    return p0[0], eh[1]


# Bx = np.zeros((len(hspace), 2))
# By = np.zeros((len(hspace), 2))
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

    # get phase offset of 45 deg
    B45 = np.dot(B[:, :2], np.sqrt(2)*np.array([1, 1])/2)
    # plt.plot(tspace, B45, label="B45")
    B45A, B45P = get_sin_params(tspace, B45)
    phis[j] = -B45P

plt.plot(hspace, np.tan(phis)*d/np.sqrt(3))
plt.show()

y = np.tan(phis)*d/np.sqrt(3)
plt.plot(hspace, y, '.')
a, b = np.polyfit(hspace, y, 1)
print(a,b)
plt.xlabel("h [m] - height of wires")
plt.ylabel("tan(phi)*d/sqrt(3)")
plt.text(np.min(hspace), np.max(y), 'y = ' + '{:.2e}'.format(b) + ' + {:.5f}'.format(a) + 'x', size=14)
plt.title("Confirm h = tan(phi)*d/sqrt(3)")
plt.show()

# plt.plot(tspace, B[:, 0], label="x-component")
# plt.plot(tspace, B[:, 1], label="y-component")
# B45 = np.dot(B[:, :2], np.sqrt(2)*np.array([1, 1])/2)
# plt.plot(tspace, B45, label="B45")
# B45A, B45P = get_sin_params(tspace, B45)
# print(B45A, np.rad2deg(B45P))
# plt.plot(tspace, B45A*np.sin(tspace*np.pi*2*50 + B45P))
# plt.legend()
# plt.show()



