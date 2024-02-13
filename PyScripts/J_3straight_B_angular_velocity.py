import numpy as np
from matplotlib import pyplot as plt

mu = 4*np.pi*1e-7
Imax = 500
d = 15

h = 10

phase_xs = [-d, 0, d]
hspace = np.arange(5, 20, 1)
tspace = np.arange(0, 0.02, 0.0005)


def get_cos_params(samples):
    N = len(samples)
    x = np.linspace(-np.pi, np.pi, N, endpoint=False)
    template = np.exp(1j * x)
    corr = 2 / N * template@samples
    R = np.abs(corr)
    phi = np.log(corr).imag
    return R, phi/(2*np.pi)


Bx = np.zeros((len(hspace), 2))
By = np.zeros((len(hspace), 2))
B45 = np.zeros((len(hspace), 2))


for h in hspace:
    # for j, h in enumerate(hspace):
    Ihat = np.array([0, 0, 1])
    B = np.zeros((tspace.shape[0], 3))
    for i, t in enumerate(tspace):
        for p, phase_x in enumerate(phase_xs):
            r = np.array([h, phase_x, 0])
            I = Imax*np.sin(2*np.pi*50*t + (p-1)*2*np.pi/3)*Ihat

            B[i, :] += (mu*np.cross(I, r))/(2*np.pi*np.sum(np.square(r)))

    # theta is angle of magnetic field
    theta = np.arctan2(B[:, 0], B[:, 1])
    # fix wrap around
    theta[np.argmin(theta):] += np.pi * 2
    # plt.title("Theta(t)")
    # plt.plot(tspace, w)
    plt.plot(tspace, np.gradient(theta, tspace), label=str(h))
    plt.legend()
plt.show()
