import numpy as np
from matplotlib import pyplot as plt
# plt.style.use("dark_background")
from scipy.optimize import fsolve

mu = 1.25e-6
I = 400  # Amperes
L = 400  # meters
H = 20  # meters
x_p = 9


def determine_a_and_yd(sag):
    """
    this function determines 'a' in the equation y = h + a*cosh(z/a) such that y=h+sag when z=L/2
    """
    # using knowledge from wikipedia,
    # def Hfunc(a):
    #     return L - 2*a*np.arccosh((sag + a)/a)
    # a = fsolve(Hfunc, np.array([100]))[0]

    def myfun(a):
        return a*np.cosh(L/(2*a)) + ((H - sag) - a) - H

    a = fsolve(myfun, np.array([100]))[0]
    # print(a, mya)

    return a, H - sag - a


def simulate(sag_space, Nz=1000, plot=False):
    B = np.zeros((3, sag_space.shape[0]))
    z = np.linspace(-L/2, L/2, Nz)

    for i, sag in enumerate(sag_space):
        a, y_d = determine_a_and_yd(sag)
        cat = lambda _z: y_d + a*np.cosh(_z/a)
        sum = np.zeros(3)
        for n in range(Nz-1):
            z_n = z[n]
            z_np1 = z[n+1]
            dl = np.array([0, cat(z_np1) - cat(z_n), z_np1-z_n])
            r_prime = np.array([0, cat(z_n), z_n])
            length_r_prime = np.sqrt(np.sum(np.square(r_prime)))
            sum += np.cross(dl, r_prime)/length_r_prime**3
        if plot:
            plt.plot(z, cat(z))
        B[:, i] = mu*I*sum/(4*np.pi)
    return B[0, :]

sag_space = np.linspace(4, 13, 100)
Bx_straight = []
for sag in sag_space:
    m = -6.371284917248984e-05*(sag) - 0.0002679703055899428
    # k = -0.0002780037137790345*(sag) + 1.0014311437466543
    k = 1
    # h = sag

    # m = 1.2058807974775086e-06*h**2 + -9.024222671699517e-05*h + -0.00012991928540464497
    # k = -1.240472197618192e-05*h**2 + -5.099830303110974e-06*h + 1.0000110327974554
    Bx_straight.append((m*(H-sag) + k)*mu*I/(2*np.pi*(H-sag)))  # using simple Biot-Savart
plt.subplots(1, 2, figsize=(9, 4), layout="constrained")
plt.subplot(1, 2, 1)
simulate(sag_space[::10], plot=True)
Bx = simulate(sag_space)

plt.plot([0], [0], 'ko')
plt.title("Different levels of sag")
plt.xlabel("z [m]")
plt.ylabel("y [m]")
plt.savefig("Graphs\\ExamineSingleCatenaryDifferentSags.pdf", format='pdf', bbox_inches='tight')
# plt.show()
# quit()
# relative difference
y = -Bx*1e6
y2 = np.array(Bx_straight)*1e6


# plt.figure(figsize=(6, 4))
plt.subplot(1, 2, 2)
plt.plot(sag_space, y)
# plt.plot(np.log(sag_space), y2)
plt.xlabel("Sag [m]")
plt.ylabel(r"$\mathbf{B} [\mu T]$")
plt.savefig("Graphs\\ExamineSingleCatenaryFixedSuspension.pdf", format='pdf', bbox_inches='tight')
# plt.show()

# rel = y/y2
# zs = np.polyfit(H - sag_space, y/y2, 2)
# fun = np.poly1d(zs)

plt.figure()
plt.plot(H-sag_space, y/y2)
# plt.plot(H-sag_space, fun(H-sag
# _space))
plt.show()
quit()
plt.figure()
# do it again for higher values of sag to see the full picture
sag_space = np.linspace(0.1, 2000, 200)
Bx = simulate(sag_space)

# relative difference
y = Bx/Bx_straight


plt.figure()
plt.plot(sag_space, y)
plt.xlabel("Sag [m]")
plt.ylabel(r"Relative change in $\mathbf{B}$")
# plt.savefig("Graphs\\ExamineSingleCatenaryFullRelationship.pdf", format='pdf', bbox_inches='tight')
plt.show()

