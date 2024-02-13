import numpy as np
from matplotlib import pyplot as plt
# plt.style.use("dark_background")
from scipy.optimize import fsolve

mu = 1.25e-6
I = 400  # Amperes
L = 300  # meters
x_p = 9


def simulate_and_find_slope(h, x_p=0, plot=False):
    def determine_a_and_yd(sag):
        """
        this function determines 'a' in the equation y = h + a*cosh(z/a) such that y=h+sag when z=L/2
        """

        # using knowledge from wikipedia,
        def Hfunc(a):
            return L - 2 * a * np.arccosh((sag + a) / a)

        a = fsolve(Hfunc, np.array([100]))[0]

        # def myfun(a):
        #     return a*np.cosh(L/(2*a)) - sag - a
        #
        # mya = fsolve(myfun, np.array([100]))[0]
        # print(a, mya)

        return a, h - a

    def simulate(sag_space, x_p=0, Nz=1000, plot=False):
        B = np.zeros((3, sag_space.shape[0]))
        z = np.linspace(-L / 2, L / 2, Nz)
        # aa =[]

        for i, sag in enumerate(sag_space):
            a, y_d = determine_a_and_yd(sag)
            # aa.append(a)
            cat = lambda _z: y_d + a * np.cosh(_z / a)
            sum = np.zeros(3)
            for n in range(Nz - 1):
                z_n = z[n]
                z_np1 = z[n + 1]
                dl = np.array([0, cat(z_np1) - cat(z_n), z_np1 - z_n])
                r_prime = np.array([x_p, cat(z_n), z_n])
                length_r_prime = np.sqrt(np.sum(np.square(r_prime)))
                sum += np.cross(dl, r_prime) / (length_r_prime ** 3)
            if plot:
                plot.plot(z, cat(z))
            B[:, i] = mu * I * sum / (4 * np.pi)
        return B[0, :]

    if plot:
        fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(6, 3), layout="constrained")

    Bx_straight = -mu * I * h / (2 * np.pi * (h ** 2 + x_p ** 2))  # using simple Biot-Savart
    sag_space = np.linspace(0.1, 20, 10)
    Bx = simulate(sag_space, x_p=x_p, plot=ax1 if plot else False)
    # fig.title(f"Different levels of sag, h={h:.2f} m")
    # ax1.xlabel("z [m]")
    # ax1.ylabel("y [m]")
    # plt.savefig("Graphs\\ExamineSingleCatenaryDifferentSags.pdf", format='pdf', bbox_inches='tight')
    # quit()
    # relative difference
    y = Bx / Bx_straight
    # best fit straight line
    ((m, k), residuals, _, _, _) = np.polyfit(sag_space, y, 1, full=True)

    if plot:
        ax1.plot([0], [0], 'ko')
        ax2.plot(sag_space, y)
        ax2.plot(sag_space, sag_space * m + k)

        plt.show()

    # calculate B with sag=20, using either method:
    print(f"Simulated Bx sith sag=20m : {Bx[-1]}"
          f" with fitted line approximation: {m:.4f}*sag + {k:.2f}"
          f" difference: {abs(Bx[-1] - Bx_straight * (20 * m + k))}")
    return m, k


hspace = np.linspace(7, 15, 5)

for x_p in np.linspace(0, 10, 10):
    ms = []
    ks = []
    for h in hspace:
        m, k = simulate_and_find_slope(h, x_p=x_p)
        ms.append(m)
        ks.append(k)

    # ((m, k), residuals, _, _, _) = np.polyfit(hspace, ms, 2, full=True)
    # print(m, k)
    # plt.plot(hspace, ms)
    # plt.plot(hspace, m*hspace + k)
    # plt.show()
    # ((m, k), residuals, _, _, _) = np.polyfit(hspace, ks, 2, full=True)
    # print(m, k)
    # plt.plot(hspace, ks)
    # plt.plot(hspace, m*hspace + k)
    # plt.show()
    plt.figure("ms")
    ((a, b, c), residuals, _, _, _) = np.polyfit(hspace, ms, 2, full=True)
    print(f"m = {a}*h**2 + {b}*h + {c}")
    plt.plot(hspace, ms)
    plt.plot(hspace, a * hspace ** 2 + b * hspace + c)
    # plt.show()

    plt.figure("ks")
    ((a, b, c), residuals, _, _, _) = np.polyfit(hspace, ks, 2, full=True)
    print(f"k = {a}*h**2 + {b}*h + {c}")
    plt.plot(hspace, ks)
    plt.plot(hspace, a * hspace ** 2 + b * hspace + c)

plt.show()
