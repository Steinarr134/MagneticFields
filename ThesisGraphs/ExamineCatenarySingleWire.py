import numpy as np
from matplotlib import pyplot as plt
# plt.style.use("dark_background")
from scipy.optimize import fsolve

mu = 1.25e-6
I = 400  # Amperes
L = 400  # meters
h = 8  # meters


def determine_a_and_yd(sag):
    """
    this function determines 'a' in the equation y = h + a*cosh(z/a) such that y=h+sag when z=L/2
    """
    # using knowledge from wikipedia,
    def Hfunc(a):
        return L - 2*a*np.arccosh((sag + a)/a)
    a = fsolve(Hfunc, np.array([100]))[0]

    # def myfun(a):
    #     return a*np.cosh(L/(2*a)) - sag - a
    #
    # mya = fsolve(myfun, np.array([100]))[0]
    # print(a, mya)

    return a, h-a


def simulate(sag_space, Nz=1000):
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
        plt.plot(z, cat(z))
        B[:, i] = mu*I*sum/(4*np.pi)
    return B[0, :]


Bx_straight = -mu*I/(2*np.pi*h)  # using simple Biot-Savart
sag_space = np.linspace(0.1, 20, 100)
plt.subplots(figsize=(6, 3), layout="constrained")
Bx = simulate(sag_space)

plt.plot([0], [0], 'ko')
plt.title("Different levels of sag")
plt.xlabel("z [m]")
plt.ylabel("y [m]")
# plt.savefig("Graphs\\ExamineSingleCatenaryDifferentSags.pdf", format='pdf', bbox_inches='tight')
# plt.show()
# quit()
# relative difference
y = Bx/Bx_straight

# best fit straight line
((m, k), residuals, _, _, _) = np.polyfit(sag_space, y, 1, full=True)

# calculate B with sag=20, using either method:
print(f"Simulated Bx sith sag=20m : {Bx[-1]}"
      f"\n with fitted line approximation: {Bx_straight*(20*m  + k)}")

# calculate R^2 value of the fit
y_hat = m*sag_space + k
y_bar = y.mean()
ss_tot = ((y-y_bar)**2).sum()
ss_res = ((y-y_hat)**2).sum()
R2 = 1 - (ss_res/ss_tot)

s = f"r = {m:.5e}*sag + {k:.5f}, \nR^2 = {R2:.8f}"

plt.figure(figsize=(6, 4))
plt.plot(sag_space, y)
plt.legend([s])
plt.xlabel("Sag [m]")
plt.ylabel(r"Relative change in $\mathbf{B}$")
# plt.savefig("Graphs\\ExamineSingleCatenaryLineFit.pdf", format='pdf', bbox_inches='tight')
# plt.show()

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

