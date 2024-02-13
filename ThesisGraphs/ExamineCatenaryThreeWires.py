import numpy as np
from matplotlib import pyplot as plt
# plt.style.use("dark_background")
from scipy.optimize import fsolve

"""
This program simulates the magnetic flux density at a location under three wires
that hang in a catenary shape. The current is assumed to be three perfect phasors

The sag is varied from practically nothing to a lot

The output are plots:

A - showing the strength of the x and y components as a function of sag
B - showing the ratio between the x and y components as a function of sag
C - showing the ratio of the two x components as a function of sag
"""

mu = 1.25e-6
I_rms = 400  # Amperes
I0 = I_rms*np.sqrt(2)
L_default = 300  # meters
H_default = 10  # meters
x_p_default = 5  # meters
f = 50  # Hertz
r = np.array([0, 0, 0])
d = 0.5  # [m], spacing between measuring locations

# The values of sag to be considered
N_sag = 50
hspace = np.linspace(5, H_default-0.5, N_sag)
sag_space = H_default - hspace

# Timespace
N_t = 500
t_space = np.linspace(0, 1/f, N_t)

# The current in each phase (right, center and left)
I_L = I0*np.sin(2*np.pi*f*t_space - 2*np.pi/3)
I_C = I0*np.sin(2*np.pi*f*t_space)
I_R = I0*np.sin(2*np.pi*f*t_space + 2*np.pi/3)


def rms(signal):
    return np.sqrt(np.mean(np.square(signal)))


def determine_a_and_yd(sag, L=L_default, H=H_default, x_p=x_p_default):
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


from joblib import Memory
memory = Memory("cache", verbose=0)

@memory.cache
def simulate(sag_space, L=L_default, H=H_default, x_p=x_p_default, Nz=1000, plot=False):
    # arrays to hold the results
    Bx = np.zeros((sag_space.shape[0]))
    By = np.zeros((sag_space.shape[0]))
    Bx2 = np.zeros((sag_space.shape[0]))

    # divide the z dimension down into smaller increments
    z = np.linspace(-L/2, L/2, Nz)

    # Loop through the values of sag to simulate
    for i, sag in enumerate(sag_space):

        # determine the 'a' constant and create the catenary function
        a, y_d = determine_a_and_yd(sag, L, H, x_p)
        cat = lambda _z: y_d + a*np.cosh(_z/a)

        # optionally plot the catenary curve
        if plot:
            plt.plot(z, cat(z))

        # Array to hold the magnetic flux density
        B1 = np.zeros((N_t, 3))
        B2 = np.zeros((N_t, 3))

        for x, I_P in zip([-x_p, 0, x_p], [I_L, I_C, I_R]):
            # Array to hold the result of the integral of  (dl X r')/|r'|^3
            integralsum = np.zeros(3)
            # Calculate the integral by summing over z axis
            for n in range(Nz-1):
                # z_n = z[n]
                # z_np1 = z[n+1]
                # dl = np.array([0, cat(z_np1) - cat(z_n), z_np1-z_n])
                # r_prime = np.array([x, cat(z_n) , z_n])
                # length_r_prime = np.sqrt(np.sum(np.square(r_prime)))
                # integralsum += np.cross(dl, r_prime)/length_r_prime**3

                z_n = np.array([x, cat(z[n]), z[n]])
                z_np1 = np.array([x, cat(z[n+1]), z[n+1]])
                dl = z_np1 - z_n
                r_prime = (z_n + z_np1)/2 - r
                # r_prime = z_n - r
                length_r_prime = np.sqrt(np.sum(np.square(r_prime)))
                r_prime_unit = r_prime/length_r_prime
                integralsum += np.cross(dl, r_prime_unit)/length_r_prime**2
            # Now with the integral solved the x, y, z components of
            # the magnetic flux density caused by this phase wire can be calculated
            # add the contribution of this phase to total.
            B1 += (mu/(4*np.pi)) * np.outer(I_P, integralsum)

            # Redo the entire thing for second measuring location,
            # d meters below the first one
            r2 = r - np.array([0, d, 0])

            # Array to hold the result of the integral of  (dl X r')/|r'|^3
            integralsum = np.zeros(3)
            # Calculate the integral by summing over z axis
            for n in range(Nz-1):
                z_n = np.array([x, cat(z[n]), z[n]])
                z_np1 = np.array([x, cat(z[n+1]), z[n+1]])
                dl = z_np1 - z_n
                r_prime = (z_n + z_np1)/2 - r2
                length_r_prime = np.sqrt(np.sum(np.square(r_prime)))
                r_prime_unit = r_prime/length_r_prime
                integralsum += np.cross(dl, r_prime_unit)/length_r_prime**2
            # Now with the integral solved the x, y, z components of
            # the magnetic flux density caused by this phase wire can be calculated
            # add the contribution of this phase to total.
            B2 += (mu/(4*np.pi)) * np.outer(I_P, integralsum)


        # B1 and B2 now hold the magnetic flux density in x, y, z, directions over time
        # Our interest is in the rms of the x and y components
        Bx[i] = rms(B1[:, 0])
        By[i] = rms(B1[:, 1])
        Bx2[i] = rms(B2[:, 0])

    # finally return the values of interest
    return Bx, By, Bx2


"""
First plot shows different levels of sag and the x1, x2 and y components as a function of sag
"""
# set up subplot figure
plt.subplots(1, 2, figsize=(9, 4), layout="constrained")
plt.subplot(1, 2, 1)
# simulate every 10th value and plot (plots the wire shapes)
simulate(sag_space[::8], plot=True)

# simulate all the values
Bx, By, Bx2 = simulate(sag_space)

Bx = Bx*1e6
Bx2 = Bx2*1e6
By = By*1e6

plt.plot([0], [0], 'ko')
plt.plot([0], [-d], 'ko')
plt.title("Different levels of sag")
plt.xlabel("z [m]")
plt.ylabel("y [m]")

# Bx = Bx*1e-6
# By = By*1e-6

# plt.figure(figsize=(6, 4))
plt.subplot(1, 2, 2)
plt.title("Magnetic Flux Density vs Sag")
plt.plot(sag_space, Bx, label='x-component')
plt.plot(sag_space, By, label='y-component')
plt.plot(sag_space, Bx2, label='x-component, second location')
# plt.plot(sag_space, y2, label='Straight line')
# plt.plot(np.log(sag_space), y2)
plt.xlabel("Sag [m]")
plt.ylabel(r"$\mathbf{B} [\mu T]$")
plt.legend()
# plt.savefig("Graphs\\CatenaryThreeWiresBalancedPhases.pdf", format='pdf', bbox_inches='tight')


"""
Second figure shows the error in the Bx/By * d_p/sqrt(3) as a function of sag, Nz=1000
"""
plt.figure()
y = By/Bx * x_p_default/np.sqrt(3)
x = H_default  - sag_space
(a, b) = np.polyfit(x, y, 1)
plt.plot(x, y-x, label="$N_z= 1000$")#label=r"$f = \frac{B_x d_p}{B_y \sqrt{3}}$")
plt.xlabel("sag [m]")
plt.ylabel("$\\frac{B_x d_p}{B_y \sqrt{3}} - h$")

"""
Third figure shows the error as a function of Nz for one value of sag
It levels out at around Nz = 1000, thus it's unneccesary to go higher.
"""
# plt.figure()
# errors = []
# n_z_s = np.logspace(2, 5, 40)
# for n_z in n_z_s:
#     print(n_z)
#     Bx, By, Bx2 = simulate(np.array([5]), Nz=int(n_z))
#     y = By/Bx * x_p_default /np.sqrt(3)
#     x = H_default  - sag_space
#     errors.append(y-x)
# plt.loglog(n_z_s, errors)
# plt.xlabel("$log(N_z)$")
# plt.ylabel("error")

"""
Fourth figure shows the distribution of errors for a wide range of input parameters
"""
errors = []

# n = 10
# for D in np.linspace(200, 400, n):
#     for d_p in np.linspace(5, 10, n):
#         for Y0 in np.linspace(10, 20, n):
#             sag_space = np.linspace(0.5, Y0-5, n)
#             Bx, By, Bx2 = simulate(sag_space, L=D, x_p=d_p, H=Y0)
#             errors.append( By/Bx * x_p_default /np.sqrt(3))

# it takes way to long to go through it like this so instead just loop through random points
np.random.seed(1234)
N = 10000

def rand(bounds):
    return bounds[0] + np.random.random()*(bounds[1]-bounds[0])

for i in range(N):
    Y0 = rand((10, 20))
    Bx, By, Bx2 = simulate(np.array([rand((1, Y0-5))]),
                           L=rand((200, 400)),
                           x_p=rand((5, 10)),
                           H=Y0)
    errors.append(By/Bx * x_p_default/np.sqrt(3))
    if i % 100 == 0:
        print(i)

stuff = np.array(errors).flatten()
plt.figure()
plt.hist(stuff)

plt.show()
quit()