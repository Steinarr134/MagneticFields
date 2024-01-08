import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, fsolve

mu = 4*np.pi*1e-7
Imax = 500
d = 9
f = 50
freq = f
h = 12
N_t = 100
Y0 = 20
D = 400
phase_order = [-1, 0, 1]

phase_xs = [-d, 0, d]
hspace = np.arange(5, Y0-3, 1)
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
    # print(eh, ehannad)
    return p0[0], eh[1]


# Bx = np.zeros((len(hspace), 2))
# By = np.zeros((len(hspace), 2))
# B45 = np.zeros((len(hspace), 2))

O = np.zeros(3)
zrange = np.arange(0, D/2, 1)
# for h in hspace:
# for Imax in np.arange(200, 500, 100):
print(f"{Imax=}")
bs = []

Ds = np.linspace(300, 500, 100)
for D in Ds:  # probably second order coupling

# Ys = np.linspace(20, 22, 10)
# for Y0 in Ys:  # linearly coupled to b

# for mu in np.linspace(mu*0.9, mu*1.1, 5):  # uncoupled

# ds = np.linspace(d*0.8, d*1.2, 10)
# for d in ds:  # linearly coupled
    # break
    phis = np.zeros(hspace.shape)
    for j, h in enumerate(hspace):
        sag = Y0-h
        print("sag: ", sag)
        """
        First work out the catenary curve describing the wires
        
        catenary curve is:
        y = a*cosh(z/a) + b  ,                (y, z)
        just need to find a and b to satisfy (Y0, D/2) and (Y0-sag, 0)
        cosh(0) = 1 so second boundary condition gives Y0-sag = a+b
        
        second boundary is not useful on it's own
        instead use H = 2*a*arcosh((h+a)/a) where H is tower spacing and h is sag
        """
        # define Hfunc which has to be solved to find a and b
        def Hfunc(a):
            return D - 2*a*np.arccosh((sag + a)/a)

        # solve Hfunc to find a and b
        a = fsolve(Hfunc, np.array([100]))[0]
        b = a - (Y0 - sag)

        # now define the function  y = catenary(z)
        def catenary(z):
            return a*np.cosh(z/a) - b

        # list to hold the results of this simulation (with this sag)
        Bs = np.zeros((N_t, 3))

        """
        B = mu*I/4pi * Integral over wire of [ dl x r_unit / r^2 ]
        """
        # first calculate the integral for each phase, B will then scale with I as I changes over time

        for P_x, p in zip(phase_xs, phase_order):

            integral_result = np.zeros(3)
            # starting point of curve
            last_P = np.array([P_x, catenary(0), 0])
            for z in zrange[1:]:
                next_P = np.array([P_x, catenary(z), z])
                dl = next_P - last_P

                OP = last_P - O
                OP_length = np.sqrt(np.sum(np.power(OP, 2)))
                OP_unit = OP/OP_length
                integral_result += np.cross(dl, OP_unit)/(OP_length**2)
                last_P = next_P
            # since integral only uses half the wire the other half is same except mirrored over XY plane
            # thus the Z component drops out but X and Y component doubles
            integral_result = 2*integral_result
            integral_result[2] = 0
            # now that the integral has been computed the B can be simulated over time

            for i, t in enumerate(tspace):

                I = Imax*np.sin(t*freq*2*np.pi +p*(120*np.pi*2/360))
                # add this to B
                Bs[i, :] += mu*I*integral_result/(4*np.pi)

        Bs = np.array(Bs)
        # get phase offset of 45 deg
        B45 = np.dot(Bs[:, :2], np.sqrt(2)*np.array([1, 1])/2)
        # plt.plot(tspace, B45, label="B45")
        B45A, B45P = get_sin_params(tspace, B45)
        phis[j] = -B45P

    # plt.plot(hspace, np.tan(phis)*d/np.sqrt(3))
    # plt.show()

    y = np.tan(phis)*d/np.sqrt(3)
    a, b = np.polyfit(hspace, y, 1)
    bs.append(b)

    # plt.plot(hspace, y, '.', label='y = ' + '{:.2e}'.format(b) + ' + {:.5f}'.format(a) + 'x')
    # plt.plot([0, hspace[-1]], [0, hspace[-1]])
    # print(a,b)
    # plt.xlabel("h [m] - height of wires")
    # plt.ylabel("tan(phi)*d/sqrt(3)")
    # plt.text(np.min(hspace), np.max(y), 'y = ' + '{:.2e}'.format(b) + ' + {:.5f}'.format(a) + 'x', size=14)
    # plt.title("Confirm h = tan(phi)*d/sqrt(3)")
    # plt.show()

# plt.plot(tspace, B[:, 0], label="x-component")
# plt.plot(tspace, B[:, 1], label="y-component")
# B45 = np.dot(B[:, :2], np.sqrt(2)*np.array([1, 1])/2)
# plt.plot(tspace, B45, label="B45")
# B45A, B45P = get_sin_params(tspace, B45)
# print(B45A, np.rad2deg(B45P))
# plt.plot(tspace, B45A*np.sin(tspace*np.pi*2*50 + B45P))
# plt.legend()
# plt.show()

plt.plot(Ds, bs)
print(np.polyfit(Ds, bs, 2, full=True))
print(repr(Ds))
print(repr(bs))
plt.show()


"""
By above it is clear that
b(D) = array([ 1.89988166e-07, -1.82173924e-04, -2.36870574e-01]   (3)
b(Y) = (array([-0.00104394, -0.25848545])                           (2)
b(d) = (array([-3.10402753e-02,  2.01125325e-15])                     (1)

since b(d) = d*(   ) + 0 it is clear that all terms of b include d and then some combination of 1, Y, D and D^2
so we can assume 

b(d, D, Y) = d*(a0 + a1*Y + a2*D + a3*Y*D + a4*D^2 + a5*Y*D^2)
(1) then gives 
(a0 + a1*Y0 + a2*D0 + a3*Y0*D0 + a4*D0^2 + a5*Y0*D0^2) = -3.10402753e-02

rearranging gives
b(d, D, Y) = Y*d*(a1 + a3*D + a5*D^2) + d*(a0 + a2*D + a4*D^2)
thus

(a1 + a3*D0 + a5*D0^2) = -0.00104394/d0  
(a0 + a2*D + a4*D^2) = -0.25848545/d0

now
b(d, D, Y) = D^2*d*(a4 + a5*Y) + D*d*(a2 + a3*Y) + d*(a1+ a1*Y)
thus
(a4 + a5*Y0) =  1.89988166e-07/d0
(a2 + a3*Y) = -1.82173924e-04/d0
(a0+ a1*Y) = -2.36870574e-01/d0

those are 6 equations for 6 unknowns. in matrix form it is

A*a = b
where a = [a0, a1, a2, a3, a4, a5]'
and:
"""
Y0 = 20
D0 = 400
d0 = 9
A = np.array([
    [1, Y0, D0, Y0*D0, D0**2, Y0*D0**2],
    [0, 1, 0, D0, 0, D0**2],
    [1, 0, D0, 0, D0**2, 0],
    [0, 0, 0, 0, 1, Y0],
    [0, 0, 1, Y0, 0, 0],
    [1, Y0, 0, 0, 0, 0]
])
b = np.array([
    [-3.10402753e-02],
    [-0.00104394/d0],
    [-0.25848545/d0],
    [1.89988166e-07/d0],
    [-1.82173924e-04/d0],
    [-2.36870574e-01/d0],
])

print(A, b)
a = np.linalg.solve(A, b)
print(a)
