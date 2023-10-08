import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

mu = 1.25e-6
I = 500  # Current in wire
Phase_offset_x = 10  # phase spacing
Y0 = 20  # starting height of wires at tower
freq = 50
Phase1_x = -Phase_offset_x
Phase2_x = 0
Phase3_x = Phase_offset_x
Ps_x = [Phase1_x, Phase2_x, Phase3_x]
O = np.zeros(3)  # Origin point  (where measurement is taken)
D = 250  # distance between towers

# stuff to plot later
maxBs = []
ground_clearances = []

# sags of interest
sags = np.arange(5, Y0-7, 0.5)
# range of z to integrate over wire
zrange = np.arange(0, D/2, 1)
for sag in sags:
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
    def Hfunc(a):
        return D - 2*a*np.arccosh((sag + a)/a)
    a = fsolve(Hfunc, np.array([100]))[0]
    b = a - (Y0 - sag)

    def catenary(z):
        return a*np.cosh(z/a) - b

    Bs = []
    for t in np.arange(0, 1/freq, 1/(freq*50)):
        I_1 = I*np.sin(t*freq*2*np.pi - 120*np.pi*2/360)
        I_2 = I*np.sin(t*freq*2*np.pi)
        I_3 = I*np.sin(t*freq*2*np.pi + 120*np.pi*2/360)

        # print(t, I_2)

        B = np.zeros(3)  # Magnetic field at measuring point (origin)

        """
        B = mu*I/4pi * Integral over wire of [ dl x r_unit / r^2 ]
        """
        for P_x, I_n in zip(Ps_x, [I_1, I_2, I_3]):

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
                # print(np.cross(dl, OP_unit)/(OP_length**2))

            # since integral only uses half the wire the other half is same except mirrored over XY plane
            # thus the Z component drops out but X and Y component doubles
            integral_result = 2*integral_result
            integral_result[2] = 0
            # print(integral_result)
            B += mu*I_n*integral_result/(4*np.pi)

        Bs.append(B[0])

    # plt.plot(Bs)
    # plt.show()
    maxB = np.max(Bs)
    maxBs.append(maxB*1e6)
    ground_clearances.append(Y0-sag)



plt.plot(ground_clearances, maxBs)
plt.title("Ground Clearance vs Magnetic strength \n Catenary curves, horizontally arranged phases \n I=" + str(I) + " A")
plt.xlabel("Ground clearance [m]")
plt.ylabel("Magnetic field strength [micro T]")
plt.grid(1)
# plt.show()


# fit a function of the form y=A*x^(-b)

from scipy.optimize import curve_fit


def func(x, A, b):
    return A*x**(-b)


popt, pcov = curve_fit(func, ground_clearances, maxBs)#, bounds=((-np.inf, 1), (np.inf, 2)))
print(popt, pcov)

plt.plot(ground_clearances, func(ground_clearances, *popt))

plt.show()