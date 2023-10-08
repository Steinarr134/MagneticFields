import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

mu = 1.25e-6
I = 440  # Current in wire
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
Bcurves = []
ground_clearances = []

# sags of interest
sags = np.arange(5, Y0-7, 0.5)
# sags = [5]
# range of z to integrate over wire
zrange = np.arange(0, D/2, 1)
alphas = np.arange(0, np.pi/16, np.pi/160)
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


    """
    function that finds the (x,y,z) coordinates of the wire given the z coordinate
    and alpha, galloping angle
    """
    def wire_coordinates(z, alpha):
        y0 = catenary(z)  # find y coordinate given no galloping
        s = Y0 - y0  # sag at this location
        return np.array([s*np.sin(alpha), y0 - np.cos(alpha), z])



    """
    B = mu*I/4pi * Integral over wire of [ dl x r_unit / r^2 ]
    """
    Bs = []
    # loop over galloping angles
    for alpha in alphas:
        integral_result = np.zeros(3)
        # starting point of curve
        last_P = wire_coordinates(0, alpha)
        for z in zrange[1:]:
            next_P = wire_coordinates(z, alpha)
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
        Bs.append(1e6*mu*I*integral_result[0]/(4*np.pi))
    Bcurves.append(Bs)
    # plt.plot(Bs)
    # plt.show()
    # maxB = np.max(Bs)
    # maxBs.append(maxB*1e6)
    ground_clearances.append(Y0-sag)

for Bcurve in Bcurves:
    plt.plot(np.rad2deg(alphas), Bcurve)
plt.title("Single wire catenary curve galloping \n different levels of ground clearance")
plt.xlabel("Galloping [°]")
plt.ylabel("Magnetic field strength [micro T]")
plt.legend(Y0 - sags)
plt.grid(1)
plt.show()

