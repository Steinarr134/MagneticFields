import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

mu = 1.25e-6
I = 500  # Current in wire
Phase_offset_y = 10  # phase spacing
Y0 = 20  # starting height of wires at tower
freq = 50
Phase1_y = Y0
O = np.zeros(3)  # Origin point  (where measurement is taken)
D = 250  # distance between towers

# stuff to plot later
B_catenaries = []
B_straights = []
ground_clearances = []

# sags of interest
sags = np.arange(5, Y0-7, 0.5)
# range of z to integrate over wire
zrange = np.arange(0, D/2, 0.1)
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

    B = np.zeros(3)  # Magnetic field at measuring point (origin)

    """
    B = mu*I/4pi * Integral over wire of [ dl x r_unit / r^2 ]
    """

    integral_result = np.zeros(3)
    # starting point of curve
    last_P = np.array([0, catenary(0), 0])
    for z in zrange[1:]:
        next_P = np.array([0, catenary(z), z])
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
    B = mu*I*integral_result/(4*np.pi)

    B_catenaries.append(-B[0] * 1e6)
    ground_clearances.append(Y0-sag)

    # if using biot-savart solution for straight line
    B_straights.append(1e6*mu*I/(2*np.pi*(Y0-sag)))


plt.plot(ground_clearances, B_catenaries)
plt.plot(ground_clearances, B_straights)
# print(repr(ground_clearances), repr(maxBs))
plt.title(f"Ground Clearance vs Magnetic strength \n Catenary curves, vertically arranged phases \n I={I} A")
plt.xlabel("Ground clearance [m]")
plt.ylabel("Magnetic field strength [micro T]")
plt.grid(1)

plt.figure()
ratios = np.array(B_straights)/np.array(B_catenaries)
offset = np.array(B_straights) - np.array(B_catenaries)
(zs, residual, _, _, _) = np.polyfit(ground_clearances, offset, 2, full=True)
print(residual)
plt.plot(ground_clearances, offset)
plt.plot(ground_clearances, np.poly1d(zs)(ground_clearances))

plt.show()



