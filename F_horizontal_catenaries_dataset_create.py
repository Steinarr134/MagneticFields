import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

mu = 1.5e-6
Phase_offset_x = 10  # phase spacing
Y0 = 20  # starting height of wires at tower
freq = 50
Phase1_x = -Phase_offset_x
Phase2_x = 0

Phase3_x = Phase_offset_x
Ps_x = [Phase1_x, Phase2_x, Phase3_x]
O1 = np.zeros(3)  # Origin point  (where measurement is taken)
O2 = np.array([0, -0.5, 0])  # second measuring point
D = 250  # distance between towers

# sags of interest
# sags = np.arange(5, Y0-7, 0.5)
sags = np.array([10])
# Currents
# Is = np.arange(200, 500, 25)
Is = np.array([500])
# data to be saved
DataB1s = np.zeros((Is.shape[0], sags.shape[0]))
DataB2s = np.zeros((Is.shape[0], sags.shape[0]))

# range of z to integrate over wire
zrange = np.arange(0, D/2, 1)
for i, I in enumerate(Is):
    maxBs = []
    ground_clearances = []
    for j, sag in enumerate(sags):
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

        B1s = []
        B2s = []
        for t in np.arange(0, 1/freq, 1/(freq*50)):
            I_1 = I*np.sin(t*freq*2*np.pi - 120*np.pi*2/360)
            I_2 = I*np.sin(t*freq*2*np.pi)
            I_3 = I*np.sin(t*freq*2*np.pi + 120*np.pi*2/360)

            # print(t, I_2)

            B1 = np.zeros(3)  # Magnetic field at measuring point (origin)
            B2 = np.zeros(3)  # Magnetic field strength d below origin

            """
            B = mu*I/4pi * Integral over wire of [ dl x r_unit / r^2 ]
            """
            for P_x, I_n in zip(Ps_x, [I_1, I_2, I_3]):

                integral_result_1 = np.zeros(3)
                integral_result_2 = np.zeros(3)
                # starting point of curve
                last_P = np.array([P_x, catenary(0), 0])
                for z in zrange[1:]:
                    next_P = np.array([P_x, catenary(z), z])
                    dl = next_P - last_P

                    OP = last_P - O1
                    OP_length = np.sqrt(np.sum(np.power(OP, 2)))
                    OP_unit = OP/OP_length
                    integral_result_1 += np.cross(dl, OP_unit)/(OP_length**2)


                    OP = last_P - O2
                    OP_length = np.sqrt(np.sum(np.power(OP, 2)))
                    OP_unit = OP/OP_length
                    integral_result_2 += np.cross(dl, OP_unit)/(OP_length**2)

                    last_P = next_P
                    # print(np.cross(dl, OP_unit)/(OP_length**2))

                # since integral only uses half the wire the other half is same except mirrored over XY plane
                # thus the Z component drops out but X and Y component doubles
                integral_result_1 = 2*integral_result_1
                integral_result_1[2] = 0
                # print(integral_result)
                B1 += mu*I_n*integral_result_1/(4*np.pi)

                integral_result_2 = 2*integral_result_2
                integral_result_2[2] = 0
                # print(integral_result)
                B2 += mu*I_n*integral_result_2/(4*np.pi)

            B1s.append(B1[0])
            B2s.append(B2[0])

        # plt.plot(Bs)
        # plt.show()
        maxB1 = np.max(B1s)
        DataB1s[i, j] = maxB1
        maxB2 = np.max(B2s)
        DataB2s[i, j] = maxB2

# import pickle
# with open("D250_Yn20_z.pkl", 'wb') as f:
#     pickle.dump((Is, sags, Y0, D, DataB1s, DataB2s),f)

print(DataB1s, DataB2s)



#
# plt.plot(ground_clearances, maxBs)
# plt.title("Ground Clearance vs Magnetic strength \n Catenary curves, horizontally arranged phases \n I=" + str(I) + " A")
# plt.xlabel("Ground clearance [m]")
# plt.ylabel("Magnetic field strength [micro T]")
# plt.grid(1)
# # plt.show()
#
#
# # fit a function of the form y=A*x^(-b)
#
# from scipy.optimize import curve_fit
#
#
# def func(x, A, b):
#     return A*x**(-b)
#
#
# popt, pcov = curve_fit(func, ground_clearances, maxBs)#, bounds=((-np.inf, 1), (np.inf, 2)))
# print(popt, pcov)
#
# plt.plot(ground_clearances, func(ground_clearances, *popt))
#
# plt.show()