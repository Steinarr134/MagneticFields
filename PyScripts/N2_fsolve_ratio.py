"""
Goal is to write a function that solves for h given a ratio and hyperparameters d_P, Y0, D

Rough structure of the function:

Use some precomputed values to create an initial guess for h

Simulate and use newton-s algorithm (or just fsolve) to refine the guess for h

-----
It works

"""
import numpy as np
from scipy.optimize import least_squares, fsolve

mu = 1.25e-6
f = 50
I0 = 500
O = np.zeros((3))

def find_h(ratio, Y0, D, d_P, d=0.5):

    def calc_ratio(h):
        zrange = np.linspace(0, D/2, 100)
        I_C = I0
        I_L = np.sin(np.deg2rad(90 +120))*I_C
        I_R = np.sin(np.deg2rad(90-120))*I_C

        ret = []
        sag = Y0 - h[0]

        # First work out the geometries of the catenary curve
        # define Hfunc which has to be solved to find a and b
        def Hfunc(a):
            return D - 2*a*np.arccosh((sag + a)/a)

        # solve Hfunc to find a and b
        a = fsolve(Hfunc, np.array([100]))[0]
        b = a - (Y0 - sag)

        # now define the function  y = catenary(z)
        def catenary(z):
            return a*np.cosh(z/a) - b

        # list to hold the results of this simulation
        B1 = np.zeros((1,3))

        """
        B = mu*I/4pi * Integral over wire of [ dl x r_unit / r^2 ]
        """

        # might be able to make these calculations faster by skipping the current and stuff
        # that cancels out when taking the ration B1/B2 however,
        # first calculate the integral for each phase, B will then scale with I as I changes over time
        for p_n, (P_x, I_P) in enumerate(zip([-d_P, 0, d_P], [I_L, I_C, I_R])):

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

            # add this to B
            B1 += mu*np.atleast_2d(I_P).T*integral_result/(4*np.pi)
        B2 = np.zeros((1,3))
        # first calculate the integral for each phase, B will then scale with I as I changes over time
        for p_n, (P_x, I_P) in enumerate(zip([-d_P, 0, d_P], [I_L, I_C, I_R])):

            integral_result = np.zeros(3)
            # starting point of curve
            last_P = np.array([P_x, catenary(0), 0])
            for z in zrange[1:]:
                next_P = np.array([P_x, catenary(z), z])
                dl = next_P - last_P
                OP = last_P - (O - np.array([0, d, 0]))
                OP_length = np.sqrt(np.sum(np.power(OP, 2)))
                OP_unit = OP/OP_length
                integral_result += np.cross(dl, OP_unit)/(OP_length**2)
                last_P = next_P
            # since integral only uses half the wire the other half is same except mirrored over XY plane
            # thus the Z component drops out but X and Y component doubles
            integral_result = 2*integral_result
            integral_result[2] = 0
            # now that the integral has been computed the B can be simulated over time

            # add this to B
            B2 += mu*np.atleast_2d(I_P).T*integral_result/(4*np.pi)

        return ratio - B2[0, 0]/B1[0, 0]

    h0 = 10

    result = least_squares(calc_ratio, [h0], bounds=((3), (20)))
    return result

if __name__ == '__main__':

    def F(sagspace, D, Y0=12, d_P=7, d=0.5):
        """
        This function will calculate the ratio between B1 and B2
        those being the amplitudes of the x components of the magnetic field strength
        at two locations, spaced d apart.
        """
        print("sdkfhj")

        zrange = np.linspace(0, D/2, 100)
        I_C = I0
        I_L = np.sin(np.deg2rad(90 +120))*I_C
        I_R = np.sin(np.deg2rad(90-120))*I_C


        ret = []
        hspace = Y0-sagspace
        for h in hspace:
            sag = Y0 - h

            # First work out the geometries of the catenary curve
            # define Hfunc which has to be solved to find a and b
            def Hfunc(a):
                return D - 2*a*np.arccosh((sag + a)/a)

            # solve Hfunc to find a and b
            a = fsolve(Hfunc, np.array([100]))[0]
            b = a - (Y0 - sag)

            # now define the function  y = catenary(z)
            def catenary(z):
                return a*np.cosh(z/a) - b

            # list to hold the results of this simulation
            B1 = np.zeros((1,3))

            """
            B = mu*I/4pi * Integral over wire of [ dl x r_unit / r^2 ]
            """

            # might be able to make these calculations faster by skipping the current and stuff
            # that cancels out when taking the ration B1/B2 however,
            # first calculate the integral for each phase, B will then scale with I as I changes over time
            for p_n, (P_x, I_P) in enumerate(zip([-d_P, 0, d_P], [I_L, I_C, I_R])):

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

                # add this to B
                B1 += mu*np.atleast_2d(I_P).T*integral_result/(4*np.pi)
            B2 = np.zeros((1,3))
            # first calculate the integral for each phase, B will then scale with I as I changes over time
            for p_n, (P_x, I_P) in enumerate(zip([-d_P, 0, d_P], [I_L, I_C, I_R])):

                integral_result = np.zeros(3)
                # starting point of curve
                last_P = np.array([P_x, catenary(0), 0])
                for z in zrange[1:]:
                    next_P = np.array([P_x, catenary(z), z])
                    dl = next_P - last_P
                    OP = last_P - (O - np.array([0, d, 0]))
                    OP_length = np.sqrt(np.sum(np.power(OP, 2)))
                    OP_unit = OP/OP_length
                    integral_result += np.cross(dl, OP_unit)/(OP_length**2)
                    last_P = next_P
                # since integral only uses half the wire the other half is same except mirrored over XY plane
                # thus the Z component drops out but X and Y component doubles
                integral_result = 2*integral_result
                integral_result[2] = 0
                # now that the integral has been computed the B can be simulated over time

                # add this to B
                B2 += mu*np.atleast_2d(I_P).T*integral_result/(4*np.pi)

            ret.append(B2[0, 0]/B1[0, 0])
        # print(np.array(ret))
        return np.array(ret)
    sagspace = np.linspace(2, 8)
    D_bounds = (200, 400)
    d_P_bounds = (5, 8)
    Y0_bounds = (10, 15)

    def rand(bounds):
        return bounds[0] + np.random.random()*(bounds[1]-bounds[0])

    # from M1_hyperparameter_research import F
    print(F(np.array([5]), 250, 15, 7))

    for i in range(10):
        D = rand(D_bounds)
        Y0 = rand(Y0_bounds)
        d_P = rand(d_P_bounds)
        sag = rand(sagspace)

        h = Y0-sag
        ratio = F(np.array([sag]),D, Y0, d_P)

        print( Y0, D, d_P, h, find_h(ratio, Y0, D, d_P)["x"])

