"""
Ok, I-ve learned a lot from writing L2!

I think I have a good process:

first create functions that estimate h from xx and xy based on Y0, D and d_P.

next loop through phase angles and use least squares to find an optimum way to combine the xx and xy results
into on master results.

Here I will write a function that does this
It will accept Y0, D and d_P and it will find xx and xy functions
and finally return a function that guesses h based on the xy and xx ratios.

The main section will then test if it works

It works! but only if Y0/d_P < 2.2
"""
import time

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.optimize import least_squares, fsolve

plt.style.use("dark_background")
mu = 1.25e-6
f = 50
I0 = 500
O = np.zeros((3))
d = 0.5
N_t = 100



def rms(signal):
    return np.sqrt(np.mean(np.square(signal)))


def calculate_ratio(h, Y0, D, d_P, phi0, phi1, e0=0.05, e1=0.02):
    zrange = np.linspace(0, D/2, 100)
    """
    Calculate the currents 
    """
    a = np.exp(np.pi*2j/3)
    T = np.array([
        [1, 1, 1],
        [1, a, a**2],
        [1, a**2, a]
    ])
    phasors = np.matmul(T, np.array([[(I0*e0)*np.exp(1j*np.deg2rad(phi0))], [I0], [(I0*e1)*np.exp(1j*np.deg2rad(phi1))]]))
    p_C = np.angle(phasors[0])
    A_C = np.abs(phasors[0])
    p_L = np.angle(phasors[1])
    A_L = np.abs(phasors[1])
    p_R = np.angle(phasors[2])
    A_R = np.abs(phasors[2])

    tspace = np.linspace(0, 1/f, N_t)
    # calcualte the current in each phase
    I_L = A_L*np.sin(tspace*f*2*np.pi + p_L)
    I_C = A_C*np.sin(tspace*f*2*np.pi + p_C)
    I_R = A_R*np.sin(tspace*f*2*np.pi + p_R)

    ret = []
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
    B1 = np.zeros((N_t, 3))

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
    B2 = np.zeros((N_t, 3))
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

    return np.max(B2[:, 0])/np.max(B1[:, 0]), np.max(B1[:,1])/np.max(B1[:, 0])


def make_xx_fun(Y0, D, d_P):
    def calc_ratio(h):
        zrange = np.linspace(0, D/2, 100)
        I_C = I0
        I_L = np.sin(np.deg2rad(90 + 120))*I_C
        I_R = np.sin(np.deg2rad(90 - 120))*I_C

        ret = []
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

        return rms(B2[:, 0])/rms(B1[:, 0])
    ratios = []
    hspace = np.linspace(5, Y0-1)
    for h in hspace:
        ratios.append(calc_ratio(h))
    zs = np.polyfit(ratios, hspace, 4)
    fun = np.poly1d(zs)
    # plt.plot(ratios, hspace)
    # mi = np.min(ratios)
    # ma = np.max(ratios)
    # ratiospace = np.linspace(mi - 0.5*(ma - mi), ma + 0.5*(ma - mi))
    # plt.plot(ratiospace, fun(ratiospace), '--')
    # plt.show()
    return fun


def make_xy_fun(Y0, D, d_P):
    def calc_ratio(h):
        Nz = 1001
        zrange = np.linspace(-D/2, D/2, Nz)
        tspace = np.linspace(0, 1/f, N_t)
        I_C = np.sin(2*np.pi*f*tspace)*I0
        I_L = np.sin(2*np.pi*f*tspace + np.deg2rad(120))*I0
        I_R = np.sin(2*np.pi*f*tspace + np.deg2rad(-120))*I0

        ret = []
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
        B1 = np.zeros((N_t,3))

        """
        B = mu*I/4pi * Integral over wire of [ dl x r_unit / r^2 ]
        """

        # might be able to make these calculations faster by skipping the current and stuff
        # that cancels out when taking the ration B1/B2 however,
        # first calculate the integral for each phase, B will then scale with I as I changes over time
        for p_n, (P_x, I_P) in enumerate(zip([-d_P, 0, d_P], [I_L, I_C, I_R])):

            integral_result = np.zeros(3)
            # starting point of curve
            last_P = np.array([P_x, catenary(zrange[0]), zrange[0]])
            for z in zrange[1:]:
                next_P = np.array([P_x, catenary(z), z])
                dl = next_P - last_P
                OP = last_P - O
                OP_length = np.sqrt(np.sum(np.power(OP, 2)))
                OP_unit = OP/OP_length
                integral_result += np.cross(dl, OP_unit)/(OP_length**2)
                last_P = next_P

            # for n in range(Nz-1):
            #     z_n = np.array([P_x, catenary(zrange[n]), zrange[n]])
            #     z_np1 = np.array([P_x, catenary(zrange[n+1]), zrange[n+1]])
            #     dl = z_np1 - z_n
            #     r_prime = (z_n + z_np1)/2
            #     length_r_prime = np.sqrt(np.sum(np.square(r_prime)))
            #     r_prime_unit = r_prime/length_r_prime
            #     integral_result += np.cross(dl, r_prime_unit)/length_r_prime**2

            # since integral only uses half the wire the other half is same except mirrored over XY plane
            # thus the Z component drops out but X and Y component doubles
            # integral_result = 2*integral_result
            # integral_result[2] = 0
            # now that the integral has been computed the B can be simulated over time

            # add this to B
            B1 += mu*np.atleast_2d(I_P).T*integral_result/(4*np.pi)

        return rms(B1[:, 1])/rms(B1[:, 0])
    offsets = []
    for h in np.linspace(5, Y0-1):
        offsets.append(calc_ratio(h)*d_P/np.sqrt(3) - h)

    # plt.plot(offsets)
    # plt.show()
    print("xy average offset: ", float(np.average(offsets)))

    def fun(ratio):
        return ratio*d_P/np.sqrt(3) - float(np.average(offsets))
    return fun


def get_master_function(Y0, D, d_P):

    print("finding xx and xy functions based on perfect balance")
    xx_fun = make_xx_fun(Y0, D, d_P)
    xy_fun = make_xy_fun(Y0, D, d_P)
    n = 10
    z_space = np.linspace(-180, 180, n)
    n_space = np.linspace(-180, 180, n)
    h_space = np.linspace(5, Y0-0.5, n)
    print("simulating all kinds of unbalances..")
    if True:
        xx = np.zeros((n, n, n))
        xy = np.zeros((n, n, n))
        correct = np.zeros((n, n, n))
        for k, h in enumerate(h_space):
            print(f"{k/n:.0%} done")
            for i, phi0 in enumerate(z_space):
                for j, phi1 in enumerate(n_space):
                    xx_ratio, xy_ratio = calculate_ratio(h, Y0, D, d_P, phi0, phi1)
                    xx_h_guess = xx_fun(xx_ratio)
                    xy_h_guess = xy_fun(xy_ratio)
                    # print(Y0, D, d_P, phi0, phi1, h, xx_h_guess, xy_h_guess)
                    xx[k, i, j] = xx_h_guess
                    xy[k, i, j] = xy_h_guess
                    correct[k, i, j] = h
        np.save("cache/xx.npy", xx)
        np.save("cache/xy.npy", xy)
        np.save("cache/correct.npy", correct)
    else:
        xx = np.load("cache/xx.npy")
        xy = np.load("cache/xy.npy")
        correct = np.load("cache/correct.npy")

    def fun(x):
        ans = correct - (x[0]*xx + x[1]*xy + x[2]*xx**2 + x[3]*xy**2 + x[4]*xx*xy)
        return ans.flatten()

    results = least_squares(fun, np.ones(5)/5)
    x = results["x"]
    i, j, k = np.unravel_index(np.argmax(np.abs(fun(x)).reshape((10, 10, 10))), (10, 10, 10))
    print(f"masterfun: biggest err = {np.max(np.abs(fun(x))):.2f} @ (h={h_space[i]},phi0={z_space[j]},phi1={n_space[k]})")

    def retfun(xx_ratio, xy_ratio):
        h_xx = xx_fun(xx_ratio)
        h_xy = xy_fun(xy_ratio)
        return x[0]*h_xx + x[1]*h_xy + x[2]*h_xx**2 + x[3]*h_xy**2 + x[4]*h_xx*h_xy
    return retfun, np.max(np.abs(fun(x)))

if __name__ == '__main__':

    sag_bounds = (2, 8)
    D_bounds = (200, 400)
    d_P_bounds = (5, 8)
    Y0_bounds = (10, 15)
    phi_bounds = (-180, 180)

    def rand(bounds):
        return bounds[0] + np.random.random()*(bounds[1]-bounds[0])

    # np.random.seed(2)
    # D = rand(D_bounds)
    # Y0 = rand(Y0_bounds)
    # d_P = rand(d_P_bounds)

    n = 3
    Y0 = 15
    D = 300
    Y0s = np.linspace(*Y0_bounds, n)
    d_Ps = np.linspace(*d_P_bounds, n)
    # Ds = np.linspace(*D_bounds, n)
    # ratios = np.linspace(2.2, 2.5, 10)
    # d_Ps = Y0/ratios
    maxerrs = np.zeros((n, n))
    for i, Y0 in enumerate(Y0s):
        for j, d_P in enumerate(d_Ps):
            print(f"{Y0=:.1f}\t{D=:.1f}\t{d_P=:.1f}\t getting master function")
            mf, maxerr = get_master_function(Y0, D, d_P)

            """
            now test the master function by randomly choosing sag, and phase angles and see if the rsults are good
            """

            err = []
            ts = []
            # print("Master function found, testing it by throwing random sag and unbalance at it")
            # for _ in range(200):
            #     sag = rand((1, Y0-5))
            #     phi0 = rand(phi_bounds)
            #     phi1 = rand(phi_bounds)
            #     h = Y0 - sag
            #     xx_ratio, xy_ratio = calculate_ratio(h, Y0, D, d_P, phi0, phi1)
            #     guess = mf(xx_ratio, xy_ratio)
            #     err.append(h-guess)

            # """
            # problem happens at 180, 180 so simulate there for different sags
            # """
            #
            # for sag in np.linspace(1, Y0-5):
            #     h = Y0 - sag
            #     xx_ratio, xy_ratio = calculate_ratio(h, Y0, D, d_P, 180, 180)
            #     guess = mf(xx_ratio, xy_ratio)
            #     err.append(h-guess)
            # print(err)

            # print(f"max={maxerr:.2f}")
            # plt.stairs(*np.histogram(err))
            # plt.show()
            maxerrs[i, j] = maxerr

    print(repr(maxerrs))
