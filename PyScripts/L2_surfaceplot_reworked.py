"""
Goal is to use the results from N2_fsolve_ratio

then loop through phi0 and phi1
for each, calculate the ratio and then work back to find h
then do similar to what I was trying in L1, find a good central point


"""

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



def find_h_xx(ratio, Y0, D, d_P, d=0.5):

    def calc_ratio(h):
        zrange = np.linspace(0, D/2, 100)
        I_C = I0
        I_L = np.sin(np.deg2rad(90 + 120))*I_C
        I_R = np.sin(np.deg2rad(90 - 120))*I_C

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

    h0 = (3 + Y0)/2

    h = least_squares(calc_ratio, [h0], bounds=((3), (Y0-0.5)))
    if abs(h["x"] - (Y0-0.5)) < 0.01:
        print(h)
    return h



def find_h_xy(ratio, Y0, D, d_P, d=0.5):

    def calc_ratio(h):
        zrange = np.linspace(0, D/2, 100)
        tspace = np.linspace(0, 1/f, N_t)
        I_C = np.sin(2*np.pi*f*tspace)*I0
        I_L = np.sin(2*np.pi*f*tspace + np.deg2rad(120))*I0
        I_R = np.sin(2*np.pi*f*tspace + np.deg2rad(-120))*I0

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

        # I don't understand what is happening here but I meant to do this:
        ans = ratio - np.max(B1[:, 1])/np.max(B1[:, 0])
        # so the upper level function could find the zero point
        # however I accidentally did this:
        # ans = B1[:, 1]/B1[:, 0]
        # which returns a vector where the instantaneous ratios over time get returned
        # and somehow it works perfectly, i.e. it is not dependent on the current unbalance
        # nevermind, it wasn't actually working, it just found a nic h value regardless of the ratio
        # so it didn't change with sag or anything.
        # print(np.sum(np.square(ans)))
        # print(h, ans)
        # print(h, ans)
        return ans

    h0 = (3 + Y0)/2

    h = least_squares(calc_ratio, [h0], bounds=((3), (Y0-0.5)))
    # if h["cost"] > 0.0001:
    #     print(h)
    return h


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

        return B2[0, 0]/B1[0, 0]
    ratios = []
    hspace = np.linspace(5, Y0-1)
    for h in hspace:
        ratios.append(calc_ratio(h))
    plt.plot(ratios, hspace)
    zs = np.polyfit(ratios, hspace, 4)
    fun = np.poly1d(zs)
    mi = np.min(ratios)
    ma = np.max(ratios)
    ratiospace = np.linspace(mi - 0.5*(ma - mi), ma + 0.5*(ma - mi))
    plt.plot(ratiospace, fun(ratiospace), '--')
    plt.show()
    return fun

def make_xy_fun(Y0, D, d_P):
    def calc_ratio(h):
        zrange = np.linspace(0, D/2, 100)
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

        return np.max(B1[:, 1])/np.max(B1[:, 0])
    offsets = []
    for h in np.linspace(5, Y0-1):
        offsets.append(calc_ratio(h)*d_P/np.sqrt(3) - h)

    def fun(ratio):
        return ratio*d_P/np.sqrt(3) - int(np.average(offsets))
    return fun



sagspace = np.linspace(2, 8)
D_bounds = (200, 400)
d_P_bounds = (5, 8)
Y0_bounds = (10, 15)
phi_bounds = (-180, 180)

def rand(bounds):
    return bounds[0] + np.random.random()*(bounds[1]-bounds[0])
if __name__ == '__main__':

    from M1_hyperparameter_research import F

    # print(calculate_ratio(10, 15, 250, 7, 0, 0, 0, 0))
    phi1 = 0
    phi0 = 0
    # np.random.seed(42)
    # D = rand(D_bounds)
    # Y0 = rand(Y0_bounds)
    # d_P = rand(d_P_bounds)
    # for i in range(10):
    #     sag = rand(sagspace)
    #     for j in range(5):
    #         phi0 = rand(phi_bounds)
    #         phi1 = rand(phi_bounds)
    #
    #         h = Y0-sag
    #
    #         xx_ratio, xy_ratio = calculate_ratio(h, Y0, D, d_P, phi0, phi1)
    #         # xx_h_guess = find_h_xx(xx_ratio, Y0, D, d_P)["x"]
    #         print(xy_ratio, Y0, D, d_P, phi0, phi1, h)
    #         xy_h_guess = find_h_xy(xy_ratio, Y0, D, d_P)["x"]
    # #
    #         print(Y0, D, d_P, phi0, phi1, h, xy_h_guess, h/xy_h_guess)
    # quit()
    n = 10
    z_space = np.linspace(-180, -150, n)
    n_space = np.linspace(-180, 180, n)
    TH0, TH1 = np.meshgrid(z_space, n_space)
    D = 300
    Y0 = 14
    d_P = 7
    sag = 5

    xx = np.zeros((n, n))
    xy = np.zeros((n, n))

    np.random.seed(1)
    D = rand(D_bounds)
    Y0 = rand(Y0_bounds)
    d_P = rand(d_P_bounds)
    sag = rand(sagspace)
    print(f"{Y0=:.2f}\t{D=:.2f}\t{d_P=:.2f}\t")

    h = Y0-sag

    xy_fun = make_xy_fun(Y0, D, d_P)
    xx_fun = make_xx_fun(Y0, D, d_P)
    # from joblib import Memory
    # memory = Memory("cache", verbose=0)
    # @memory.cache
    def do_the_stuff(phi0, phi1):
        # xx_h_guess = find_h_xx(xx_ratio, Y0, D, d_P)
        # xy_h_guess = find_h_xy(xy_ratio, Y0, D, d_P)
        return xx_h_guess, xy_h_guess

    for i, phi0 in enumerate(z_space):
        for j, phi1 in enumerate(n_space):
            xx_ratio, xy_ratio = calculate_ratio(h, Y0, D, d_P, phi0, phi1)
            xy_h_guess = xy_fun(xy_ratio)
            xx_h_guess = xx_fun(xx_ratio)
            print( phi0, phi1, xx_ratio, xy_ratio, h, xx_h_guess, xy_h_guess)
            xx[i, j] = xx_h_guess
            xy[i, j] = xy_h_guess


    print(repr(xx))
    print(repr(xy))

    fig2, (ax21, ax22) = plt.subplots(ncols=2, subplot_kw={"projection": "3d"})

    surf1 = ax21.plot_surface(TH0, TH1, xx.T-h, cmap=cm.coolwarm,
                              linewidth=0, antialiased=False)
    surf2 = ax22.plot_surface(TH0, TH1, xy.T-h, cmap=cm.coolwarm,
                              linewidth=0, antialiased=False)
    fig3, ax3 = plt.subplots(subplot_kw={"projection": "3d"})

    # def fun(x):
    #     ans = h - (x[0]*xx + x[1]*xy + x[2]*xx**2 + x[3]*xy**2 + x[4]*xx*xy)
    #     return ans.flatten()
    #
    # results = least_squares(fun, np.ones(5)/5)
    # print(results["cost"])
    # print(results["x"])
    #
    # x = results["x"]
    # A = (x[0]*xx + x[1]*xy + x[2]*xx**2 + x[3]*xy**2 + x[4]*xx*xy)
    A = 0.5*xx + 0.5*xy
    # A = -0.27275753*xxErrors + 1.49172883*xyErrors + 0.29223873*xxErrors**2 - 0.52753101*xyErrors**2
    # A = 0.8972132 * xx + 1.04751483 * xy - 0.20115701 * xx ** 2 + -0.27398632 * xy ** 2 - 0.47011142 * xx * xy
    # A = 0.8972132 * xx + 1.04751484 * xy -0.02235078 * xx ** 2 + -0.03044292 * xy ** 2 -0.0522346 * xx * xy

    surf = ax3.plot_surface(TH0, TH1, A-h, cmap=cm.coolwarm,
                            linewidth=0, antialiased=False)

    for ax in [ax21, ax22, ax3]:
        ax.set_xlabel("Zero phase")
        ax.set_ylabel("Negative phase")
    plt.show()
