from HyParI_curve_fit import *
from scipy.optimize import fsolve


mu = 1.25e-6
f = 50
N_t = 1000
I0 = 500
phase_order = [-1, 0, 1]
tspace = np.arange(0, 1/f, (1/f)/N_t)

d = 0.5  # spacing between x-coils [m]
O = np.zeros(3)
from joblib import Memory
memory = Memory("cache", verbose=0)

@memory.cache
def F(sagspace, D, Y0=12, d_P=7, d=0.5):
    """
    This function will calculate the ratio between B1 and B2
    those being the amplitudes of the x components of the magnetic field strength
    at two locations, spaced d apart.
    """
    print("sdkfhj")

    zrange = np.linspace(0, D/2, 100)
    I_C = I0
    I_L = np.sin(np.deg2rad(120))*I_C
    I_R = np.sin(np.deg2rad(-120))*I_C


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


if __name__ == '__main__':

    D_bounds = (250, 400)
    d_P_bounds = (5, 8)
    Y0_bounds = (12, 15)
    # d_bounds = ()

    sag_space = np.linspace(2, 8, 20)


    def func(x, a, b, c, d):
        # return a/np.cosh(x+c) + b
        # return a*x**3 + b*x**2 + c*x + d
        return a + b*x**-1 + c*x**-2 + d*x**-3
    #
    Z = n_hyperparameters(F, F_m=func, hyper_order=4,
                          bounds=[D_bounds, Y0_bounds]#, d_P_bounds ]
                          , xspace=sag_space, n=20)

    print(repr(Z))
    print(get_ndimquation(Z, letters=["sag", "D", "Y0", "d_P"]))

