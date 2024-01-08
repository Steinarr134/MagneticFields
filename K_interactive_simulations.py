import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve, curve_fit

mu = 1.25e-6
f = 50
h = 12
N_t = 100
I0 = 500
phase_order = [-1, 0, 1]
tspace = np.arange(0, 1/f, (1/f)/N_t)

d = 0.5  # spacing between x-coils [m]
O = np.zeros(3)


# function to get sag from
def func(x):
    # return 2.09385992/np.sinh(x - 1.0043624) +30.68953421
    # return 4.66487133e+03/np.cosh(x -8.58502502e-01) -4.65016916e+03
    return 4.66487133e+03/np.cosh(x -8.58502502e-01) -4.65016916e+03

# the parameters to set
params1 = {

    "negphase": {"init": 0, "min": -180, "max": 180, "dispname": r"NP"},
    "negampl": {"init": 0, "min": 0, "max": 5, "dispname": r"NA"},
    "zerophase": {"init": 0, "min": -180, "max": 180, "dispname": r"ZP"},
    "zeroampl": {"init": 0, "min": 0, "max": 5, "dispname": r"ZA"},
    # "phi_L": {"init": 0, "min": -45, "max": 45, "dispname": r"$\phi_L$"},
    # "phi_R": {"init": 0, "min": -45, "max": 45, "dispname": r"$\phi_R$"},
    # "epsilon_L": {"init": 0, "min": -0.72, "max": 0.72, "dispname": r"$\epsilon_L$"},
    # "epsilon_R": {"init": 0, "min": -0.72, "max": 0.72, "dispname": r"$\epsilon_R$"},
    # "I_0": {"init": 500, "min": 100, "max": 1000, "dispname": r"$I_0$"},
    "Y0": {"init": 20, "min": 10, "max": 20, "dispname": r"Y0"},
    "sag": {"init": 8, "min": 2, "max": 10, "dispname": r"sag"},
    "d_P": {"init": 10, "min": 5, "max": 11, "dispname": r"$d_P$"},
    "D": {"init": 250, "min": 250, "max": 450, "dispname": r"$D$"},
}

Axes = {}


def get_sin_params(tspace, samples):
    def sin(t, A, p, w):
        return A*np.sin(w*t + p)

    p0 = [(np.max(samples) - np.min(samples))/2, np.pi, np.pi*2*f]
    eh, ehannad = curve_fit(sin, tspace, samples, p0)
    # print(eh, ehannad)
    return p0[0], eh[1]


def estimate_h(B1, B2, d_P, Y0):
    """
    This is my function that, using the magnetic field at the two locations, estimates h, the height.

    This version assumes only measuring with 2 coils, one at 45° and one pointing in the x direction
    """

    # First reduce the B1 measurement to a 45 deg coil
    B45 = np.dot(np.array([[1, 1, 0]])/np.sqrt(2), B1.T)[0]
    # and reduce second measurement to only x component
    Bx = B2[:, 0]

    # Computation method starts here.
    # B1 is the 45 deg coil measurement and B2 is the x component d meters below
    # the measurements are in sync and cover exactly 1 period of the 50 Hz signal,
    # although they are not necessarily aligned with the period

    # first, find the phase and amplitude of the two signals
    a1, p1 = get_sin_params(tspace, B45)
    a2, p2 = get_sin_params(tspace, Bx)

    phi = np.pi/2 - (p1-p2)

    h1 = (1/np.tan(phi))*d_P/np.sqrt(3) + 0.38
    h2 = Y0 - func(a2/(a1*np.sqrt(2)*np.sin(phi))) + 0.05

    print(a1*np.sqrt(2)*np.sin(phi), B1[25, 0])

    return h1, h2



def get_unbalance(plot=True):
    posampl = I0
    negampl = params1['negampl']['slider'].val*I0/100
    zerampl = params1['zeroampl']['slider'].val*I0/100
    negphase = np.deg2rad(params1['negphase']['slider'].val)
    zerphase = np.deg2rad(params1['zerophase']['slider'].val)
    a = np.exp(np.pi*2j/3)
    T = np.array([
        [1, 1, 1],
        [1, a, a**2],
        [1, a**2, a]
    ])
    phasors = np.matmul(T, np.array([[zerampl*np.exp(1j*zerphase)], [posampl], [negampl*np.exp(1j*negphase)]]))
    p_C = np.angle(phasors[0])
    A_C = np.abs(phasors[0])
    p_L = np.angle(phasors[1])
    A_L = np.abs(phasors[1])
    p_R = np.angle(phasors[2])
    A_R = np.abs(phasors[2])

    if plot:
        # Draw the stuff
        ax = Axes["Sequence"]
        ax.clear()
        colors = 'bry'
        def myplot(p1, p2, c):
            ax.plot(np.array([p1[0], p2[0]]), np.array([p1[1], p2[1]]), c)
        for p in range(3):
            myplot([-I0*1.5, 0], [-I0*1.5 + np.real(phasors[p, 0]), np.imag(phasors[p, 0])], colors[p])
            myplot([0, 0], [negampl*np.cos(negphase - p*np.pi*2/3), negampl*np.sin(negphase - p*np.pi*2/3)], colors[p])
            myplot([I0*0.5 + p, p], [I0*0.5 + p + zerampl*np.sin(zerphase), p+zerampl*np.cos(zerphase)], colors[p])

        ax.set_aspect('equal', adjustable='box')
    return p_L, p_R, p_C, A_C, A_L, A_R


def update(slidername, val):
    # get params
    p_L, p_R, p_C, A_C, A_L, A_R = get_unbalance()
    # phi_L = np.deg2rad(params1['phi_L']['slider'].val)
    # phi_R = np.deg2rad(params1['phi_R']['slider'].val)
    # epsilon_L = params1['epsilon_L']['slider'].val
    # epsilon_R = params1['epsilon_R']['slider'].val
    # I0 = params1['I_0']['slider'].val
    Y0 = params1['Y0']['slider'].val
    sag = params1['sag']['slider'].val
    d_P = params1['d_P']['slider'].val
    D = params1['D']['slider'].val
    # I0 = params1['I_0']['slider'].val

    zrange = np.linspace(0, D/2, 100)
    # Calculate the amplitude and phase of each current
    # A_L = I0*(1 + epsilon_L)
    A_C = I0
    # A_R = I0*(1 + epsilon_R)
    # p_L = phase_order[0]*np.pi*2/3 + phi_L
    p_C = 0
    # p_R = phase_order[2]*np.pi*2/3 + phi_R

    # calcualte the current in each phase
    I_L = A_L*np.sin(tspace*f*2*np.pi + p_L)
    I_C = A_C*np.sin(tspace*f*2*np.pi + p_C)
    I_R = A_R*np.sin(tspace*f*2*np.pi + p_R)

    # plot the currents
    ax_current = Axes["Current"]
    ax_current.clear()
    ax_current.plot(tspace, I_L)
    ax_current.plot(tspace, I_C)
    ax_current.plot(tspace, I_R)

    # This creates a plot of unbalanced symmetric components
    # I just found it to be very unintresting to see, so I stopped
    # showing it
    # ax_unbalanced = Axes["Unbalanced"]
    # ax_unbalanced.clear()
    # ax_unbalanced.plot([0, A_L*np.cos(p_L)], [0, A_L*np.sin(p_L)])
    # ax_unbalanced.plot([0, A_L*np.cos(p_C)], [0, A_L*np.sin(p_C)])
    # ax_unbalanced.plot([0, A_L*np.cos(p_R)], [0, A_L*np.sin(p_R)])
    # ax_unbalanced.set_aspect('equal', adjustable='box')

    """
    Now to calculate the Magnetic flux density at the measuring location (origin point)
    first solve the geometry of the catenaries
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

    # list to hold the results of this simulation
    B1 = np.zeros((N_t, 3))

    """
    B = mu*I/4pi * Integral over wire of [ dl x r_unit / r^2 ]
    """
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

    # The ellipse can now be plotted
    ax_ellipse = Axes["Ellipse"]
    ax_ellipse.clear()
    ax_ellipse.plot(*B1[:, :2].transpose())
    ax_ellipse.set_aspect('equal', adjustable='box')
    # TODO: label these axes

    """
    Do pretty muh the same thing except for a second point at a location d lower on the y axis
    """

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

    """
    Add the results of the two developed methods to the text box
    """
    # correct result
    lines = []
    h = Y0 - sag
    lines.append(f"$h = Y_0 - sag = {h:.2f} $")

    # with By/Bx method, equivalent to tan(phi) method
    B1ymax = np.max(B1[:, 1])
    B1xmax = np.max(B1[:, 0])
    h_ByBx = B1ymax*d_P/(np.sqrt(3)*B1xmax) + 0.4
    lines.append(r"$h = \frac{B_y d_P}{B_x \sqrt{3}} =$" + f"{h_ByBx:.2f}, err={h - h_ByBx:.2f}")

    # fitted function to slope of B1x/B2x
    B2xmax = np.max(B2[:, 0])
    h_funky = Y0 - func(B2xmax/B1xmax)+0.1
    lines.append(r"$h = func(\frac{B_{1|x}}{B_{2|x}}) = " + f"{h_funky:.2f}$, err={h - h_funky:.2f}")
    # print(f"{B2xmax=}\n{B1xmax=}")


    Axes["text"][1].set_text("\n".join(lines))

    # update radioplot as well if thing doesn't mess with the geomtery
    radio = Axes["Radio"][1]
    if radio.value_selected in ["ZA", "ZP", "NA", "NP"]:
        plot_as_a_function_of(radio.value_selected)



def plot_as_a_function_of(thing):
    # ax = Axes["chosen"]
    # ax.clear()
    # ax.plot(np.random.random_integers(1, 10, 20))
    # plt.draw()
    # return
    n = 21
    integral_will_stay_the_same = False
    integrals_1 = [None, None, None]
    integrals_2 = [None, None, None]
    # this plots as a function of thing
    # First creat the spaces to iterate over, all contain only the current value except for
    # thing which has a space with n values
    if thing == params1['negphase']["dispname"]:
        NP_space = np.linspace(params1['negphase']['min'], params1['negphase']['max'], n)
        chosen = NP_space
        integral_will_stay_the_same = True
    else:
        NP_space = [params1['negphase']['slider'].val]

    if thing == params1['negampl']["dispname"]:
        NA_space = np.linspace(params1['negampl']['min'], params1['negampl']['max'], n)
        chosen = NA_space
        integral_will_stay_the_same = True
    else:
        NA_space = [params1['negampl']['slider'].val]

    if thing == params1['zerophase']["dispname"]:
        ZP_space = np.linspace(params1['zerophase']['min'], params1['zerophase']['max'], n)
        chosen = ZP_space
        integral_will_stay_the_same = True
    else:
        ZP_space = [params1['zerophase']['slider'].val]

    if thing == params1['zeroampl']["dispname"]:
        ZA_space = np.linspace(params1['zeroampl']['min'], params1['zeroampl']['max'], n)
        chosen = ZA_space
        integral_will_stay_the_same = True
    else:
        ZA_space = [params1['zeroampl']['slider'].val]

    if thing == params1['Y0']["dispname"]:
        Y0_space = np.linspace(params1['Y0']['min'], params1['Y0']['max'], n)
        chosen = Y0_space
    else:
        Y0_space = [params1['Y0']['slider'].val]

    if thing == params1['sag']["dispname"]:
        sag_space = np.linspace(params1['sag']['min'], params1['sag']['max'], n)
        chosen = sag_space
    else:
        sag_space = [params1['sag']['slider'].val]

    if thing == params1['d_P']["dispname"]:
        d_P_space = np.linspace(params1['d_P']['min'], params1['d_P']['max'], n)
        chosen = d_P_space
    else:
        d_P_space = [params1['d_P']['slider'].val]

    if thing == params1['D']["dispname"]:
        D_space = np.linspace(params1['D']['min'], params1['D']['max'], n)
        chosen = D_space
    else:
        D_space = [params1['D']['slider'].val]

    # most of these for loops only have one thing that they loop over so don't freak out over the incessant looping
    h_real = []
    h_1 = []
    h_2 = []
    h_3 = []
    h_tans = []
    h_funky_tans = []
    print(thing, NP_space, NA_space, ZP_space, ZA_space, Y0_space, d_P_space, D_space, sag_space)

    for NP in NP_space:
        for NA in NA_space:
            for ZP in ZP_space:
                for ZA in ZA_space:
                    for Y0 in Y0_space:
                        for d_P in d_P_space:
                            for sag in sag_space:
                                for D in D_space:

                                    """
                                    Calculate the currents 
                                    """
                                    a = np.exp(np.pi*2j/3)
                                    T = np.array([
                                        [1, 1, 1],
                                        [1, a, a**2],
                                        [1, a**2, a]
                                    ])
                                    phasors = np.matmul(T, np.array([[(I0*(ZA/100))*np.exp(1j*np.deg2rad(ZP))], [I0], [(I0*(NA/100))*np.exp(1j*np.deg2rad(NP))]]))
                                    p_C = np.angle(phasors[0])
                                    A_C = np.abs(phasors[0])
                                    p_L = np.angle(phasors[1])
                                    A_L = np.abs(phasors[1])
                                    p_R = np.angle(phasors[2])
                                    A_R = np.abs(phasors[2])

                                    # calcualte the current in each phase
                                    I_L = A_L*np.sin(tspace*f*2*np.pi + p_L)
                                    I_C = A_C*np.sin(tspace*f*2*np.pi + p_C)
                                    I_R = A_R*np.sin(tspace*f*2*np.pi + p_R)


                                    """
                                    Now to calculate the Magnetic flux density at the measuring location (origin point)
                                    first solve the geometry of the catenaries
                                    """
                                    zrange = np.linspace(0, D/2, 100)

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
                                    # first calculate the integral for each phase, B will then scale with I as I changes over time
                                    # Each location has to be done seperately,
                                    for p_n, (P_x, I_P) in enumerate(zip([-d_P, 0, d_P], [I_L, I_C, I_R])):

                                        if integral_will_stay_the_same and integrals_1[2] is not None:
                                            integral_result = integrals_1[p_n]
                                        else:
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
                                            integrals_1[p_n] = integral_result
                                        # now that the integral has been computed the B can be simulated over time

                                        # add this to B
                                        B1 += mu*np.atleast_2d(I_P).T*integral_result/(4*np.pi)

                                    """
                                    Do pretty muh the same thing except for a second point at a location d lower on the y axis
                                    """

                                    B2 = np.zeros((N_t, 3))
                                    # first calculate the integral for each phase, B will then scale with I as I changes over time
                                    for p_n, (P_x, I_P) in enumerate(zip([-d_P, 0, d_P], [I_L, I_C, I_R])):

                                        if integral_will_stay_the_same and integrals_2[2] is not None:
                                            integral_result = integrals_2[p_n]
                                        else:
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
                                            integrals_2[p_n] = integral_result
                                        # now that the integral has been computed the B can be simulated over time

                                        # add this to B
                                        B2 += mu*np.atleast_2d(I_P).T*integral_result/(4*np.pi)

                                    """
                                    Add the results of the two developed methods to the text box
                                    """
                                    # correct result
                                    lines = []
                                    h = Y0 - sag
                                    # lines.append(f"$h = Y_0 - sag = {h:.2f} $")
                                    h_real.append(h)

                                    h_tan, h_funky_tan = estimate_h(B1, B2, d_P, Y0)
                                    h_tans.append(h_tan)
                                    h_funky_tans.append(h_funky_tan)
                                    # # with By/Bx method, equivalent to tan(phi) method
                                    B1ymax = np.max(B1[:, 1])
                                    B1xmax = np.max(B1[:, 0])
                                    h_ByBx = B1ymax*d_P/(np.sqrt(3)*B1xmax) + 0.4
                                    # lines.append(r"$h = \frac{B_y d_P}{B_x \sqrt{3}} =$" + f"{h_ByBx:.2f}, err={h - h_ByBx:.2f}")
                                    h_1.append(h_ByBx)

                                    # fitted function to slope of B1x/B2x
                                    B2xmax = np.max(B2[:, 0])
                                    h_funky = Y0 - func(B2xmax/B1xmax) + 0.1
                                    # lines.append(r"$h = func(\frac{B_{1|x}}{B_{2|x}}) = " + f"{h_funky:.2f}$, err={h - h_funky:.2f}")
                                    h_2.append(h_funky)
                                    # print(f"{B2xmax=}\n{B1xmax=}")

                                    # shape of ellipse regardless of tilt
                                    ellipse = np.linalg.norm(B1, axis=1)
                                    h_ellipse = np.max(ellipse)*d_P/(np.min(ellipse)*np.sqrt(3))+0.4
                                    h_3.append(h_ellipse)
    print("chosen plot update")
    ax = Axes["chosen"]
    ax.clear()
    h_1 = np.array(h_1); h_2 = np.array(h_2); h_real=np.array(h_real); h_3 = np.array(h_3)
    ax.plot(chosen, h_real, 'y')
    ax.plot(chosen, h_1, 'b')
    ax.plot(chosen, h_2, 'r')
    h_avg = (h_1 + h_2)/2
    ax.plot(chosen, h_avg)
    ax.plot(chosen, h_3, 'g')
    ax.set_ylabel("h[m]")
    ax.set_xlabel(thing)
    ax.legend(["real",
               f"By/Bx | maxerr={max(abs(h_real-h_1)):.2f}",
               f"func(B1/B2) | maxerr={max(abs(h_real-h_2)):.2f}",
               f"avg | maxerr={max(abs(h_real-h_avg)):.2f}"])
    plt.draw()



