import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve, curve_fit

"""
This is a simulation of one specific sag and current combination  of horizontally spaced phases
but looks at the measurement inaccuracy of the coil caused by the field strength being curved
"""

mu = 1.25e-6
I = 500  # Current in wire
Phase_offset_x = 10  # phase spacing
Y0 = 20  # starting height of wires at tower
freq = 50
Pspacing = 10
Phase1_x = -Phase_offset_x
Phase2_x = 0
Phase3_x = Phase_offset_x
Ps_x = [Phase1_x, Phase2_x, Phase3_x]
D = 250  # distance between towers
R = 0.17/2 # radius of coil in m
# R = 1
R2 = R**2


# stuff to plot later
Bs = []
Ds = []

# sags of interest
# sags = np.arange(5, Y0-7, 0.5)
sag = 15
# range of z to integrate over wire
zrange = np.arange(0, D/2, 1)
yrange = np.arange(-5, 9, 0.25)
# for sag in sags:
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


# def integral_of_half_circle_between(x1, x2):
#     return 0.5*x1*np.sqrt(R**2 - x1**2) + 0.5*R**2*np.arctan(x1/np.sqrt(R**2 - x**2)) \
#      - (0.5*x2*np.sqrt(R**2 - x2**2) + 0.5*R**2*np.arctan(x2/np.sqrt(R**2 - x**2)))
#
def integral_of_half_circle_between(x1, x2):
    return 0.5*np.sqrt(R2 - x2**2)*x2 + 0.5*R2*np.arcsin(x2/R) - (0.5*np.sqrt(R2 - x1**2)*x1 + 0.5*R2*np.arcsin(x1/R))

def is_inside_circle(x_, y_):
    return (x_**2 + y_**2) < R2

def rect_area_inside_circle(center, stepsize):
    # always use positive quadrant
    center = np.abs(center)
    # since it's symmetric always use x<=y
    center = np.sort(center)
    x1 = center[0] - stepsize/2
    x2 = center[0] + stepsize/2
    y1 = center[1] - stepsize/2
    y2 = center[1] + stepsize/2
    if not is_inside_circle(x1, y1): # completely outside circle
        return 0
    if is_inside_circle(x2, y2): # completely inside circle
        return stepsize**2

    """
    now Let's say we have this rectangle
  y2 ________________________
    |a      \                |             
    |         \              |             
    |           \            |             
    |            \           |             
    |             \          |             
    |              \         |             
    |               \       b|               
  y1 ------------------------
    x1                       x2
    
    we already know that (x1, y1) is inside circle and (x2, y2) is outside
    that leaves 4 possibilities depending on (x1, y2) and (x2, y1),
    both inside, both outside, and either one inside but other outside
    let's call the corners in question a and b
    
    simplest case is a outside, b inside, then it is the integral of half circle between x1 and x2
    and then subtract th rectangle between x1 and x2 with height y1
    
    also simple is a inside b outside (pictured above) in that case it is about flipping the axis and
    using same logic, it's the integral of half circle between y1 and y2 minus the rectangle
    
    both outside is trickier, then find the intersection of half circle and x=x1 , call it y_i
    the area will be the integarl of half circle from y1 to y_i minus x1*(y_i-y1)
    
    both inside is then simple just do the integral minus rectangle below minus the rectangle above
    """
    a = is_inside_circle(x1, y2)
    b = is_inside_circle(x2, y1)

    # print(a,b)
    if a:
        if b:
            # print(a,b)
            return integral_of_half_circle_between(x1, x2) - y1*stepsize \
                - rect_area_inside_circle(np.array(center) + np.array([0, stepsize]), stepsize)
        else:
            #
            return integral_of_half_circle_between(y1, y2) - x1*stepsize
    else:
        if b:
            return integral_of_half_circle_between(x1, x2) - y1*stepsize
        else:
            y_i = np.sqrt(R**2 - x1**2)
            # print(y_i)
            return integral_of_half_circle_between(y1, y_i) - (y_i-y1)*x1


def calculate_B_at(O):
    # for t in np.arange(0, 1/freq, 1/(freq*50)):
    I_1 = I*np.cos(- 120*np.pi*2/360)
    I_2 = I*np.cos(0)
    I_3 = I*np.cos(120*np.pi*2/360)

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
    return B[0]


"""
Here I need to integrate over the area of the coil, so need to split a circle into smaller areas and loop over those
but how to split a circle? The whole thing is symmetric over the y-axis so I can use that
perhaps just split into rectangles and multiply the field strength by the area inside the circle?
"""

stepsize = R/32
grid = np.arange(-R+stepsize/2, R, stepsize)
from itertools import product
flux = 0
for x, y in product(grid, grid):

    """ now to calculate how much of the rectangle is inside the coil"""
    A = rect_area_inside_circle((x, y), stepsize)

    B = calculate_B_at((x, y, 0))

    flux += B*A

print(f"Flux based on area of coil and B at (0, 0, 0): {calculate_B_at((0, 0, 0))*R2*np.pi}")
print(f"Flux based on simulation of field curvature: {flux}")

"""
results with I = 500, sag = 15 (excessive), R=0.17/2 and stepsize = R/32:
Flux based on area of coil and B at (0, 0, 0): -3.9870481908076406e-07
Flux based on simulation of field curvature: -3.987097859512047e-07

negligible difference!
"""

