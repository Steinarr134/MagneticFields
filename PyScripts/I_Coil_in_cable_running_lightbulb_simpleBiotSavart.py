import numpy as np
from matplotlib import pyplot as plt

V_rms = 220  # rms voltage
Lightbulb_wattage = 100
I_rms = Lightbulb_wattage/V_rms
mu = 1.25e-6

wire_y = 0.2  # the two wires are at plus and minus 0.2 m


# use biot-savart simplified solution assuming infinite, straight wires
def biot_savart(r):
    return mu*I_rms/(2*np.pi*r)


# construct grid
y_space = np.linspace(-2*wire_y, 2*wire_y, 100)

B = []
for y in y_space:
    # they two wires constructively interfere
    B.append(biot_savart(y - wire_y) - biot_savart(y + wire_y))

B = np.array(B)
plt.plot(y_space, B*1e6)
plt.show()


def calculate_B_at(y):
    return np.interp(y, y_space, B)


R = 0.17/2  # radius of coil in meters
R2 = R**2


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
    if not is_inside_circle(x1, y1):  # completely outside circle
        return 0
    if is_inside_circle(x2, y2):  # completely inside circle
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


stepsize = R/32
grid = np.arange(-R+stepsize/2, R, stepsize)
from itertools import product
flux = 0
for x, y in product(grid, grid):

    """ now to calculate how much of the rectangle is inside the coil"""
    A = rect_area_inside_circle((x, y), stepsize)

    B_here = calculate_B_at(y)
    print(B_here)
    quit()
    flux += B_here*A

print("flux with integral vs flux with center value times area: ", flux, np.pi*R2*calculate_B_at(0))
I2flux = flux/I_rms
print("constant to get flux from current:", I2flux)
"""
the flux_rms through the coil at I_rms is thus known, the emf in the coil will be N*dflux/dt.
however, I should be considered as the sum of frequency components, measured via fft
given by arrays I_fs anf I_as for frequency and amplitude of each component.

I(t) is then sum of I_as[i]*sin(2pi*I_fs[i]*t) and the derivative is then
emf(t) = sum of N*I_as[i]*2pi*I_fs[i]*cos(2pi*I_fs[i]*t)

"""
