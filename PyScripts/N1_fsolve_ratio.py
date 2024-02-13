"""
Goal is to write a function that uses fsolve to work backwards from B2/B1 to h.

However I got a bit distracted and did this for straight wires, still an interesting result.
The curve is completely opposite to what happens with catenary curves.

"""

import numpy as np
from scipy.optimize import fsolve
from matplotlib import  pyplot as plt

def Zx(h, d_P, e0, e1, th0, th1):
    a = 1 + (3*h**2)/(d_P**2)
    Cx = 1 + e1*np.cos(th1) + e0*a*np.cos(th0)
    Dx = e1*np.sin(th1) + e0*a*np.sin(th0)
    print(Cx, Dx)
    return np.sqrt(np.square(Cx)+np.square(Dx))

def Bx(h, d_P, e0, e1, th0, th1):
    z1 = Zx(h, d_P, e0, e1, th0, th1)
    return (d_P**2 * z1)/(h*(h**2 + d_P**2))

def ratio(h, d_P, e0=0, e1=0, th0=0, th1=0, d=0.5):
    B1 = Bx(h, d_P, e0, e1, th0, th1)
    B2 = Bx(h+d, d_P, e0, e1, th0, th1)
    return B2/B1

def phi(h, d_P, e0, e1, th0, th1):
    a = 1 + (3*h**2)/(d_P**2)
    Cx = 1 + e1*np.cos(th1) + e0*a*np.cos(th0)
    Dx = e1*np.sin(th1) + e0*a*np.sin(th0)
    # print(Cx, Dx)
    return np.arctan(Cx/Dx)

def rand(bounds):
    return bounds[0] + np.random.random()*(bounds[1]-bounds[0])

if __name__ == '__main__':
    d_P = 7
    hspace = np.linspace(5, 10)
    ratios = []
    for h in hspace:
        ratios.append(ratio(h, d_P))

    plt.plot(ratios, hspace)
    plt.show()

    # Secondary wondering, are the X components in phase?
    # The phase angle will be:

    for i in range(10):
        e0 = 0.05
        e1 = 0.02
        th1 = rand((-np.pi, np.pi))
        th0 = rand((-np.pi, np.pi))
        h = rand((5, 10))
        d_P = rand((6, 8))

        args = [h, d_P, e0, e1, th0, th1]

        print(args, phi(*args), phi(h+0.5, *args[1:]))