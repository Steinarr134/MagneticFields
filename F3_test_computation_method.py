import numpy as np
import pickle
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

with open("D250_Yn20_z.pkl", 'rb') as f:
    Is, sags, Y0, D, Bs1, Bs2 = pickle.load(f)

Bs1 = Bs1*1e6
Bs2 = Bs2*1e6

I_slope = 0.8801664872840699

sag_zs = [-6397.04420309,  15075.57886288, -11784.0619688,  3068.56372453]
sag_poly = np.poly1d(sag_zs)

def func(x):
    # return 2.09385992/np.sinh(x - 1.0043624) +30.68953421
    return 4.66487133e+03/np.cosh(x -8.58502502e-01) -4.65016916e+03

"""
Test out computation methods to get I from B1 and B2
"""
for j, sag in enumerate(sags):
    for i, I in enumerate(Is):
        B1 = Bs1[i, j]
        B2 = Bs2[i, j]

        # # to find current first find the intercept of a line with slope I_slope
        # # B2 = I_slope*B1 + intercept
        # intercept = B2 - I_slope*B1
        #
        # # current can then be found with
        # I_estimate = 4639.698799504301*intercept
        # print(I, I_estimate, abs(I-I_estimate))
        # # it's not very good :(
        # # I think it's because the current slopes had large residuals
        # # the sag lines had negligable residuals so let's try them out

        # to find sag first find slope of line through (B1, B2) and (0, 0)
        slope = B2/B1
        sag_estimate = func(slope)
        print(sag, sag_estimate)
