import numpy as np
import pickle
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy.optimize import curve_fit

with open("D250_Yn20_z.pkl", 'rb') as f:
    Is, sags, Y0, D, Bs1, Bs2 = pickle.load(f)

Bs1 = Bs1*1e6
Bs2 = Bs2*1e6




"""
plot B2 as a function of B1, grouped by sag or current
note that the lines do not overlap!
suggesting that 
"""
colors = "bgrcmykbgrcmykbgrcmyk"
plt.subplot(1, 2, 1)
for i in range(Bs1.shape[1]):
    plt.plot(Bs1[:, i], Bs2[:, i], color=colors[i])
plt.title("B1 vs B2 grouped by sag")
plt.xlabel("B1 [uT]")
plt.ylabel("B2 [uT]")
plt.subplot(1, 2, 2)
for i in range(Bs1.shape[0]):
    plt.plot(Bs1[i,:], Bs2[i, :], color=colors[i])
plt.title("B1 vs B2 grouped by current")
plt.xlabel("B1 [uT]")
plt.ylabel("B2 [uT]")

# investigate the equations for the lines
slopes = []
ks = []
for i in range(Bs1.shape[1]):
    sag = sags[i]
    ((m, k), residuals, _, _, _) = np.polyfit(Bs1[:, i], Bs2[:, i], 1, full=True)
    # print(stuff)
    print(f"for sag = {sag}, equation is B2 = {m:.5f}*B1 + {k:.2e}, residual = {residuals[0]:.3e}")
    slopes.append(m)
    ks.append(k)
    # quit()

logsags = np.log(sags)
logslopes = np.log(slopes)
plt.figure()
plt.plot(slopes, sags)
plt.title("sag as a function of the slope of B2 as a function of B1 ")
plt.ylabel("sag")
plt.xlabel("slope")
(z, residuals, _, _, _) = np.polyfit(slopes, logsags, 2, full=True)
# print(f"sag = {m}*slope + {k}, residual = {residuals[0]}")
print(z, residuals)
print("\n"*5)
# plt.plot(Is, ks)

def func(x, a, b, c):
    return a/np.cosh(x+ c) + b
    # return -((x+c)**a) + b

popt, pcov = curve_fit(func, slopes, sags, bounds=((-np.inf, -np.inf, -np.min(slopes)), np.inf))
ys = func(np.array(slopes), *popt)
print(np.array(sags) - ys)
print("residual: ", np.sum(np.square(np.array(sags) - ys)))
plt.figure()
plt.title("diff of sag(slopes)")
plt.plot(slopes, sags)
xs = np.arange(np.min(slopes), np.max(slopes), 0.001)

plt.plot(xs, func(xs, *popt))
# print((xs, func(xs, *popt)))
print(popt)
plt.show()
quit()

slopes = []
ks = []
for i in range(Bs1.shape[0]):
    I = Is[i]
    ((m, k), residuals, _, _, _) = np.polyfit(Bs1[i, :], Bs2[i, :], 1, full=True)
    # print(stuff)
    print(f"for I = {I}, equation is B2 = {m:.5f}*B1 + {k:.2e}, residual = {residuals[0]:.3e}")
    slopes.append(m)
    ks.append(k)
    # quit()

print(f"Average slope is {np.average(slopes)}")
plt.figure()
((m, k), residuals, _, _, _) = np.polyfit(ks, Is, 1, full=True)
plt.title("Current as a function of the intercept of B2 as a function of B1 \n"
          + f"current = {m:.5e}*intercept + {k:.5e}, residual = {residuals[0]:.3e}")
print(f"current = {m}*intercept + {k}, residual = {residuals[0]}")
plt.ylabel("current")
plt.xlabel("intercept")
plt.plot(ks, Is)
plt.show()
