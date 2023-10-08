import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt


Y = [15.0, 14.5, 14.0, 13.5, 13.0, 12.5, 12.0, 11.5, 11.0, 10.5, 10.0, 9.5, 9.0, 8.5, 8.0, 7.5]
X = [2.1140308065658755, 2.2927928929931847, 2.4913555216901986, 2.712400621819753, 2.9590401418257875, 3.2348971449831923, 3.544205276939524, 3.8919318131231186, 4.283931388295731, 4.7271402686210475, 5.229825116544704, 5.8019063783290745, 6.455385942427054, 7.204923675194322, 8.068631454510323, 9.069192777893173]
X = np.array(X)
Y = np.array(Y)

def func(x, a, b, c):
    # return (d*x+c)**a + b
    return a*(x+b)**(-0.25) + c

popt, pcov = curve_fit(func, X, Y)

ys = func(X, *popt)
print(popt)
print("residuals: ", np.sum(np.square(ys - Y)))

xs = np.arange(np.min(X), np.max(X), 0.1)
plt.plot(X, Y)
plt.plot(xs, func(xs, *popt))
plt.show()