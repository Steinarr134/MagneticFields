import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt


Y = [15.0, 14.5, 14.0, 13.5, 13.0, 12.5, 12.0, 11.5, 11.0, 10.5, 10.0, 9.5, 9.0, 8.5, 8.0, 7.5]
X = [17.77776511775287, 18.6819014767805, 19.663030940132533, 20.73055702250977, 21.8954217656369, 23.170428265480496, 24.57064715948238, 26.113933601724696, 27.821591255581417, 29.719234277392086, 31.837919453110146, 34.21565226660584, 36.89941874458452, 39.94796951163061, 43.45440950057679, 47.5309755679478]
X = np.array(X)
Y = np.array(Y)

def func(x, a, b, c):
    # return (d*x+c)**a + b
    return a*(x+b)**(-1) + c

popt, pcov = curve_fit(func, X, Y)

ys = func(X, *popt)
print(popt)
print("residuals: ", np.sum(np.square(ys - Y)))

xs = np.arange(np.min(X), np.max(X), 0.1)
plt.plot(X, Y)
plt.plot(xs, func(xs, *popt))
plt.show()