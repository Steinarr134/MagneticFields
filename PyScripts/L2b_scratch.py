from HyParI import *
from scipy.optimize import least_squares

xx = np.load("cache/xx.npy")
xy = np.load("cache/xy.npy")
correct = np.load("cache/correct.npy")


def fun2(x):
    ans = correct - (x[0]*xx + x[1]*xy + x[2]*xx**2 + x[3]*xy**2 + x[4]*xx*xy)
    return ans.flatten()
def fun3(x):
    ans = correct - (1 + x[0]*xx + x[1]*xx**2)*(1 + x[2]*xy + x[3]*xy**2)
    return ans.flatten()
def fun1(x):
    ans = correct - (x[0]*xx + x[1]*xy)
    return ans.flatten()


# results = least_squares(fun2, np.random.random(5))
results = least_squares(fun1, np.random.random(2))
x = results["x"]
print(x, np.max(np.abs(fun1(x))))
