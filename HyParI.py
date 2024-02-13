"""
These are two functions (developed in EqDiscovery/A4.py)
That find the relations of hyperparameters, assuming the function is modelled
as a polynomial

"""
import itertools
import numpy as np
from inspect import signature
from scipy.optimize import curve_fit
from matplotlib import  pyplot as plt


def get_ndimquation(Z, letters="xabcdefg", cutoff=1e-6):
    def _get_pow(d, p):
        if p == 0:
            return ""
        elif p == 1:
            return letters[d]
        else:
            return letters[d] + "^" + str(p)

    def rec_inner(Z, depth=0):
        # it's recursive
        ret = []
        if Z.ndim == 0:
            if abs(Z-1) < cutoff:
                return [""]
            elif abs(Z) < cutoff:
                return []
            else:
                return [f"{Z:.3f}"]
        else:
            for i, item in enumerate(Z):
                p = Z.shape[0] - i - 1
                inner = rec_inner(item, depth=depth+1)
                if len(inner) == 0:
                    continue
                elif len(inner) == 1:
                    ret.append(inner[0] + _get_pow(depth, p))
                else:
                    if inner[-1] == "":
                        inner[-1] = "1"
                    inner = " + ".join(inner)
                    inner = inner.replace("+ -", "- ")
                    ret.append("(" + inner + ")" + _get_pow(depth, p))
            if depth > 0:
                return ret
            else:
                ret = " + ".join(ret)
                ret = ret.replace("+ -", "- ")
                return ret

    return rec_inner(Z)


def n_hyperparameters(F, F_m_order, hyper_order, bounds, xspace, n=10):
    # bounds define the number of hyperparamters
    n_hyper = len(bounds)

    boundspaces = [np.linspace(*bounds[i], n) for i in range(n_hyper)]

    # This is supposed to work as a recursive function,
    # however the outermost run is slightly different
    # So first I'll define this inner recursive function

    def rec_inner(Z, depth=0):
        # At this point Z holds the values of a Z constant as a function
        # of m hyperparameters, thus Z is an m-dimensional array, each
        # dimension of length n

        # First the exit condition, if we have reached the final depth
        if depth == n_hyper - 1:
            zs = np.polyfit(boundspaces[depth], Z, hyper_order)
            from matplotlib import pyplot as plt
            plt.plot(boundspaces[depth], Z)
            plt.plot(boundspaces[depth], np.poly1d(zs)(boundspaces[depth]))
            plt.show()
            return zs
        else:
            # Z is an m-dimensional array holding the value of Z as a function of
            # the last m boundspaces.
            # Here we assume Z is a polynomial with parameters Z_a(b, c, ....):
            # Z(a, b, c, ...) = Z_a2(b, c, ...)*a^2 + Z_a1(b, c, ...)*a + Z_a0(b, c, ...)
            Z_a = np.zeros([n]*(n_hyper - depth-1) + [hyper_order + 1])
            for idx in itertools.product(*[list(range(n)) for _ in range(n_hyper-depth-1)]):
                Z_a[idx] = np.polyfit(boundspaces[depth], Z[*[slice(None)] + list(idx)], hyper_order)

                # from matplotlib import pyplot as plt
                # poly = np.poly1d(Z_a[idx])
                # plt.plot(boundspaces[depth], Z[*[slice(None)] + list(idx)])
                # plt.plot(boundspaces[depth], poly(boundspaces[depth]))
                # plt.show()

            inret = []
            for i_h in range(hyper_order + 1):
                inret.append(rec_inner(Z_a[*[slice(None)]*(n_hyper-depth-1) + [i_h]], depth=depth+1))
            return inret


    # first work out all the values of F as a function of x and all hyperparameters
    Z_all = np.zeros([n]*n_hyper + [F_m_order + 1])
    idxs = itertools.product(*[list(range(n)) for _ in range(n_hyper)])
    N = n**(n_hyper)
    last_percent = 0
    for i, idx in enumerate(idxs):
        hypers = [boundspaces[i][idx[i]] for i in range(n_hyper)]
        y = F(xspace, *hypers)
        zs = np.polyfit(xspace, y, F_m_order)
        Z_all[idx] = zs
        if int(95*i/N) > last_percent:
            last_percent = int(95*i/N)
            print(f"{last_percent}%")


        # from matplotlib import pyplot as plt
        # plt.plot(xspace, y)
        # plt.plot(xspace, np.poly1d(zs)(xspace))
        # plt.show()

    # in the case of two hyperparameters, Z_all is a nxnx3 matrix

    # the return value is a multidimensional array that we'll build with lists
    ret = []
    # loop through the each first order parameter
    # F(x, a, b) = Z_2(a,b)*x^2 + Z_1(a,b)x + Z_0(a,b)
    for i_m in range(F_m_order + 1):
        # extract the Z_n
        Z_i = Z_all[*[slice(None)]*n_hyper + [i_m]]
        ret.append(rec_inner(Z_i))

    ret = np.array(ret)
    # ret[np.abs(ret) < 1e-6] = 0
    return ret

if __name__ == '__main__':
    count = 0
    def F(x, a, b):
        # print(a, b)
        global count
        count += 1
        return (a**2 + a*b + b**2 + 2)*x**2 + (b*a**2 + 3*b**2*a - 1)*x + a**2*b**2 + b*a + a + b + 5

    order = 2
    Z = n_hyperparameters(F, order, order, ((3, 10), (1, 5)), np.linspace(5, 15, 20))
    print(repr(Z))
    print(get_ndimquation(Z))
    print(count)

    # Test with 3 parameters
    count = 0
    def F(x, a, b, c):
        # print(a, b)
        global count
        count += 1
        return (a**2*c**2 + a*b*c + b**2*c + 3*c**2 + 2)*x**2 \
            + (b*a**2*c + 3*b**2*a*(c**2 + 1) + c*4 - 1)*x \
            + c*a**2*b**2 + b*a + a + b*c**2 + 5

    order = 2
    Z = n_hyperparameters(F, order, order, ((3, 10), (1, 5), (-5, 1)), np.linspace(5, 15, 20))
    print(repr(Z))
    print(get_ndimquation(Z))
    print(count)



