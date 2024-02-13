import numpy as np
from matplotlib import pyplot as plt
plt.style.use("dark_background")


Zresult = np.array([[ 4.98090117e-07, -5.95282992e-04,  1.10246665e-01,
                       2.00928394e+02, -1.82573096e+05],
                     [-1.92970113e-06,  2.31998914e-03, -4.68947765e-01,
                      -7.21608870e+02,  6.54813166e+05],
                     [ 2.78805049e-06, -3.36815725e-03,  7.28114120e-01,
                       9.72605511e+02, -8.81227134e+05],
                     [-1.78192012e-06,  2.16120865e-03, -4.92690343e-01,
                      -5.83048623e+02,  5.27340255e+05],
                     [ 4.25348354e-07, -5.17567616e-04,  1.23169118e-01,
                       1.31157723e+02, -1.18377700e+05]])

from M1results import Z as Zresult

def get_F_m(D, Y0, d_P):
    stuff = [1, D, Y0, d_P]

    def _get_pow(depth, p):
        return stuff[depth]**p

    def rec_inner(Z, depth=0):
        # it's recursive
        # print(depth, Z.ndim)
        ret = []
        if Z.ndim == 0:
            return Z
        else:
            for i, item in enumerate(Z):
                p = Z.shape[0] - i - 1
                inner = rec_inner(item, depth=depth+1)
                ret.append(inner*_get_pow(depth, p))
                # if len(inner) == 0:
                #     continue
                # elif len(inner) == 1:
                #     ret.append(inner[0] + _get_pow(depth, p))
                # else:
                #     if inner[-1] == "":
                #         inner[-1] = "1"
                #     inner = " + ".join(inner)
                #     inner = inner.replace("+ -", "- ")
                #     ret.append("(" + inner + ")" + _get_pow(depth, p))
            if depth > 0:
                return np.sum(ret)
            else:
                return ret
    return rec_inner(Zresult)


if __name__ == '__main__':
    sagspace = np.linspace(2, 8)
    D_bounds = (200, 400)
    d_P_bounds = (5, 8)
    Y0_bounds = (10, 15)

    def rand(bounds):
        return bounds[0] + np.random.random()*(bounds[1]-bounds[0])

    from M1_hyperparameter_research import F

    while True:
        D = rand(D_bounds)
        Y0 = rand(Y0_bounds)
        d_P = rand(d_P_bounds)

        Y0 = 12
        d_P = 7


        ratios = F(sagspace, D, Y0, d_P)
        zs = get_F_m(D, Y0, d_P)
        F_m = np.poly1d(zs)
        print(F_m)
        print(np.poly1d(np.polyfit(ratios, sagspace, 4)))
        plt.plot(ratios, F_m(ratios))
        plt.plot(ratios, sagspace)
        plt.figure()
        plt.plot(ratios, F_m(ratios) - sagspace)
        plt.show()
        # break


