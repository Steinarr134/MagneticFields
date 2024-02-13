from L2_surfaceplot_reworked import *
np.random.seed(42)
D = rand(D_bounds)
Y0 = rand(Y0_bounds)
d_P = rand(d_P_bounds)
sag = rand(sagspace)
guesses = []
for h in np.linspace(5, 10):
    xyratio = calculate_ratio(h, Y0, D, d_P, 0, 0, 0, 0)[1]
    print(xyratio*d_P/np.sqrt(3) - h)