import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

mu = 1.25e-6
I = 500  # Current in wire
Phase_offset_x = 10  # phase spacing
Y0 = 20  # starting height of wires at tower
freq = 50
Phase1_x = -Phase_offset_x
Phase2_x = 0
Phase3_x = Phase_offset_x
Ps_x = [Phase1_x, Phase2_x, Phase3_x]
O = np.zeros(3)  # Origin point  (where measurement is taken)
D = 250  # distance between towers

# given B = 5uT find ground clearance


