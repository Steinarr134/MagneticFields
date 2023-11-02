import numpy as np
from matplotlib import pyplot as plt

mu = 1.25e-6
I = 440
Phase_offset_x = 10
Y0 = 15
freq = 50
Phase1_x = -Phase_offset_x
Phase2_x = 0
Phase3_x = Phase_offset_x
Ps_x = [Phase1_x, Phase2_x, Phase3_x]
O = np.zeros(2)

maxBs = []
ground_clearances = []
sags = np.arange(0, 5, 0.1)
for sag in sags[0:1]:
    Bs = []
    for t in np.arange(0, 1/freq, 1/(freq*50)):
        I_1 = I*np.sin(t*freq*2*np.pi - 120*np.pi*2/360)
        I_2 = I*np.sin(t*freq*2*np.pi)
        I_3 = I*np.sin(t*freq*2*np.pi + 120*np.pi*2/360)


        B = np.zeros(2)  # Magnetic field at measuring point (origin)

        for P_x, I_n in zip(Ps_x, [I_1, I_2, I_3]):
            P = np.array([P_x, Y0 - sag])
            OP = P-O
            OP_length = np.sqrt(np.sum(np.power(OP, 2)))
            OP_unit = OP/OP_length
            OP_unit_perp = np.array([OP_unit[1], -OP_unit[0]])
            B_length = mu*I_n/(2*np.pi*OP_length)
            # print(OP_unit_perp*B_length)
            B += OP_unit_perp*B_length
        Bs.append(np.sqrt(np.sum(np.power(B, 2))))
        print(B, I_1, I_2, I_3)
    maxB = np.max(Bs)
    print(maxB)
    maxBs.append(maxB*1e6)
    ground_clearances.append(Y0-sag)

plt.plot(ground_clearances, maxBs)
plt.title("Ground Clearance vs Magnetic strength")
plt.xlabel("Ground clearance [m]")
plt.ylabel("Megnetic field strength [micro T]")
plt.grid(1)
plt.show()

