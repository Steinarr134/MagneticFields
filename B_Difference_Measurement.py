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
d = 1
O1 = np.zeros(2)
O2 = np.array([0, -d])

maxBs1 = []
maxBs2 = []
ground_clearances = []
sags = np.arange(0, 5, 0.1)
for sag in sags:
    Bs1 = []
    Bs2 = []
    for t in np.arange(0, 1/freq, 1/(freq*50)):
        I_1 = I*np.sin(t*freq*2*np.pi - 120*np.pi*2/360)
        I_2 = I*np.sin(t*freq*2*np.pi)
        I_3 = I*np.sin(t*freq*2*np.pi + 120*np.pi*2/360)

        # print(t, I_2)

        B1 = np.zeros(2)  # Magnetic field at measuring point (origin)
        B2 = np.zeros(2)  # magnetic field at measuring point 2 (d below origin)
        for P_x, I_n in zip(Ps_x, [I_1, I_2, I_3]):
            P = np.array([P_x, Y0 - sag])
            OP = P-O1
            OP_length = np.sqrt(np.sum(np.power(OP, 2)))
            OP_unit = OP/OP_length
            OP_unit_perp = np.array([OP_unit[1], -OP_unit[0]])
            B_length = mu*I_n/(2*np.pi*OP_length)
            # print(OP_unit_perp*B_length)
            B1 += OP_unit_perp*B_length

            OP = P-O2
            OP_length = np.sqrt(np.sum(np.power(OP, 2)))
            OP_unit = OP/OP_length
            OP_unit_perp = np.array([OP_unit[1], -OP_unit[0]])
            B_length = mu*I_n/(2*np.pi*OP_length)
            # print(OP_unit_perp*B_length)
            B2 += OP_unit_perp*B_length
        Bs1.append(np.sqrt(np.sum(np.power(B1, 2))))
        Bs2.append(np.sqrt(np.sum(np.power(B2, 2))))
    maxB1 = np.max(Bs1)
    maxBs1.append(maxB1*1e6)
    maxB2 = np.max(Bs2)
    maxBs2.append(maxB2*1e6)
    ground_clearances.append(Y0-sag)


plt.plot(ground_clearances, maxBs1)
plt.plot(ground_clearances, maxBs2)
plt.title("Ground Clearance vs Magnetic strength")
plt.xlabel("Ground clearance [m]")
plt.ylabel("Magnetic field strength [micro T]")
plt.grid(1)
plt.show()

