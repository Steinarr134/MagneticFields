import os

from matplotlib import pyplot as plt
import numpy as np
import scipy.fftpack
from scipy import signal
plt.style.use('dark_background')
# 4k sps of just minimal noise
# filename = "Data\\magnetic_raw_2023_10_23___13_08_21.txt"

# recorded under otl at 32 ksps
# filename = "Data\\magnetic_raw_2023_10_30___21_05_13.txt"
#
# recorded on my living room floor
# filename = "Data\\magnetic_raw_2023_10_31___14_13_21.txt"

# recorded out in an empty lot approx. 1 km from otl
# filename = "Data\\magnetic_raw_2023_10_31___16_12_28.txt"

# newest file:
files = os.listdir("Data")
for filename in files[-4:]:
    filename = "Data\\" + filename

    print(filename)

    with open(filename, 'rb') as f:
        lines = f.readlines()

    A = []
    B = []
    T = []
    t = 0

    for i, line in enumerate(lines):
        # using float format'
        # tabs = line.split(b'\t')
        if b'was' in line:
            # plt.figure(filename)
            plt.subplot(1, 2, 1)
            plt.plot(T, A)
            # plt.show()
            #
            plt.subplot(1, 2, 2)
            N = len(A)
            yf = scipy.fftpack.fft(A)
            xf = np.linspace(0.0, 32000/2, N//2)

            # fig, ax = plt.subplots()
            ybla =  2.0/N * np.abs(yf[:N//2])
            plt.plot(xf, ybla)
            # plt.show()
            break

        A.append(float(line.strip()))
        # B.append(float(tabs[1]))
        T.append(i/32000.)
plt.show()