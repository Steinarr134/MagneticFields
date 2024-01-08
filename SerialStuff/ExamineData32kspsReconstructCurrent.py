from matplotlib import pyplot as plt
import numpy as np
import scipy.fftpack
from scipy import signal
plt.style.use('dark_background')
# 4k sps of just minimal noise
# filename = "Data\\magnetic_raw_2023_10_23___13_08_21.txt"

# recorded under otl at 32 ksps
filename = "Data\\magnetic_raw_2023_10_30___21_05_13.txt"
#
# recorded on my living room floor
# filename = "Data\\magnetic_raw_2023_10_31___14_13_21.txt"

# recorded out in an empty lot approx. 1 km from otl
# filename = "Data\\magnetic_raw_2023_10_31___16_12_28.txt"


with open(filename, 'rb') as f:
    lines = f.readlines()

last_i = 0
missing_count = 0
A = []
B = []
T = []
t = 0

def find_first_period_start(stuff):
    ready = False
    for i in range(len(stuff)):
        if stuff[i] < 0:
            ready = True
        if ready and stuff[i]>0 and stuff[i+1]>0:
            return i


def find_last_period_end(stuff):
    ready = False
    for i in range(1, len(stuff)):
        if stuff[-i] > 0:
            ready = True
        if ready and stuff[-i] < 0:
            return len(stuff)-i


def twos_comp(val, bits):
    """compute the 2's complement of int value val"""
    if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
        val = val - (1 << bits)        # compute negative value
    return val


def raw_signed(h):
    raw = int(h, base=16)
    return twos_comp(raw, 24)

r2uvA = 1.1064*2.048/(32*8.388608)
def uV_A(h):
    return raw_signed(h)*r2uvA
r2uvB = 1.1054*2.048/(32*8.388608)
def uV_B(h):
    return raw_signed(h)*r2uvB

start = 105000
# for s in range(10):
#     T = []
#     A = []
#     B = []
#     t = 0
rms = []
for i, line in enumerate(lines):
    # using float format'
    # tabs = line.split(b'\t')
    if b'was' in line:
        start = find_first_period_start(A)
        end = find_last_period_end(A)
        print(start, end)
        A = np.array(A[start:end])
        A = A - np.average(A)
        T = T[start:end]
        plt.figure("original")
        plt.plot(T, A)
        # plt.show()
        #
        plt.figure("fft")
        N = len(A)
        yf = scipy.fftpack.fft(A)
        xf = scipy.fftpack.fftfreq(N, 1/32000)


        # fig, ax = plt.subplots()
        y_fftplot =  2.0/N * np.abs(yf[:N//2])
        x_fftplot = np.linspace(0.0, 32000/2, N//2)
        # plt.plot(yf)
        # plt.show()
        # plt.plot(yf)
        fixed_yf = []
        for i in range(len(yf)):
            if -100 < xf[i] < 100:
                fixed_yf.append(yf[i])
                continue
            fixed_yf.append(yf[i]*50/xf[i])
        # plt.plot(xf, fixed_yf)
        plt.plot(x_fftplot, y_fftplot)
        plt.plot(x_fftplot, 2.0/N * np.abs(fixed_yf[:N//2]))
        plt.xlabel("Frequency [Hz]")
        plt.title("FFT")
        plt.legend(["emf in coil", "Reconstructed current"])
        # plt.show()
        plt.figure("org fft")
        plt.plot(xf, yf)
        plt.plot(xf, fixed_yf)
        plt.figure("reconstructed")
        plt.plot(A)
        recon = scipy.fftpack.ifft(fixed_yf)
        plt.plot(recon)
        plt.show()

        # # apply a butterworth lowpass filter of different degrees
        # for j in range(1, 20, 3):
        #     sos = signal.butter(j, 150, 'lp', fs=32000, output='sos')
        #     filtered = signal.sosfilt(sos, A)
        #     plt.plot(filtered)
        # plt.show()

        A_use = A[:4400]
        print(np.mean(A_use))
        # plt.plot(A_use)
        # plt.show()
        rms.append(np.sqrt(np.mean(np.square(A_use - np.mean(A_use)))))

        A = []
        T = []
        continue
    A.append(float(line.strip()))
    # B.append(float(tabs[1]))
    T.append(i/32000.)
plt.plot(rms)
plt.show()
plt.plot(T, A)

A_rms = []
B_rms = []
t_rms = range(1950, len(lines)-3000, 500)
sag = []
def func(x):
    # return 2.09385992/np.sinh(x - 1.0043624) +30.68953421
    return 4.66487133e+03/np.cosh(x -8.58502502e-01) -4.65016916e+03
for i in t_rms:
    A_rms.append(np.sqrt(np.mean(np.square(A[i-1950:i+1950]))))
    B_rms.append(np.sqrt(np.mean(np.square(B[i-1950:i+1950]))))
    sag.append(func(A_rms[-1]/B_rms[-1])*1000)
plt.plot(t_rms, A_rms)
plt.plot(t_rms, B_rms)
plt.plot(t_rms, sag)
plt.show()