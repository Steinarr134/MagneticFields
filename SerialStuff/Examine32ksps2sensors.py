import os

from matplotlib import pyplot as plt
import numpy as np
import scipy.fftpack
from scipy import signal
plt.style.use('dark_background')

# first of what I measured under Suðurnesjalína
filename = "Data\\magnetic_raw_2024_01_10___11_37_37.txt"

# newest file:
# files = os.listdir("Data")
# filename = "Data\\" + files[-1]

print(filename)

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


# function to get sag from
def func(x):
    # return 2.09385992/np.sinh(x - 1.0043624) +30.68953421
    # return 4.66487133e+03/np.cosh(x -8.58502502e-01) -4.65016916e+03
    return 4.66487133e+03/np.cosh(x -8.58502502e-01) -4.65016916e+03



w = 1/32000
for dump_nr in range(1):
    data = []
    A = []
    B = []
    for i, line in enumerate(lines):
        # using float format'
        # tabs = line.split(b'\t')
        if b'was' in line:
            B = np.array(data)
            data = []
            tspace = np.arange(0, A.shape[0]*w, w)

            B -= np.average(B)
            A -= np.average(A)

            # start = find_first_period_start(A)
            # end = find_last_period_end(A)
            # print(start, end)
            # A = np.array(A[start:end])
            # A = A - np.average(A)
            # T = T[start:end]

            plt.figure("Raw waveform")
            plt.plot(tspace, A)
            plt.plot(tspace, B)
            plt.title(r"32k sps")
            plt.xlabel("time [ms]")
            plt.ylabel(r"emf [$\mu V$]")

            plt.figure("Ratio")
            plt.plot(tspace, A/B)
            # plt.show()
            # continue
            #
            # plt.subplot(1, 2, 2)
            N = len(A)
            # yf = scipy.fftpack.fft(scipy.signal.windows.flattop(N)*np.array(A))
            yf = scipy.fftpack.fft(np.array(A))
            xf = np.linspace(0.0, 32000/2, N//2)

            # fig, ax = plt.subplots()
            ybla =  2.0/N * np.abs(yf[:N//2])
            # import pickle
            # with open(r"C:\Users\Notandi\PycharmProjects\MagneticFields\LabTestFrequencyResponseOfCoil\fft_coil_emf.pkl", 'wb+') as pf:
            #     pickle.dump((xf, ybla), pf)
            # with open(r"C:\Users\Notandi\PycharmProjects\MagneticFields\LabTestFrequencyResponseOfCoil\coil_emf.pkl", 'wb+') as pf:
            #     pickle.dump((A, T), pf)
            plt.figure("fft")
            plt.plot(xf, ybla)
            plt.xlabel("Frequency [Hz]")
            plt.ylabel(r"emf [$\mu V$]")
            plt.title("FFT")
            # plt.show()
            # quit()

            plt.figure("butterworth")
            # # apply a butterworth lowpass filter of different degrees
            for j in range(1, 5):
                sos = signal.butter(j, 150, 'lp', fs=32000, output='sos')
                filteredA = signal.sosfilt(sos, A)
                filteredB = signal.sosfilt(sos, B)
                plt.plot(tspace, filteredA, label=f"{j}-order | A")
                plt.plot(tspace, filteredB, label=f"{j}-order | B")
                print("ratio: ", np.max(filteredA)/np.max(filteredB))
            plt.legend()

            a1 = np.max(filteredA)
            a2 = np.max(filteredB)

            print(func(a2/a1))

            plt.show()

            A_use = A[:4400]
            print(np.mean(A_use))
            # plt.plot(A_use)
            # plt.show()
            # rms.append(np.sqrt(np.mean(np.square(A_use - np.mean(A_use)))))

            A = []
            T = []
            continue

        elif b'A' in line:
            pass
        elif b'B' in line:
            A = -np.array(data)
            data = []
        elif line.strip():
            data.append(float(line.strip()))
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