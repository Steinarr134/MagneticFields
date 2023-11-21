import numpy as np
from matplotlib import pyplot as plt

plt.style.use('dark_background')

with open("MA1.CSV", 'r') as f:
    lines = f.readlines()

t = []
m = []
for line in lines[1:]:
    stuff = line.split(',')
    t.append(float(stuff[0]))
    m.append(float(stuff[1].strip()))

import scipy
N = len(m)
sampling_f = 1/np.average(np.diff(t))
print(sampling_f)
yf = scipy.fftpack.fft(m)
xf = scipy.fftpack.fftfreq(N, 1/sampling_f)
xbla = np.linspace(0.0, sampling_f/2, N//2)

ybla =  2.0/N * np.abs(yf[:N//2])
#
import pickle
with open("fft_current.pkl", 'wb+') as f:
    pickle.dump((xbla, ybla), f)
plt.plot(xbla, ybla)
plt.show()

