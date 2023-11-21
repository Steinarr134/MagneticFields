import pickle
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
plt.style.use("dark_background")

with open("fft_coil_emf.pkl", 'rb') as f:
    (coil_x, coil_y) = pickle.load(f)
with open("fft_current.pkl", 'rb') as f:
    (current_x, current_y) = pickle.load(f)


plt.plot(coil_x, coil_y/np.max(coil_y), 'r')
plt.plot(current_x, current_y/np.max(current_y), 'g')

def norm(x):
    return x/np.max(x)

coil_y_max = np.max(coil_y)
coil_peaks_i = find_peaks(norm(coil_y), threshold=0.001)[0]
plt.plot(coil_x[coil_peaks_i], coil_y[coil_peaks_i]/coil_y_max, '*r')


current_y_max = np.max(current_y)
current_peaks_i = find_peaks(norm(current_y), threshold=0.0001)[0]
plt.plot(current_x[current_peaks_i], current_y[current_peaks_i]/current_y_max, '*g')

frequencies_to_compare = [50, 150, 250, 450, 550, 650, 750, 850, 950, 1150 ]
plt.figure()
print(current_x[current_peaks_i])
print(coil_x[coil_peaks_i])

coil_fs = []
current_fs = []
for f in frequencies_to_compare:
    coil_i = coil_peaks_i[np.abs(coil_x[coil_peaks_i] - f).argmin()]
    current_i = current_peaks_i[np.abs(current_x[current_peaks_i] - f).argmin()]
    print(coil_x[coil_i], current_x[current_i])
    coil_fs.append(coil_y[coil_i])
    current_fs.append(current_y[current_i])

coil_fs = np.array(coil_fs)
current_fs = np.array(current_fs)
plt.plot(frequencies_to_compare, current_fs/coil_fs)
plt.show()