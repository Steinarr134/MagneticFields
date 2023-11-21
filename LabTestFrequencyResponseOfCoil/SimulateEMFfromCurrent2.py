import numpy as np
import scipy.signal
from matplotlib import pyplot as plt
plt.style.use("dark_background")
"""Based of the simulation done in I_Coil_in_Cable_running_lightbulb"""

# I2flux = -4.7420568253394683e-08
I2flux = -4.8572698433768217e-08
Resistor = 20.1

""" import the measured current """
with open("MA1.CSV", 'r') as f:
    lines = f.readlines()

t_space = []
current = []
for line in lines[1:]:
    stuff = line.split(',')
    t_space.append(float(stuff[0]))
    current.append(float(stuff[1].strip())/Resistor)

dI = np.gradient(current, t_space)
emf_sim = 126 * I2flux * dI * 1e6  # 126 windings, 1e6 to express in micro volts
emf_sim = np.roll(emf_sim, -270)
# plt.plot(t_space, np.gradient(current, t_space))
# plt.plot(t_space, emf)


N = len(emf_sim)
sampling_f = 1/np.average(np.diff(t_space))
yf = scipy.fftpack.fft(emf_sim)
xf = scipy.fftpack.fftfreq(N, 1/sampling_f)
xbla = np.linspace(0.0, sampling_f/2, N//2)

ybla = 2.0/N * np.abs(yf[:N//2])

plt.plot(xbla, ybla)
import pickle
""" load measured emf fft from file """
with open("fft_coil_emf.pkl", 'rb') as pf:
    (coil_x_fft, coil_y_fft) = pickle.load(pf)

plt.plot(coil_x_fft*0.995, coil_y_fft)
plt.legend(["Simulated emf from current", "Measured emf"])
plt.xlabel("Frequency [Hz]")
plt.ylabel(r"emf [$\mu V]$")


""" calculate the An constant for different frequency values"""
plt.figure()
with open("coil_emf.pkl", 'rb') as pf:
    (coil_emf, coil_t) = pickle.load(pf)

coil_t = np.array(coil_t)*1.00595  # fix oscillator error
t_space = np.array(t_space)
coil_t = np.array(coil_t)
t_space -= t_space[0]
coil_t -= coil_t[0]
emf_sim = -emf_sim
plt.plot(t_space, emf_sim)
plt.plot(coil_t, coil_emf)
plt.figure()
# quit()
Bn_coil = []
Bn_sim = []
odd_frequencies = [50, 150, 250, 350, 450, 550, 650, 750, 850, 950, 1050, 1150]
for f in odd_frequencies:
    Bn_coil.append(np.trapz(coil_emf*np.sin(2*np.pi*f*coil_t), coil_t)*2/(coil_t[-1]-coil_t[0]))
    Bn_sim.append(np.trapz(emf_sim * np.sin(2 * np.pi * f * t_space), t_space) * 2 / (t_space[-1] - t_space[0]))
    # plt.plot(coil_t, coil_emf*np.sin(2*np.pi*f*coil_t))
    # plt.plot(t_space, emf*np.sin(2*np.pi*f*t_space))
    # print(f, Bn_sim, Bn_coil)
    # plt.show()
plt.plot(odd_frequencies, np.abs(Bn_sim))
plt.plot(odd_frequencies, np.abs(Bn_coil))
plt.legend(["Simulated emf from current", "Measured emf"])
plt.title("Fourier series component of the frequencies")
plt.figure()

"""apply a butterworth filter on both signals and compare them"""

coil_fs = 1/np.average(np.diff(coil_t))
sim_fs = 1/np.average(np.diff(t_space))
factor = sim_fs/coil_fs

width = 3
emf_sim_filtered = [emf_sim[0]]
for i in range(1, len(coil_t)):
    c = round(i*factor)
    emf_sim_filtered.append(np.average(emf_sim[c-width-1:c+width]))

emf_sim_filtered = scipy.ndimage.uniform_filter1d(emf_sim_filtered, 3)

plt.plot(coil_t, emf_sim_filtered)
plt.plot(coil_t, coil_emf)
plt.show()
