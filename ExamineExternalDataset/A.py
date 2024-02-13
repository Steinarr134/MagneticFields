import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy import fft
plt.style.use("dark_background")

df = pd.read_csv("sigId-893.csv")
f = 60

"""
Time,Current.Ia,Current.Ib,Current.Ic,Current.In,Voltage.Va,Voltage.Vb,Voltage.Vc
"""


def get_sin_params0(tspace, samples):
    tspace = tspace[::20]
    samples = samples[::20]
    def sin(t, A, p, w):
        return A*np.sin(w*t + p)

    p0 = [(np.max(samples) - np.min(samples))/2, -0.5, np.pi*2*f]
    eh, ehannad = curve_fit(sin, tspace, samples, p0)
    # print(eh, ehannad)
    return eh[0], eh[1]


from pylab import matrix, cos, sin, lstsq, norm, pi
from math import atan2
def get_sin_params1(tList,yList):
    tList = tList[::20]
    yList = yList[::20]
    '''
        freq in Hz
        tList in sec
    returns
        phase in degrees
    '''
    freq = 60
    b = matrix(yList).T
    rows = [ [sin(freq*2*pi*t), cos(freq*2*pi*t), 1] for t in tList]
    A = matrix(rows)
    (w,residuals,rank,sing_vals) = lstsq(A,b)
    phase = atan2(w[1,0],w[0,0])*180/pi
    amplitude = norm([w[0,0],w[1,0]],2)
    bias = w[2,0]
    return (phase,amplitude)


T = 65e-6

def get_sin_params2(x, y):
    N = x.size
    yf = np.fft.fft(y) # to normalize use norm='ortho' as an additional argument

    # Where is a 200 Hz frequency in the results?
    freq = np.fft.fftfreq(x.size, d=T)
    index, = np.where(np.isclose(freq, 60, atol=1/(T*N)))

    # Get magnitude and phase
    magnitude = np.abs(yf[index[0]])
    phase = np.angle(yf[index[0]])
    print(magnitude, phase)
    plt.plot(freq[0:N//2], 2/N*np.abs(yf[0:N//2]), label='amplitude spectrum')   # in a conventional form
    plt.plot(freq[0:N//2], np.angle(yf[0:N//2]), label='phase spectrum')
    plt.legend()
    plt.grid()
    plt.show()
    return magnitude, phase

def get_sin_params(t, y):
    # This one is stupid and only works when the signal is very close to being a sine wave
    magnitude = (np.max(y) - np.min(y))/2

    # fos phase just find crossover point
    ready = False
    crossover = 0
    for i in range(len(y)):
        if y[i] < 0:
            ready = True
        if ready and y[i]>0 and y[i+1]>0:
            crossover = i
            break
    phase = -t[crossover]*2*np.pi*60
    # plt.plot(t, y)
    # plt.plot(t[crossover], y[crossover], 'r+')
    # plt.plot(t, magnitude*np.sin(2*np.pi*60*))
    return magnitude, phase

def get_unbalance(Ia, Ib, Ic, plot=False):

    Aa, Pa = get_sin_params(tspace, Ia)
    Ab, Pb = get_sin_params(tspace, Ib)
    Ac, Pc = get_sin_params(tspace, Ic)

    if Aa < 0:
        Pa += np.pi
        Aa *= -1
    if Ab < 0:
        Pb += np.pi
        Ab *= -1
    if Ac < 0:
        Pc += np.pi
        Ac *= -1


    Pb = (Pb -Pa)%(np.pi*2)
    Pc = (Pc -Pa)%(np.pi*2)
    Pa_org = Pa
    Pa = 0
    print(f"{Aa=:.2f}, Pa={np.rad2deg(Pa):.2f}")
    print(f"{Ab=:.2f}, Pb={np.rad2deg(Pb):.2f}")
    print(f"{Ac=:.2f}, Pc={np.rad2deg(Pc):.2f}")




    unbalanced = np.array([[
        Aa*np.exp(Pa*1j),
        Ab*np.exp(Pb*1j),
        Ac*np.exp(Pc*1j),
    ]]).T
    a = np.exp(np.pi*2j/3)
    T = np.array([
        [1, 1, 1],
        [1, a**2, a],
        [1, a, a**2]
    ])

    symmetrics = 1/3*np.matmul(T, unbalanced)

    if not ((abs(abs(Pa - Pb)%(np.pi*2)-np.pi*2/3) < np.deg2rad(10))
            or (abs(abs(Pa - Pb)%(np.pi*2)-np.pi*4/3) < np.deg2rad(10))):
        plot = True

    print(symmetrics)
    A0 = np.real(symmetrics[0])
    An = np.real(symmetrics[1])
    Ap = np.real(symmetrics[2])
    P0 = np.imag(symmetrics[0])
    Pn = np.imag(symmetrics[1])
    Pp = np.imag(symmetrics[2])

    if Ap < 0:
        Pp += np.pi
        Ap *= -1
    if A0 < 0:
        P0 += np.pi
        A0 *= -1
    if An < 0:
        Pn += np.pi
        An *= -1

    P0 = (P0-Pp)%(np.pi*2)
    Pn = (Pn-Pp)%(np.pi*2)
    Pp -= 0

    # if Pn < 3:
    #     plot = True

    if plot:
        plt.plot(tspace, Ia, 'b')
        plt.plot(tspace, Ib, 'g')
        plt.plot(tspace, Ic, 'r')
        plt.plot(tspace, Aa*np.sin(f*np.pi*2*tspace + Pa_org), 'b')
        plt.plot(tspace, Ab*np.sin(f*np.pi*2*tspace + Pb + Pa_org), 'g')
        plt.plot(tspace, Ac*np.sin(f*np.pi*2*tspace + Pc + Pa_org), 'r')



        # plt.figure()
        # colors = 'bry'
        # def myplot(p1, p2, c):
        #     plt.plot(np.array([p1[0], p2[0]]), np.array([p1[1], p2[1]]), c)
        # I0 = 400
        #
        # As = [Aa, Ab, Ac]
        # Ps = [Pa, Pb, Pc]
        # for p in range(3):
        #     myplot([-I0*1.5, 0], [-I0*1.5 + As[p]*np.cos(Ps[p]),  As[p]*np.sin(Ps[p])], colors[p])
        #     # myplot([-I0*1.5, 0], [-I0*1.5 + np.real(phasors[p, 0]), np.imag(phasors[p, 0])], colors[p])
        #     # myplot([0, 0], [negampl*np.cos(negphase - p*np.pi*2/3), negampl*np.sin(negphase - p*np.pi*2/3)], colors[p])
        #     # myplot([I0*0.5 + p, p], [I0*0.5 + p + zerampl*np.sin(zerphase), p+zerampl*np.cos(zerphase)], colors[p])
        plt.show()


    return (Ap, Pp), (An, Pn), (A0, P0)



N = 500
Ans = []
Pns = []
A0s = []
P0s = []
bigtime = []
for i in range(4000):
    start = i*N//10
    stop = i*N//10 + N
    tspace = np.array(df["Time"][start:stop]*1e-6)
    Ia = np.array(df["Current.Ia"][start:stop])
    Ib = np.array(df["Current.Ib"][start:stop])
    Ic = np.array(df["Current.Ic"][start:stop])

    (Ap, Pp), (An, Pn), (A0, P0) = get_unbalance(Ia, Ib, Ic)#, plot=tspace[0]>0.115)
    print(An, A0)
    Ans.append(An)
    Pns.append(Pn)
    A0s.append(A0)
    P0s.append(P0)
    bigtime.append(tspace[0])

plt.subplot(1, 2, 1)
plt.plot(bigtime, A0s, label="A0")
plt.plot(bigtime, Ans, label="An")
plt.legend()
plt.subplot(1, 2, 2)
plt.plot(bigtime, P0s, label="P0")
plt.plot(bigtime, Pns, label="Pn")
plt.legend()

plt.show()