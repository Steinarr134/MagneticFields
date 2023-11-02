from matplotlib import pyplot as plt
import numpy as np
plt.style.use('dark_background')
# 4k sps of just minimal noise
# filename = "Data\\magnetic_raw_2023_10_23___13_08_21.txt"

# recorded under otl at 1000 sps
filename = "Data\\magnetic_raw_2023_10_29___15_32_22.txt"
with open(filename, 'rb') as f:
    lines = f.readlines()

last_i = 0
missing_count = 0
A = []
B = []
T = []
t = 0

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
# for i, line in enumerate(lines[start + s*3900:(start+(s+1)*3900)]):
for i, line in enumerate(lines[1000:]):
    # using float format'

    # tabs = line.split(b'\t')
    # A.append(float(tabs[0]))
    # B.append(float(tabs[1]))
    # T.append(i*0.00097)

    # for hexaformat
    header = line[0]
    pps = header & 0x40
    header_i = header & 0x3f
    t += (header_i - last_i)%0x40
    T.append(t)
    A.append(uV_A(line[1:7]))
    B.append(uV_B(line[7:13]))
    if (last_i+1)%0x40 != header_i:
        print("missing line")
        missing_count += 1
        if last_i+2 != header_i:
            print("double trouble", last_i, header_i)
        print(f"{i=}, last was {last_i}, expecting {(last_i+1)%0x40}, got {header_i}")
    last_i = header_i%0x40

print(f"{missing_count=}, total={len(lines)}")

plt.plot(T, B)
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
# plt.plot(t_rms, sag)
plt.show()