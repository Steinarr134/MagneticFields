from serial import Serial

last_i = 0
with Serial("COM10", 115200) as serin:
    with Serial("COM50", 115200) as serout:
        for i in range(1000):
            serin.readline()
        for j in range(2000):
            incoming = serin.readline()
            # print(incoming)
            # for b in incoming[:8]:
            #     print(hex(b))
            header = incoming[0]
            first_bit = 0b10000000
            isA = header & 0b01000000
            isB = not isA
            isTimestamp = header & 0b00100000
            pps = header & 0b00010000
            i = header & 0x0F

            print(f"header={bin(header)}, it's from {'A' if isA else 'B'},{i=}, {last_i=}")
            if not first_bit:
                print("\nNo first bit")

            if not ((last_i + 1) % 16 == i):
                print("\nMessage lost")

            last_i = i
            number = 0
            for i, b in enumerate(incoming[1:7]):
                # print(f"{b=}, {b-0x41=} => {hex(b-0x41)[2]}")
                number |= ((b - 0x41) << (i*4))

            if isB:
                serout.write(str(number).encode('ascii') + b'\n')
                # print("\r" + str(number), end="")



"""

⸮GDHLPP	-18634	16758582
"""