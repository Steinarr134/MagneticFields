from serial import Serial
from datetime import datetime
import threading
import queue
import time
import os


def get_filename(now=None):
    if not now:
        now = datetime.now()
    return now.strftime("Data" + os.sep + "magnetic_raw_%Y_%m_%d___%H_%M_%S.txt")


class WritingThread(threading.Thread):
    def __init__(self, q):
        self.q = q
        self.filename = None
        self.keep_going = True
        threading.Thread.__init__(self)
        self.daemon = True
        self.start()

    def run(self) -> None:
        self.filename = get_filename()
        with open(self.filename, 'wb+') as f:
            while self.keep_going:
                line = self.q.get()
                if not line:
                    break
                f.write(line)

    def stop(self):
        self.keep_going = False
        self.q.put(b'')
        self.join()


class ReadingThread(threading.Thread):
    def __init__(self, q, s):
        self.q = q
        self.s = s
        self.i = 0
        self.last_speed_t = time.time() - 0.1
        self.last_speed_i = 0
        self.keep_going = True
        threading.Thread.__init__(self)
        self.daemon = True
        self.start()

    def run(self) -> None:
        while self.keep_going:
            line = self.s.readline()
            if not line:
                break
            self.q.put(line)
            self.i += 1

    def stop(self):
        self.keep_going = False
        self.join()

    def speed(self):
        t = time.time()
        i = int(self.i)
        speed = (i - self.last_speed_i)/(t - self.last_speed_t)
        self.last_speed_i = i
        self.last_speed_t = t
        return speed


class DreamingThread(threading.Thread):
    def __init__(self, q, s):
        self.q = q
        self.i = 0
        self.last_speed_t = time.time() - 0.1
        self.last_speed_i = 0
        self.keep_going = True
        threading.Thread.__init__(self)
        self.daemon = True
        self.start()

    def run(self) -> None:
        while self.keep_going:
            time.sleep(0.0005)
            self.i += 1
            line = str(self.i) + "\n"
            if not line:
                break
            self.q.put(line)

    def stop(self):
        self.keep_going = False
        self.join()

    def speed(self):
        t = time.time()
        i = int(self.i)
        speed = (i - self.last_speed_i)/(t - self.last_speed_t)
        self.last_speed_i = i
        self.last_speed_t = t
        return speed


if __name__ == '__main__':

    SecondsPerFile = 60
    # Queue that stores things to be written
    Q = queue.Queue()

    # start generating stuff (read serial or in this case just dream up numbers)

    with Serial("COM10", 2000000) as serial:
        reader = ReadingThread(Q, serial)

        while True:
            # spin up a new writing thread, it starts a new file
            writer = WritingThread(Q)

            t0 = time.time()
            print("sfdlkj")
            try:

                # this portion should take SecondsPerFile to execute
                # all the while displaying (slowly) the following info:
                # * amount of lines read from arduino per second over the last few seconds
                # * amount of lines presently in the queue
                # * current filename
                # * how long until new file
                # It should look like this:
                #  | (filename) | XX s | XXXX,X l/s | XXX in q |
                s = ""
                for i in range(SecondsPerFile):
                    s = str(f"| {writer.filename} | {SecondsPerFile - int(time.time() - t0)} s "
                            f"| {reader.speed():.1f} msg/s | {Q.qsize()} msg in Q |")
                    print("\r" + s, end="")
                    time.sleep(1)
                print(s)
            except KeyboardInterrupt:
                print("stopping reader")
                reader.stop()
                break
            finally:
                print("\nstopping writer")
                writer.stop()
    print("ok bye")
