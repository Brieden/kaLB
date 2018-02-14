import matplotlib.pyplot as plt
import numpy as np
import h5py
import os


def system_test(args):
    """

    :param args:
    :return:
    """
    f = h5py.File("./output/system_test.hdf5", 'r')
    velocity_names = list(f['raw data output configuration']['velocity'])
    velocity_value = (f['raw data output configuration']['velocity'][velocity_names[-1]])
    speed = velocity_value[0, -3, :]
    x = range(len(speed))
    coefs, residuals, _, _, _ = np.polyfit(x, speed, 2, full=True)
    if residuals < speed.mean() * 1e-4:
        print("The system test flow in the pipe was passed.")
    else:
        print("The system test flow in the pipe was not passed.")
    if args.verbose:
        plt.plot(speed, label="Values of the simulation.")
        b = np.poly1d(coefs)
        plt.plot(x, b(x), "x", label="parabolic fit as a reference")
        plt.title("Velocity profile at the end of a tube")
        plt.xlabel("Position x: at right angles to the wall")
        plt.ylabel("Speed")
        plt.legend()
        plt.show()
