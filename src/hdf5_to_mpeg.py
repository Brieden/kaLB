import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
import argparse
import h5py
import os

plt.style.use('bmh')
plt.rcParams['figure.figsize'] = (16.0, 9.0)


def parse_arguments():
    """
    Parse commandline arguments.
    :return: args
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input', required=True, type=str,
        help="Specify path to input hdf5 file."
    )
    return parser.parse_args()


def make_density_pictures(number):
    plt.imshow(np.array(density_values[number]).T, origin="lower")
    plt.title("t = %i" % velocity_names[number])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("./temp_png_to_mp4/density_%i.png" % number)
    plt.cla()


def make_velocity_pictures(number):
    plt.imshow(np.sqrt(np.array(velocity_values[number][0])**2 +
                np.array(velocity_values[number][1])**2).T, origin="lower")
    plt.title("t = %i" % velocity_names[number])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("./temp_png_to_mp4/velocity_%i.png" % number)
    plt.cla()


args = parse_arguments()

f = h5py.File(args.input, 'r')
a_group_key = list(f.keys())

if 'raw data output configuration' in a_group_key:

    if 'density' in list(f['raw data output configuration']):
        density_values = []
        density_names = list(f['raw data output configuration']['density'])
        for i, name in enumerate(density_names):
            density_names[i] = int(name)
        for i in f['raw data output configuration']['density']:
            density_values.append(f['raw data output configuration']['density'][i])
        density_values = [x for y, x in sorted(zip(density_names, density_values))]
        density_names.sort()

    if 'velocity' in list(f['raw data output configuration']):
        velocity_values = []
        velocity_names = list(f['raw data output configuration']['velocity'])
        for i, name in enumerate(velocity_names):
            velocity_names[i] = int(name)
        for i in f['raw data output configuration']['velocity']:
            velocity_values.append(f['raw data output configuration']['velocity'][i])
        velocity_values = [x for y, x in sorted(zip(velocity_names, velocity_values))]
        velocity_names.sort()

if not os.path.exists("temp_png_to_mp4"):
    os.makedirs("temp_png_to_mp4")

pool = Pool()
pool.map(make_density_pictures, range(len(density_values)))
pool.map(make_velocity_pictures, range(len(velocity_values)))

if os.path.isfile("density.mp4"):
    print("file density.mp4 already exists. Please rename it and start again.")
else:
    os.system(
        "ffmpeg -r 30 -i ./temp_png_to_mp4/density_%01d.png -vb 10M ./density.mp4")

if os.path.isfile("velocity.mp4"):
    print("file velocity.mp4 already exists. Please rename it and start again.")
else:
    os.system(
        "ffmpeg -r 30 -i ./temp_png_to_mp4/velocity_%01d.png -vb 10M ./velocity.mp4")

for root, dirs, files in os.walk("temp_png_to_mp4", topdown=False):
    for name in files:
        os.remove(os.path.join(root, name))
os.rmdir("temp_png_to_mp4")
