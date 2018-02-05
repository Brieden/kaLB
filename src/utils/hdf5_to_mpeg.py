from multiprocessing import Pool
import matplotlib.pyplot as plt
import numpy as np
import h5py
import os

plt.style.use('bmh')
plt.rcParams['figure.figsize'] = (16.0, 9.0)


def make_pressure_pictures(number):
    plt.imshow(np.array(pressure_values[number]), origin="lower")
    plt.tight_layout()
    plt.savefig("./temp_png_to_mp4/pressure_%i.png" % number)
    plt.cla()


def make_velocity_pictures(number):
    plt.imshow(np.array(velocity_values[number][0]) +
               np.array(velocity_values[number][1]), origin="lower")
    plt.tight_layout()
    plt.savefig("./temp_png_to_mp4/velocity_%i.png" % number)
    plt.cla()


filename = "LBM_raw_data.hdf5"
f = h5py.File(filename, 'r')
a_group_key = list(f.keys())

if 'raw data output configuration' in a_group_key:

    if 'pressure' in list(f['raw data output configuration']):
        pressure_values = []
        pressure_names = list(f['raw data output configuration']['pressure'])
        for i, name in enumerate(pressure_names):
            pressure_names[i] = int(name)
        for i in f['raw data output configuration']['pressure']:
            pressure_values.append(f['raw data output configuration']['pressure'][i])
        pressure_values = [x for y, x in sorted(zip(pressure_names, pressure_values))]

    if 'velocity' in list(f['raw data output configuration']):
        velocity_values = []
        velocity_names = list(f['raw data output configuration']['velocity'])
        for i, name in enumerate(velocity_names):
            velocity_names[i] = int(name)
        for i in f['raw data output configuration']['velocity']:
            velocity_values.append(f['raw data output configuration']['velocity'][i])
        velocity_values = [x for y, x in sorted(zip(velocity_names, velocity_values))]

if not os.path.exists("temp_png_to_mp4"):
    os.makedirs("temp_png_to_mp4")

pool = Pool()
pool.map(make_pressure_pictures, range(len(pressure_values)))
pool.map(make_velocity_pictures, range(len(velocity_values)))

if os.path.isfile("pressure.mp4"):
    print("file pressure.mp4 already exists. Please rename it and start again.")
else:
    os.system(
        "ffmpeg -f image2 -r 30 -i ./temp_png_to_mp4/pressure_%01d.png -pix_fmt yuv420p -c:v libx264 -vcodec mpeg4 -y ./pressure.mp4")

if os.path.isfile("velocity.mp4"):
    print("file velocity.mp4 already exists. Please rename it and start again.")
else:
    os.system(
        "ffmpeg -f image2 -r 30 -i ./temp_png_to_mp4/velocity_%01d.png -pix_fmt yuv420p -c:v libx264 -vcodec mpeg4 -y ./velocity.mp4")

for root, dirs, files in os.walk("temp_png_to_mp4", topdown=False):
    for name in files:
        os.remove(os.path.join(root, name))
os.rmdir("temp_png_to_mp4")
