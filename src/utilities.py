import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
import h5py
import sys
import os


def simulation_parameters_definition(sim, simulation_parameters):
    sim.name = simulation_parameters["simulation name"]
    sim.id = simulation_parameters["simulation id"]
    sim.timesteps = simulation_parameters["time steps"]
    sim.step_offset = simulation_parameters["step offset"]
    sim.n_x = simulation_parameters["lattice points x"]
    sim.n_y = simulation_parameters["lattice points y"]
    sim.tau = simulation_parameters["tau"]

    sim.shape = (sim.n_x, sim.n_y)
    sim.rho = np.empty(shape=sim.shape)
    sim.vel = np.empty(shape=(2, sim.shape[0], sim.shape[1]))
    sim.f_in = np.empty(shape=(9, sim.shape[0], sim.shape[1]))
    sim.f_eq = np.empty_like(sim.f_in)
    sim.f_out = np.empty_like(sim.f_in)


def obstacles_definition(sim, obstacle_parameters):
    """
    Function to combine the management and association of different obstacles.
    In the loop over the obstacle can be added further types of them.

    :param sim:
    :param obstacle_parameters:
    :return: A boolean matrix in the dimension of the simulation grid with a "true" at
    the points where there is an obstacle.
    """
    sim.obstacle = np.full(shape=sim.shape, fill_value=False)
    for obstacle_parameter in obstacle_parameters:
        if obstacle_parameter["type"] == "cylindrical obstacle":
            sim.obstacle = np.logical_or(sim.obstacle,
                                         cylinder_function(sim, obstacle_parameter))
        elif obstacle_parameter["type"] == "recktangle obstacle":
            sim.obstacle = np.logical_or(sim.obstacle,
                                         recktangle_function(sim, obstacle_parameter))
        elif obstacle_parameter["type"] == "png import":
            sim.obstacle = np.logical_or(sim.obstacle,
                                         png_importer(sim, obstacle_parameter["file name"]))
        else:
            print("obstacle ", obstacle_parameter, "not recognised")
            quit()

    if sim.args.verbose:
        plt.imshow(sim.obstacle.T, origin='lower', cmap='Greys', interpolation='nearest')
        plt.show()


def cylinder_function(sim, obstacle_parameter):
    """

    :param obstacle_parameter:
    :return:
    """
    cylinder_y = obstacle_parameter["y-position"]
    cylinder_x = obstacle_parameter["x-position"]
    cylinder_r = obstacle_parameter["radius"]

    return np.fromfunction(
        lambda x, y: (cylinder_x - x) ** 2 + (cylinder_y - y) ** 2 < cylinder_r ** 2,
        sim.shape)


def recktangle_function(sim, obstacle_parameter):
    """

    :param obstacle_parameter:
    :return:
    """
    point1 = obstacle_parameter["bottom_left"]
    point2 = obstacle_parameter["top_right"]
    recktangle = np.full(shape=sim.shape, fill_value=False)
    recktangle[point1[0]:point2[0] + 1, point1[1]:point2[1] + 1] = True
    return recktangle


def png_importer(sim, png_path):
    """
    :param png_path:
    :return:
    """
    try:
        image = img.imread(png_path)[:, :, :-1].sum(axis=2)
        image = np.rot90(np.flipud(np.fliplr(image)))
    except IOError as error:
        print("ERROR: It was not possible to read the picture: " + str(png_path))
        quit()
    if image.shape == sim.shape:
        return image[:] < 1
    print("The shape of the picture %s does not match the shape of the simulation. "
          "Simulation was aborted." % png_path)
    quit()


def set_boundary_conditions(sim, boundary_conditions):
    directions = np.array(["N", "E", "S", "W"])
    sim.opposite_directions = {"N": "S", "E": "W", "S": "N", "W": "E"}
    sim.last_indices = {"N": (-1, -2), "E": (-1, -2), "S": (0, 1), "W": (0, 1)}
    sim.boundarys = {}
    sim.zou_he_conditions = {}

    # TODO zeug
    for direction in directions:
        bc = boundary_conditions[direction]
        if bc["type"] == "periodic":
            oppo_bc = boundary_conditions[sim.opposite_directions[direction]]
            if oppo_bc["type"] == "periodic":
                sim.boundarys[direction] = "periodic"
            else:
                print("ERROR: Periodic boundary conditions do not match")
                quit()
        elif bc["type"] == "bounce_back":
            sim.boundarys[direction] = "bounce_back"
            if direction == "N" or direction == "S":
                sim.obstacle[:, sim.last_indices[direction][0]] = True
            elif direction == "E" or direction == "W":
                sim.obstacle[sim.last_indices[direction][0], :] = True
            else:
                print("ERROR: This state should be impossible!")
                quit()
        elif bc["type"] == "outflow":
            sim.boundarys[direction] = "outflow"
        elif bc["type"] == "zou-he":
            sim.boundarys[direction] = "zou-he"
            sim.zou_he_conditions[direction] = (bc["v_x"], bc["v_y"])
        else:
            print("ERROR: Boundary condition does not exist or is missing")
            quit()


def initialize_output(sim, output_parameters):
    """
    :param output_parameters:
    :return:
    """
    sim.raw_output = False
    sim.picture_output = False
    sim.snapshot = False

    if not os.path.exists(sim.args.output):
        os.makedirs(sim.args.output)

    if "snapshot" in output_parameters:
        if not os.path.exists(sim.args.output + "snapshots"):
            os.makedirs(sim.args.output + "snapshots")
        sim.snapshot = True
        sim.snapshot_frequency = output_parameters["snapshot"]["output frequency"]

    if "picture output configuration" in output_parameters:
        picture_parameter = output_parameters["picture output configuration"]
        sim.picture_output = True
        sim.picture_output_frequency = picture_parameter["output frequency"]
        sim.picture_output_typ = picture_parameter["file type"]
        sim.picture_output_name = picture_parameter["file name"]

    if "raw data output configuration" in output_parameters:
        raw_parameter = output_parameters["raw data output configuration"]
        sim.raw_output = True
        sim.raw_output_frequency = raw_parameter["output frequency"]
        h5file = h5py.File(sim.args.output + raw_parameter["file name"]
                           + ".hdf5", "w")
        h5_output = h5file.create_group("raw data output configuration")
        sim.h5_velocity = h5_output.create_group("velocity")
        sim.h5_density = h5_output.create_group("density")


def store_output(sim, step):
    """
    :param step:
    :return:
    """
    if sim.snapshot:
        if step % sim.snapshot_frequency == 0:
            np.save(sim.args.output + "snapshots/snap_%05i" % (step + sim.step_offset), sim.f_in)
    if sim.raw_output:
        if step % sim.raw_output_frequency == 0:
            sim.h5_velocity.create_dataset("%i" % (step + sim.step_offset), data=sim.vel)
            sim.h5_density.create_dataset("%i" % (step + sim.step_offset), data=sim.rho)
    if sim.picture_output:
        if step % sim.picture_output_frequency == 0:
            plt.imshow((sim.vel[0] * sim.vel[0] + sim.vel[1] * sim.vel[1]).T, origin='lower', vmin=0,
                       vmax=0.004)
            plt.title("t = %i" % (step + sim.step_offset))
            plt.xlabel("x")
            plt.ylabel("y")
            plt.savefig(sim.args.output + sim.picture_output_name + "%05i." % (step + sim.step_offset) + sim.picture_output_typ)
            plt.cla()


def progress_bar(value, endvalue, bar_length=50):
    """
    :param value:
    :param endvalue:
    :param bar_length:
    :return:
    """
    if value % 100 == 0:
        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length) - 1) + '>'
        spaces = ' ' * (bar_length - len(arrow))

        sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()


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
        print("\nThe system test flow in the pipe was passed.")
    else:
        print("\nThe system test flow in the pipe was not passed.")
    if args.verbose:
        plt.plot(speed, label="Values of the simulation.")
        b = np.poly1d(coefs)
        plt.plot(x, b(x), "x", label="parabolic fit as a reference")
        plt.title("Velocity profile at the end of a tube")
        plt.xlabel("Position x: at right angles to the wall")
        plt.ylabel("Speed")
        plt.legend()
        plt.show()
