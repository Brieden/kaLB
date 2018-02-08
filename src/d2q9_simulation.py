# -*- coding: utf-8 -*-
r"""
kaLB = kaum ausgereiftes Lattice Boltzmann
-------------------------------------------
Corefunction to calculate fluid dynamics
"""
import numpy as np
# matplotlib nur zum debugging
# eigentlich sollte in dieser datei keine notwendigkeit zum plotten bestehen
# -> am besten alle funktionen zum plotten & ansehen der ergebnisse auslagern!
import matplotlib.pyplot as plt
import matplotlib.image as img
import sys
import os
import h5py


class Simulation():
    """
    Docstring
    """
    # directions
    e = np.array([
        [0, 0],
        [1, 0], [0, 1], [-1, 0], [0, -1],
        [1, 1], [-1, 1], [-1, -1], [1, -1]
    ])

    e_inverse = np.asarray([0, 3, 4, 1, 2, 7, 8, 5, 6])

    direction_sets = {
        "N": np.asarray([i for i, e_i in enumerate(e) if e_i[1] > 0]),
        "E": np.asarray([i for i, e_i in enumerate(e) if e_i[0] > 0]),
        "S": np.asarray([i for i, e_i in enumerate(e) if e_i[1] < 0]),
        "W": np.asarray([i for i, e_i in enumerate(e) if e_i[0] < 0])
    }
    
    vertical_indices = np.asarray([i for i, e_i in enumerate(e) if e_i[0] == 0])
    horizontal_indices = np.asarray([i for i, e_i in enumerate(e) if e_i[1] == 0])
    
    # weights
    w = np.array([
        4 / 9,
        1 / 9, 1 / 9, 1 / 9, 1 / 9,
        1 / 36, 1 / 36, 1 / 36, 1 / 36
    ])

    def __init__(self):
        """
        """
        # bisher unbenutzt, soll aber zur kontrolle dienen, ob alle
        # Simulationsparameter zulässig sind.
        # kann nach prüfung dann auf True gesetzt werden
        self.is_initialized = False
        # tau ist hier erstmal nur ein platzhalter und sollte sattdessen
        # in der initialize_from_json funktion aus reynolds berechnet werden.
        self.tau = 10

    def cylinder_function(self, obstacle_parameter):
        """

        :param obstacle_parameter:
        :return:
        """
        cylinder_y = obstacle_parameter["y-position"]
        cylinder_x = obstacle_parameter["x-position"]
        cylinder_r = obstacle_parameter["radius"]

        return np.fromfunction(
            lambda x, y: (cylinder_x - x) ** 2 + (cylinder_y - y) ** 2 < cylinder_r ** 2,
            self.shape)
    
    def recktangle_function(self, obstacle_parameter):
        """

        :param obstacle_parameter:
        :return:
        """
        point1 = obstacle_parameter["bottom_left"]
        point2 = obstacle_parameter["top_right"]
        recktangle = np.full(shape=self.shape, fill_value=False)
        recktangle[point1[0]:point2[0] + 1, point1[1]:point2[1] + 1] = True
        return recktangle

    def png_importer(self, png_path):
        """
        :param png_path:
        :return:
        """
        try:
            image = np.rot90(np.flipud(img.imread(png_path)[:, :, :-1].sum(axis=2)))
        except IOError as error:
            print("ERROR: It was not possible to read the picture: " + str(png_path))
            quit()
        if image.shape == self.shape:
            return image[:] < 1
        print("The shape of the picture %s does not match the shape of the simulation. "
              "Simulation was aborted." % png_path)
        quit()

    def obstacles_definition(self, obstacle_parameters):
        """
        Function to combine the management and association of different obstacles.
        In the loop over the obstacle can be added further types of them.

        :param self:
        :param obstacle_parameters:
        :return: A boolean matrix in the dimension of the simulation grid with a "true" at
        the points where there is an obstacle.
        """
        self.obstacle = np.full(shape=self.shape, fill_value=False)
        for obstacle_parameter in obstacle_parameters:
            if obstacle_parameter["type"] == "cylindrical obstacle":
                self.obstacle = np.logical_or(self.obstacle,
                                              self.cylinder_function(obstacle_parameter))
            elif obstacle_parameter["type"] == "recktangle obstacle":
                self.obstacle = np.logical_or(self.obstacle,
                                              self.recktangle_function(obstacle_parameter))
            elif obstacle_parameter["type"] == "png import":
                self.obstacle = np.logical_or(self.obstacle,
                                              self.png_importer(obstacle_parameter["file name"]))
            else:
                print("obstacle ", obstacle_parameter, "not recognised")
                quit()

        if self.args.verbose:
            plt.imshow(self.obstacle.T, origin='lower', cmap='Greys', interpolation='nearest')
            plt.show()

    def initialize_output(self, output_parameters):
        """
        :param output_parameters:
        :return:
        """
        self.raw_output = False
        self.picture_output = False
        self.snapshot = False

        if not os.path.exists(self.args.output):
            os.makedirs(self.args.output)

        if "snapshot" in output_parameters:
            if not os.path.exists(self.args.output + "snapshots"):
                os.makedirs(self.args.output + "snapshots")
            self.snapshot = True
            self.snapshot_frequency = output_parameters["snapshot"]["output frequency"]

        if "picture output configuration" in output_parameters:
            picture_parameter = output_parameters["picture output configuration"]
            self.picture_output = True
            self.picture_output_frequency = picture_parameter["output frequency"]
            self.picture_output_typ = picture_parameter["file type"]
            self.picture_output_name = picture_parameter["file name"]

        if "raw data output configuration" in output_parameters:
            raw_parameter = output_parameters["raw data output configuration"]
            self.raw_output = True
            self.raw_output_frequency = raw_parameter["output frequency"]
            h5file = h5py.File(self.args.output + raw_parameter["file name"]
                               + ".hdf5", "w")
            h5_output = h5file.create_group("raw data output configuration")
            if raw_parameter["output data"]["velocity"] == 1:
                self.h5_velocity = h5_output.create_group("velocity")
            if raw_parameter["output data"]["pressure"] == 1:
                self.h5_pressure = h5_output.create_group("pressure")

    def store_output(self, step):
        """
        :param step:
        :return:
        """
        if self.snapshot:
            if step % self.snapshot_frequency == 0:
                np.save(self.args.output + "snapshots/snap_%05i.npz" % step, self.f_in)
        if self.raw_output:
            if step % self.raw_output_frequency == 0:
                self.h5_velocity.create_dataset("%i" % step, data=self.vel)
                self.h5_pressure.create_dataset("%i" % step, data=self.rho)
        if self.picture_output:
            if step % self.picture_output_frequency == 0:
                plt.imshow((self.vel[0] * self.vel[0] + self.vel[1] * self.vel[1]).T, origin='lower')
                plt.savefig(self.args.output + self.picture_output_name + "%05i." % step + self.picture_output_typ)

    def initialize_from_json(self, inputfile, args):
        """
        """
        # hier soll alles aus der json datei ausgelsen
        # und in die simulation übertragen werden.
        #
        # ggf. könnte die funktion zur überschaubarkeit
        # in unterfunktionen aufgeteilt werden.
        #
        # zusätzlich könnte hier dann auch auf korrektheit überprüft werden
        # aber vermutlich ist eine separate funktion dafür sinnvoller,
        # oder kann mit obiger idee kombiniert werden.
        #
        # zunächst werden simulationsparameter gelesen und gesetzt.
        # achtung, json datei ist leicht umgebaut!
        simulation_parameters = inputfile["simulation parameters"]
        self.name = simulation_parameters["simulation name"]
        self.id = simulation_parameters["simulation id"]
        self.timesteps = simulation_parameters["time steps"]
        self.n_x = simulation_parameters["lattice points x"]
        self.n_y = simulation_parameters["lattice points y"]
        self.reynolds = simulation_parameters["Reynolds number"]
        self.args = args

        # für die Simulation benötigten Felder initialisieren
        self.shape = (self.n_x, self.n_y)
        self.rho = np.empty(shape=self.shape)
        self.vel = np.empty(shape=(2, self.shape[0], self.shape[1]))
        self.f_in = np.empty(shape=(9, self.shape[0], self.shape[1]))
        self.f_eq = np.empty_like(self.f_in)
        self.f_out = np.empty_like(self.f_in)

        self.obstacles_definition(inputfile["obstacle parameters"])
        
        self.set_boundary_conditions(inputfile["boundary conditions"])

        # umgang mit den output-parametern:
        # also alles was mit speichern, plotten,
        # oder ausgabe während des programms zu tun hat
        self.initialize_output(inputfile["output configuration"])

    def set_boundary_conditions(self, boundary_conditions):
        directions = np.array(["N", "E", "S", "W"])
        self.opposite_directions = {"N": "S", "E": "W", "S": "N", "W": "E"}
        self.last_indices = {"N": (-1,-2), "E": (-1,-2), "S": (0,1), "W": (0,1)}
        self.boundarys = {}
        self.zou_he_conditions = {}
        
        # TODO zeug
        for direction in directions:
            bc = boundary_conditions[direction]
            if bc["type"] == "periodic":
                oppo_bc = boundary_conditions[self.opposite_directions[direction]]
                if oppo_bc["type"] == "periodic":
                    self.boundarys[direction] = "periodic"
                else:
                    print("ERROR: Periodic boundary conditions do not match")
                    quit()
            elif bc["type"] == "bounce_back":
                self.boundarys[direction] = "bounce_back"                
                if direction == "N" or direction == "S":
                    self.obstacle[:, self.last_indices[direction][0]] = True
                elif direction == "E" or direction == "W":
                    self.obstacle[self.last_indices[direction][0], :] = True
                else:
                    print("ERROR: This state should be impossible!")
                    quit()
            elif bc["type"] == "outflow":
                self.boundarys[direction] = "outflow"
            elif bc["type"] == "zou-he":
                self.boundarys[direction] = "zou-he"
                self.zou_he_conditions[direction] = (bc["v_x"], bc["v_y"])
            else:
                print("ERROR: Boundary condition does not exist or is missing")
                quit()

    def progress_bar(self, value, endvalue, bar_length=50):
        """
        :param value:
        :param endvalue:
        :param bar_length:
        :return:
        """
        percent = float(value) / endvalue
        arrow = '-' * int(round(percent * bar_length) - 1) + '>'
        spaces = ' ' * (bar_length - len(arrow))

        sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()

    def calc_macroscopic(self):
        """
        Calculate macroscopic density and velocity.

        Take the distribution function :math:`f_i`
        and direction-vectors :math:`\\vec{e}_i`
        and compute related density :math:`\\rho`
        and velocity :math:`\\vec{u}`:

            .. math::
                \\rho(\\vec{x}) = \sum_{i} f_{i}(\\vec{x})
            .. math::
                \\vec{u}(\\vec{x}) = \\frac{1}{\\rho}\sum_{i} f_{i}\\vec{e}_{i}
        """
        self.rho = np.sum(self.f_in, axis=0)
        self.vel = np.dot(self.e.T, self.f_in.transpose((1, 0, 2))) / self.rho

    def calc_equilibrium(self):
        """
        Calculate equilibrium distribution fuction :math:`f_{eq}`.

        Take macroscopic density :math:`\\rho` and velocity :math:`\\vec{u}`
        together with direction-vectors :math:`\\vec{e}_i`
        and their according weight-vectors :math:`w_i`
        to calculate equilibrium distribution fuction :math:`f_{eq}`
        using Bhatnagar-Gross-Krook (BGK) collision, with :math:`c = 1`:

            .. math::
                f_i^{eq}(\\vec{x}) = w_i\\rho(\\vec{x}) * \\left(
                1
                + 3 \\frac{(\\vec{e}_i * \\vec{u})}{c}
                + \\frac{9}{2} \\frac{(\\vec{e}_i * \\vec{u})^2}{c^2}
                - \\frac{3}{2} \\frac{(\\vec{u} * \\vec{u})}{c^2}
                \\right)
        """
        vel_sqared = self.vel[0] * self.vel[0] + self.vel[1] * self.vel[1]
        e_n_x_vel = np.dot(self.e, self.vel.transpose(1, 0, 2))
        for i in range(9):
            self.f_eq[i] = self.rho * self.w[i] * (
                    1
                    + (3 * e_n_x_vel[i])
                    + (4.5 * e_n_x_vel[i] ** 2)
                    - (1.5 * vel_sqared)
            )

    def collision_step(self):
        """
        Perform collision-step to update distribution function f_out.

        Use distribution function :math:`f_i`,
        equilibrium distribution function :math:`f_i^{eq}`
        and relaxation parameter :math:`\\tau`
        to calculate the updated distribution function :math:`f_i^*`:

            .. math::
                f_i^* = f_i - \\frac{1}{\\tau} (f_i - f_i^{eq})
        """
        self.f_out = self.f_in - (self.f_in - self.f_eq) / self.tau

    def stream_step(self):
        """
        Perform streaming-step and produce the shifted distribution function.

        Use distribution function :math:`f_i`
        and direction-vectors :math:`\\vec{e}_i`
        to calculate shifted distribution function :math:`f_i^*`:

            .. math::
                f_i^*(\\vec{x}) = f_i(\\vec{x} + \\vec{e}_i)
        """
        for i in range(9):
            self.f_in[i] = np.roll(self.f_out[i], shift=self.e[i], axis=(0, 1))

    def correct_macroscopic(self):
        """
        Correct incoming density and velocity at a given border

        Take macroscopic density :math:`\\rho` and velocity :math:`\\vec{u}`
        and set velocity at the given border to given_instream
        and calculate according density at given border.
        """
        # TODO fulfill docstring
        for direction, condition in self.boundarys.items():
            if condition == "zou-he":
                
                last = self.last_indices[direction][0]
                v_x, v_y = self.zou_he_conditions[direction]
                direction_indices = self.direction_sets[direction]

                if direction == "N" or direction == "S":
                    self.vel[0, :, last] = v_x
                    self.vel[1, :, last] = v_y

                    f1 = np.sum(self.f_in[self.horizontal_indices, :, last], axis=0)
                    f2 = np.sum(self.f_in[direction_indices, :, last], axis=0)

                    self.rho[:, last] = f1 + (2 * f2) / (1 - self.vel[1, :, last])

                elif direction == "E" or direction == "W":
                    self.vel[0, last, :] = v_x
                    self.vel[1, last, :] = v_y

                    f1 = np.sum(self.f_in[self.vertical_indices, last, :], axis=0)
                    f2 = np.sum(self.f_in[direction_indices, last, :], axis=0)

                    self.rho[last, :] = f1 + (2 * f2) / (1 - self.vel[0, last, :])

    def correct_distr_func(self):
        """
        Correct incoming components of distribution function at Zou-He Border

        Take the current distribution function :math:`f_i`
        and calculate corrected values
        for incoming components of the distribution function using Zou-He.
        """
        # TODO fulfill docstring
        for direction, condition in self.boundarys.items():
            if condition == "zou-he":

                last = self.last_indices[direction][0]                
                opposite_direction = self.opposite_directions[direction]
                opposite_direction_indices = self.direction_sets[opposite_direction]
                matching_direction_indices = self.e_inverse[opposite_direction_indices]

                if direction == "N" or direction == "S":
                    self.f_in[opposite_direction_indices, :, last] = (
                        self.f_eq[opposite_direction_indices, :, last]
                        + self.f_in[matching_direction_indices, :, last]
                        - self.f_eq[matching_direction_indices, :, last]
                    )
                elif direction == "E" or direction == "W":
                    self.f_in[opposite_direction_indices, last, :] = (
                        self.f_eq[opposite_direction_indices, last, :]
                        + self.f_in[matching_direction_indices, last, :]
                        - self.f_eq[matching_direction_indices, last, :]
                    )

    def bounce_back(self):
        """
        Apply bounce back boundary condition (BC).

        Take the current distribution function :math:`f_i`
        and the distribution function after collision :math:`f_i^*`
        together with an array that specifies all obstacles
        and a list of indices *e_inverse*
        that correspond to the inverse directions of :math:`\\vec{e}_i`.\n
        Compute a distribution function :math:`f_i^{BC}` that satisfies BC
        by copying components of the distribution function :math:`f_i`
        - which would stream into the obstacle -
        into their inverse direction
        of the post-collision distribution function.
        """
        for i, j in enumerate(self.e_inverse):
            self.f_out[i, self.obstacle] = self.f_in[j, self.obstacle]

    def correct_outflow(self):
        """
        Correct distribution function

        Take post-streaming distribution function :math:`f_i`
        and replace incorrectly incoming components after streaming
        with second to last entries from that same distribution function.
        """
        # TODO fulfill docstring
        for direction, condition in self.boundarys.items():
            if condition == "outflow":

                opposite_direction = self.opposite_directions[direction]
                opposite_direction_indices = self.direction_sets[opposite_direction]
                last, second_to_last = self.last_indices[direction]

                if direction == "N" or direction == "S":
                    self.f_in[opposite_direction_indices, :, last] = self.f_in[opposite_direction_indices, :, second_to_last]
                elif direction == "E" or direction == "W":
                    self.f_in[opposite_direction_indices, last, :] = self.f_in[opposite_direction_indices, second_to_last, :]

    def do_simulation_step(self, step):
        self.calc_macroscopic()
        self.correct_macroscopic()
        self.calc_equilibrium()
        self.correct_distr_func()
        self.collision_step()
        self.bounce_back()
        self.stream_step()
        self.correct_outflow()
        self.store_output(step)

    def run_simulation(self):
        # initialisiere ein anfangs-geschwindigkeitsfeld.
        # das solle der anwender vorgeben!
        # wenn keines vorgegeben ist, sollte ein default generiert werden.
        # die 3 zeilen sollten also langfristig nicht hier,
        # sondern beim initialisieren stehen....
        self.vel[:] = 0
#        self.vel[0, 1, :] = 0.04 # links nach rechts
#        self.vel[0, -1, :] = -0.04 # rechts nach links
        self.vel[1, :, 0] = 0.04 # unten nach oben
#        self.vel[1, :, -1] = -0.04 # oben nach unten
        self.rho[:] = 1

        # errechne ein initiales f_in.
        # hier sollte langfristig auch die möglichkeit existieren
        # eine rehcnung mit einem existierenden f_in
        # fortzusetzen oder zu starten!
        # für den moment ist das aber erstmal okay so.
        self.calc_equilibrium()
        self.f_in = self.f_eq

        for step in range(self.timesteps):
            self.do_simulation_step(step)
            if self.args.verbose:
                self.progress_bar(step, self.timesteps)
