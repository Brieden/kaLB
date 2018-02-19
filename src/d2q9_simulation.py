# -*- coding: utf-8 -*-
"""
This file holds the Simulation class.\n
The Simulation is not intendet to run without its utilities, tho.
"""
import time
import numpy as np
from src import utilities


class Simulation():
    """
    Simulation class

    Class to hold Simulation state and perform simulation.
    It also holds calss variables that are necessary for a d2q9 simulation;
    these are shared by all instances.
    initialization methods for json-input are outsourced to utilities.

    All method-calls for the simulation are performed on the object instance.
    """

    #: directions-array
    e = np.array([
        [0, 0],
        [1, 0], [0, 1], [-1, 0], [0, -1],
        [1, 1], [-1, 1], [-1, -1], [1, -1]
    ])

    #: indices to get inverse directions: " e[e_inverse] = -e "
    e_inverse = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

    #: indices grouped by directions -> 3 indices per cardinal direction
    direction_sets = {
        "N": np.asarray([i for i, e_i in enumerate(e) if e_i[1] > 0]),
        "E": np.asarray([i for i, e_i in enumerate(e) if e_i[0] > 0]),
        "S": np.asarray([i for i, e_i in enumerate(e) if e_i[1] < 0]),
        "W": np.asarray([i for i, e_i in enumerate(e) if e_i[0] < 0])
    }

    #: indices for directions-array that have no x component
    vertical_indices = np.asarray([i for i, e_i in enumerate(e) if e_i[0] == 0])

    #: indices for directions-array that have no y component
    horizontal_indices = np.asarray([i for i, e_i in enumerate(e) if e_i[1] == 0])

    #: weights according to directions
    w = np.array([
        4 / 9,
        1 / 9, 1 / 9, 1 / 9, 1 / 9,
        1 / 36, 1 / 36, 1 / 36, 1 / 36
    ])

    def __init__(self, inputfile=None, args=None):
        """
        Initialize an instance of Simulation

        Checks for arguments.
        Start initialization if both arguments are given.
        Create mockup object for testing purposes if nothing is given.

        :param inputfile:
            an already opened .json file that specifies simulation parameters

        :param args:
            commandline arguments that are given to kaLB
        """

        # inputfile and argumends are given -> initialize properly
        if (inputfile is not None) and (args is not None):
            self.args = args
            utilities.simulation_parameters_definition(self, inputfile["simulation parameters"])
            utilities.obstacles_definition(self, inputfile["obstacle parameters"])
            utilities.set_boundary_conditions(self, inputfile["boundary conditions"])
            utilities.initialize_output(self, inputfile["output configuration"])

        # create mockup Simulation
        elif (inputfile is None) and (args is None):
            pass

        # Error
        else:
            print("ERROR: inputfile and args should be given for initialization")
            quit()

    def prepare_simulation(self):
        """
        pre-iteration: Set initial distribution function.

        If snapshot is loaded:
            Load initial distribution function from snapshot
        Else:
            Set initial macroscopic values and
            calculate initial distribution function
        """

        if self.args.snapshot:
            try:
                self.f_in = np.load(self.args.snapshot)
                if self.step_offset == 0:
                    print(
                        """
                        WARNING: Starting from a snapshot but 'step offset' was not changed.\n
                        This will cause duplication in naming the steps in your output.\n
                        """
                    )
            except IOError as error:
                print("Error: could not open snapshot-file. Simulation was aborted.\n" + str(error))
                quit()
        else:
            self.vel[:] = 0
            self.rho[:] = 1
            self.calc_equilibrium()
            self.f_in = self.f_eq

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

    def correct_macroscopic(self):
        """
        Correct incoming density and velocity at *zou-he*-borders

        Take macroscopic density :math:`\\rho` and velocity :math:`\\vec{u}`
        and set velocity at borders with *zou-he* boundary condition
        to fixed values specified in inputfile;
        and calculate according density at these borders.
        """

        # find all borders that are set to zou-he condition
        for direction, condition in self.boundarys.items():
            if condition == "zou-he":

                # setup
                last = self.last_indices[direction][0]
                v_x, v_y = self.zou_he_conditions[direction]
                direction_indices = self.direction_sets[direction]

                # horizontal case
                if direction == "N" or direction == "S":
                    # set velocities at border to fixed values (from inputfile)
                    self.vel[0, :, last] = v_x
                    self.vel[1, :, last] = v_y

                    # calculate according density at border
                    f1 = np.sum(self.f_in[self.horizontal_indices, :, last], axis=0)
                    f2 = np.sum(self.f_in[direction_indices, :, last], axis=0)
                    self.rho[:, last] = f1 + (2 * f2) / (1 - self.vel[1, :, last])

                # vertical case
                elif direction == "E" or direction == "W":
                    # set velocities at border to fixed values (from inputfile)
                    self.vel[0, last, :] = v_x
                    self.vel[1, last, :] = v_y

                    # calculate according density at border
                    f1 = np.sum(self.f_in[self.vertical_indices, last, :], axis=0)
                    f2 = np.sum(self.f_in[direction_indices, last, :], axis=0)
                    self.rho[last, :] = f1 + (2 * f2) / (1 - self.vel[0, last, :])

    def calc_equilibrium(self):
        """
        Calculate equilibrium distribution fuction :math:`f_{eq}`.

        Take macroscopic density :math:`\\rho` and velocity :math:`\\vec{u}`
        together with direction-vectors :math:`\\vec{e}_i`
        and their according weight-vectors :math:`w_i`
        to calculate equilibrium distribution fuction :math:`f_{eq}`
        using Bhatnagar-Gross-Krook (BGK) collision, with :math:`c = 1`:

            .. math::
                f_i^{eq}(\\vec{x}) = w_i\\rho(\\vec{x}) \\left(
                1
                + 3 \\frac{(\\vec{e}_i \cdot \\vec{u})}{c}
                + \\frac{9}{2} \\frac{(\\vec{e}_i \cdot \\vec{u})^2}{c^2}
                - \\frac{3}{2} \\frac{(\\vec{u} \cdot \\vec{u})}{c^2}
                \\right)
        """

        vel_sqared = self.vel[0] * self.vel[0] + self.vel[1] * self.vel[1]
        for i in range(9):  # loop is actually faster than NumPy for large arrays
            e_n_x_vel = self.e[i, 0] * self.vel[0] + self.e[i, 1] * self.vel[1]
            self.f_eq[i] = self.rho * self.w[i] * (
                1 +
                (3 * e_n_x_vel) +
                (4.5 * e_n_x_vel * e_n_x_vel) -
                (1.5 * vel_sqared)
            )

    def correct_distr_func(self):
        """
        Correct unknown components of distribution function at *zou-he* borders

        Take the current distribution function :math:`f_i`
        and calculate corrected values
        for incoming components using zou-he:

            .. math::
                f_i^{corrected} := f_i^{eq} + f_j - f_j^{eq}

        for all 3 incoming :math:`\\vec{e}_i` at the border,
        with *i*, *j* satisfying :math:`\\vec{e}_i = -\\vec{e}_j`.
        """

        # find all borders that are set to zou-he condition
        for direction, condition in self.boundarys.items():
            if condition == "zou-he":

                # setup
                last = self.last_indices[direction][0]
                opposite_direction = self.opposite_directions[direction]
                opposite_direction_indices = self.direction_sets[opposite_direction]
                matching_direction_indices = self.e_inverse[opposite_direction_indices]

                # horizontal case
                if direction == "N" or direction == "S":
                    self.f_in[opposite_direction_indices, :, last] = (
                        self.f_eq[opposite_direction_indices, :, last] +
                        self.f_in[matching_direction_indices, :, last] -
                        self.f_eq[matching_direction_indices, :, last]
                    )

                # vertical case
                elif direction == "E" or direction == "W":
                    self.f_in[opposite_direction_indices, last, :] = (
                        self.f_eq[opposite_direction_indices, last, :] +
                        self.f_in[matching_direction_indices, last, :] -
                        self.f_eq[matching_direction_indices, last, :]
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

    def bounce_back(self):
        """
        Apply bounce back (no-slip) boundary condition at obstacles.

        Take the current distribution function :math:`f_i`
        together with the obstacle-array and a list of indices
        that correspond to the inverse directions of :math:`\\vec{e}_i`.

        Correct the distribution function after collision :math:`f_i^*`
        by copying components of the distribution function :math:`f_i`
        - which would stream into the obstacle (:math:`\\vec{x} = obstacle`) -
        into their inverse direction
        of the post-collision distribution function :math:`f_i^*`:

            .. math::
                f_i^*(\\vec{x} = obstacle) := f_j(\\vec{x} = obstacle)

        with *i*, *j* satisfying :math:`\\vec{e}_i = -\\vec{e}_j`.
        """

        for i, j in enumerate(self.e_inverse):  # problem with NumPy. see Issue #9
            self.f_out[i, self.obstacle] = self.f_in[j, self.obstacle]

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

    def correct_outflow(self):
        """
        Correct distribution function

        Take post-streaming distribution function :math:`f_i`
        and replace incorrectly incoming components (due to streaming)
        at an *outflow*-border with second-to-last entries
        (one gridpoint away from the border)
        from that same distribution function.
        """

        # find all borders that are set to outflow condition
        for direction, condition in self.boundarys.items():
            if condition == "outflow":

                # setup
                opposite_direction = self.opposite_directions[direction]
                opposite_direction_indices = self.direction_sets[opposite_direction]
                last, second_to_last = self.last_indices[direction]

                # horizontal case
                if direction == "N" or direction == "S":
                    self.f_in[opposite_direction_indices, :, last] = (
                        self.f_in[opposite_direction_indices, :, second_to_last]
                    )
                # vertical case
                elif direction == "E" or direction == "W":
                    self.f_in[opposite_direction_indices, last, :] = (
                        self.f_in[opposite_direction_indices, second_to_last, :]
                    )

    def do_simulation_step(self):
        """
        Perform an iteration step with the current Simulation
        """

        self.calc_macroscopic()
        self.correct_macroscopic()
        self.calc_equilibrium()
        self.correct_distr_func()
        self.collision_step()
        self.bounce_back()
        self.stream_step()
        self.correct_outflow()

    def run_simulation(self):
        """
        Start a simulation
        """

        self.prepare_simulation()

        # main simulation loop
        t0 = time.time()
        for step in range(1, self.timesteps + 1):
            self.do_simulation_step()
            utilities.store_output(self, step)
            if not self.args.no_progessbar:
                utilities.progress_bar(step - 1, self.timesteps)
        t1 = time.time()

        # performance feedback
        if self.args.performance_feedback:
            utilities.give_powerfeedback(self, t0, t1)
