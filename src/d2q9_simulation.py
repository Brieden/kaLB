# -*- coding: utf-8 -*-
r"""
kaLB = kaum ausgereiftes Lattice Boltzmann
-------------------------------------------
Corefunction to calculate fluid dynamics
"""
import numpy as np
from src import utilities
import time


class Simulation():
    """
    Simulation class
    
    Class to hold Simulation state and perform simulation
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

    def __init__(self, inputfile, args):

        self.args = args
        utilities.simulation_parameters_definition(self, inputfile["simulation parameters"])
        utilities.obstacles_definition(self, inputfile["obstacle parameters"])
        utilities.set_boundary_conditions(self, inputfile["boundary conditions"])
        utilities.initialize_output(self, inputfile["output configuration"])

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
        e_n_x_vel = np.dot(self.e, self.vel.transpose(1, 0, 2))
        for i in range(9):
            self.f_eq[i] = self.rho * self.w[i] * (
                    1
                    + (3 * e_n_x_vel[i])
                    + (4.5 * e_n_x_vel[i] ** 2)
                    - (1.5 * vel_sqared)
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
        for i, j in enumerate(self.e_inverse):
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
        for direction, condition in self.boundarys.items():
            if condition == "outflow":

                opposite_direction = self.opposite_directions[direction]
                opposite_direction_indices = self.direction_sets[opposite_direction]
                last, second_to_last = self.last_indices[direction]

                if direction == "N" or direction == "S":
                    self.f_in[opposite_direction_indices, :, last] = self.f_in[opposite_direction_indices, :,
                                                                     second_to_last]
                elif direction == "E" or direction == "W":
                    self.f_in[opposite_direction_indices, last, :] = self.f_in[opposite_direction_indices,
                                                                     second_to_last, :]

    def do_simulation_step(self, step):
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
        utilities.store_output(self, step)

    def run_simulation(self):
        """
        Start simulation
        """
        self.vel[:] = 0
        self.rho[:] = 1
        self.calc_equilibrium()
        self.f_in = self.f_eq
        if self.args.snapshot:
            try:
                self.f_in = np.load(self.args.snapshot)
            except IOError as error:
                print("Error: could not open snapshot-file. Simulation was aborted.")
                quit()
        t_0 = time.time()
        for step in range(self.timesteps):
            self.do_simulation_step(step)
            utilities.progress_bar(step, self.timesteps)

        if self.args.verbose:
            print("\n Performance feedback: \n"
                  "%i steps in %2.f seconds with %.2f M su/s"
                  % (self.timesteps, time.time() - t_0,
                     self.n_x * self.n_y * self.timesteps * 1e-6 / (time.time() - t_0)))
