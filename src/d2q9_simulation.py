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
from PIL import Image
import matplotlib.image as img


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

    e_inverse = [0, 3, 4, 1, 2, 7, 8, 5, 6]

    direction_sets = {
        "N": np.arange(9)[np.asarray([e_i[1] > 0 for e_i in e])],
        "E": np.arange(9)[np.asarray([e_i[0] > 0 for e_i in e])],
        "S": np.arange(9)[np.asarray([e_i[1] < 0 for e_i in e])],
        "W": np.arange(9)[np.asarray([e_i[0] < 0 for e_i in e])]
    }

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
            lambda x, y: (x - cylinder_x) ** 2 + (y - cylinder_y) ** 2 < cylinder_r ** 2,
            (self.n_x, self.n_y)).T

    def png_importer(self, png_path):
        """

        :param png_path:
        :return:
        """
        try:
            image = img.imread(png_path)[:, :, :-1].sum(axis=2)
        except IOError as error:
            print("ERROR: It was not possible to read the picture: " + str(png_path))

        if image.shape == self.shape:
            return image[:] < 1
        print("The shape of the picture %s does not match the shape of the simulation. "
              "Simulation was aborted." %png_path)
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
            elif obstacle_parameter["type"] == "png import":
                self.obstacle = np.logical_or(self.obstacle,
                                              self.png_importer(obstacle_parameter["file name"]))
            else:
                print("obstacle ", obstacle_parameter, "not recognised")
                quit()

        self.obstacle = np.flip(self.obstacle, 0)
        if self.args.verbose:
            plt.imshow(self.obstacle, origin='lower', cmap='Greys',  interpolation='nearest')
            plt.show()

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
        self.name = simulation_parameters["simulation_name"]
        self.id = simulation_parameters["simulation_id"]
        self.timesteps = simulation_parameters["time_steps"]
        self.n_x = simulation_parameters["lattice points x"]
        self.n_y = simulation_parameters["lattice points y"]
        self.reynolds = simulation_parameters["Reynolds number"]
        self.args = args

        # für die Simulation benötigten Felder initialisieren
        self.shape = (self.n_y, self.n_x)
        self.rho = np.empty(shape=self.shape)
        self.vel = np.empty(shape=(2, self.shape[0], self.shape[1]))
        self.f_in = np.empty(shape=(9, self.shape[0], self.shape[1]))
        self.f_eq = np.empty_like(self.f_in)
        self.f_out = np.empty_like(self.f_in)

        self.obstacles_definition(inputfile["obstacle parameters"])

        # als nächstes müssen die einstellungen für die boundary-condition
        # ausgelesen werden.
        # hier fehlt mir allerdings noch eine gute idee wie.
        # 1) passen die boundary conditions zusammen? (periodic)
        # 2) anpassen der obstacle-matrix für bounce-back?
        # 3) falls Zou-He: wie soll die kante behandelt werden?
        boundary_conditions = inputfile["boundary conditions"]

        # umgang mit den output-parametern:
        # also alles was mit speichern, plotten,
        # oder ausgabe während des programms zu tun hat

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
        pass

    def correct_distr_func(self):
        """
        Correct incoming components of distribution function at Zou-He Border

        Take the current distribution function :math:`f_i`
        and calculate corrected values
        for incoming components of the distribution function using Zou-He.
        """
        # TODO fulfill docstring
        pass

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
        self.f_in[(6, 3, 7), -1] = self.f_in[(6, 3, 7), -2]

    def do_simulation_step(self):
        self.calc_macroscopic()
        self.correct_macroscopic()
        self.calc_equilibrium()
        self.correct_distr_func()
        self.collision_step()
        self.bounce_back()
        self.stream_step()
        self.correct_outflow()

    def run_simulation(self):
        # initialisiere ein anfangs-geschwindigkeitsfeld.
        # das solle der anwender vorgeben!
        # wenn keines vorgegeben ist, sollte ein default generiert werden.
        # die 3 zeilen sollten also langfristig nicht hier,
        # sondern beim initialisieren stehen....
        self.vel[:] = 0
        self.vel[0, 0, :] = 0.04
        self.rho[:] = 1

        # errechne ein initiales f_in.
        # hier sollte langfristig auch die möglichkeit existieren
        # eine rehcnung mit einem existierenden f_in
        # fortzusetzen oder zu starten!
        # für den moment ist das aber erstmal okay so.
        self.calc_equilibrium()
        self.f_in = self.f_eq
        # die eigentliche Schleife der Simulation
        for step in range(self.timesteps):
            # print-befehl nur dazu da um zu sehen, dass schleife läuft...
            print(step)
            # nur zum debugging:
            # plt.imshow(np.linalg.norm(self.vel, axis=0), origin='lower')
            # plt.show()
            if step%100==0:
                plt.imshow((self.vel[0]*self.vel[0]*+self.vel[1]*self.vel[1]),vmin =0, vmax=1e-19, origin='lower')
                #plt.savefig("Bild_%i.jpeg" % step)
                plt.show()
            self.do_simulation_step()
