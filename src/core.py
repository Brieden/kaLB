# -*- coding: utf-8 -*-
r"""
kaLB = kaum ausgereiftes Lattice Boltzmann
-------------------------------------------
Corefunction to calculate fluid dynamics
"""
import numpy as np
import matplotlib.pyplot as plt
import time


def initialize():
    """
    Function
    """
    print("hello world")
    return


def calc_macroscopic(f, e):
    """
    Calculate macroscopic density and velocity.

    Take a distribution function :math:`f_i`
    and direction-vectors :math:`\\vec{e}_i`
    and compute related density :math:`\\rho` and velocity :math:`\\vec{u}`:

        .. math::
            \\rho(\\vec{x}) = \sum_{i} f_{i}(\\vec{x})
        .. math::
            \\vec{u}(\\vec{x}) = \\frac{1}{\\rho}\sum_{i} f_{i}\\vec{e}_{i}

    :param f: distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param e: directions: numpy-array with shape = (*i*, *2*),
        with *i* = 9 = Number of directions.

    :returns: **rho**: macroscopic density:
        numpy-array with shape = (*n_x*, *n_y*).

    :returns: **vel**: macroscopic velocity:
        numpy-array with shape = (*2*, *n_x*, *n_y*).
    """
    rho = np.sum(f, axis=0)
    vel = np.dot(e.transpose(), f.transpose((1, 0, 2))) / rho
    return rho, vel


def calc_equilibrium(rho, vel, e, w):
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

    :param rho: macroscopic density:
        numpy-array with shape = (*n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param vel: macroscopic velocity:
        numpy-array with shape = (*2*, *n_x*, *n_y*).
        with *i* = 9 = Number of directions.

    :param e: directions:
        numpy-array with shape = (*i*, *2*),
        with *i* = 9 = Number of directions.

    :param w: weights (according to directions):
        numpy-array with shape = (*i*),
        with *i* = 9 = Number of directions.

    :returns: **f_eq**: equilibrium distribution fuction:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.
    """
    vel_sqared = vel[0]*vel[0] + vel[1]*vel[1]
    e_n_x_vel = np.dot(e, vel.transpose(1, 0, 2))
    f_eq = np.empty(shape=(e.shape[0], vel.shape[1], vel.shape[2]))
    for i in range(9):
        f_eq[i] = rho * w[i] * (
            1
            + (3 * e_n_x_vel[i])
            + (4.5 * e_n_x_vel[i]**2)
            - (1.5 * vel_sqared)
            )
    return f_eq


def collision_step(f_in, f_eq, tau):
    """
    Perform collision-step and return the updated distribution function f_out.

    Use distribution function :math:`f_i`,
    equilibrium distribution function :math:`f_i^{eq}`
    and relaxation parameter :math:`\\tau`
    to calculate the updated distribution function :math:`f_i^*`:

        .. math::
            f_i^* = f_i - \\frac{1}{\\tau} (f_i - f_i^{eq})

    :param f_in: distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param f_eq: equilibrium distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param tau: relaxation parameter:
        float :math:`\\in` (.5, :math:`\\infty`].

    :returns: **f_out**: updated distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.
    """
    return f_in - (f_in - f_eq)/tau


def stream_step(f_out, e):
    """
    Perform streaming-step and return the shifted distribution function.

    Use distribution function :math:`f_i`
    and direction-vectors :math:`\\vec{e}_i`
    to calculate shifted distribution function :math:`f_i^*`:

        .. math::
            f_i^*(\\vec{x}) = f_i(\\vec{x} + \\vec{e}_i)

    :param f_out: distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param e: directions:
        numpy-array with shape = (*i*, *2*),
        with *i* = 9 = Number of directions.

    :returns: **f_in**: shifted distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.
    """
    f_in = np.empty_like(f_out)
    for i in range(9):
        f_in[i] = np.roll(f_out[i], shift=e[i], axis=(0, 1))
    return f_in


def correct_macroscopic(rho, vel, given_instream, direction_sets, direction):
    """
    Correct incoming density and velocity at a given border

    Take macroscopic density :math:`\\rho` and velocity :math:`\\vec{u}`
    and set velocity at the given border to given_instream
    and calculate according density at given border.

    :param rho: macroscopic density:
        numpy-array with shape = (*n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param vel: macroscopic velocity:
        numpy-array with shape = (*2*, *n_x*, *n_y*).
        with *i* = 9 = Number of directions.

    :param given_instream: known velocity at given border:
        numpy-array with shape = (*2*, *1*, *n_y*).

    :param direction_sets: dict to get sets of direction indices
        for corresponding direction specifier.

    :param direction: list of direction specifiers:
        List of direction specifier strings (e.g.: ["N","E"])

    :returns: **rho**: corrected macroscopic density:
        numpy-array with shape = (*n_x*, *n_y*).

    :returns: **vel**: corrected macroscopic velocity:
        numpy-array with shape = (*2*, *n_x*, *n_y*).
    """
    # TODO do what you should do!
    return rho, vel


def correct_distr_func(f_in, f_eq, direction_sets, direction):
    """
    Correct incoming components of distribution function at Zou-He Border

    Take the current distribution function :math:`f_i`
    and calculate corrected values
    for incoming components of the distribution function using Zou-He.

    :param f_in: distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param f_eq: equilibrium distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param direction_sets: dict to get sets of direction indices
        for corresponding direction specifier.

    :param direction: list of direction specifiers:
        List of direction specifier strings (e.g.: ["N","E"])

    :returns: **f_in**: distribution function that satisfies BC:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.
    """
    # TODO do what you should do!
    return f_in


def bounce_back(f_in, f_out, obstacle, e_inverse):
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
    into their inverse direction of the post-collision distribution function.

    :param f_in: distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param f_out: post-collision distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param obstacle: boolean obstacle matrix:
        numpy-array with shape = (*n_x*, *n_y*),
        with *True* for an obstacle.

    :param e_inverse: list of inverse direction indices (e[e_inverse] = -e):
        numpy-array with shape = (*i*),
        with *i* = 9 = Number of directions.

    :returns: **f_out**: distribution function that satisfies BC:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.
    """
    for i, j in enumerate(e_inverse):
        f_out[i, obstacle] = f_in[j, obstacle]
    return f_out


def correct_outflow(f_in, direction_sets, direction):
    """
    Correct distribution function

    Take post-streaming distribution function :math:`f_i`
    and replace incorrectly incoming components after streaming
    with second to last entries from that same distribution function.

    :param f_in: distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.

    :param direction_sets: dict to get sets of direction indices
        for corresponding direction specifier.

    :param direction: list of direction specifiers:
        List of direction specifier strings (e.g.: ["N","E"])

    :returns: **f_in**: corrected distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.
    """
    # TODO fulfill docstring
    f_in[(6, 3, 7), -1] = f_in[(6, 3, 7), -2]
    return f_in


def save_state(step, f_in, intervall=1):
    """
    Save distribution function for specified steps.

    :param step: current timestep:
        int

    :param intervall: intervalls in wich to save distribution functions:
        int

    :param f_in: current distribution function:
        numpy-array with shape = (*i*, *n_x*, *n_y*),
        with *i* = 9 = Number of directions.
    """
    # TODO do something useful and acually save it
    if step % intervall == 0:
        pass
    return


def main():
    """ Function"""

    timesteps = 30
    n_x = 200
    n_y = 100
    tau = 1

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
        4/9,
        1/9, 1/9, 1/9, 1/9,
        1/36, 1/36, 1/36, 1/36
    ])

    # TODO: muss noch outgesourced werden!
    obstacle = np.full(shape=(n_x, n_y), fill_value=False)
    vel = np.full(shape=(2, n_x, n_y), fill_value=0.04)
    vel[1] = 0
    vel[0, 1:] = 0

    f_in = calc_equilibrium(1, vel, e, w)
    v_max, t_0, t_n = 0, time.time(), time.time()
    for step in range(timesteps):   # main loop

        print(
            "Schritt: %i mit %.2f M su/s"
            % (step, (n_x*n_y*1e-6)/(time.time()-t_n))
        )
        t_n = time.time()
        rho, vel = calc_macroscopic(f_in, e)
        rho, vel = correct_macroscopic(rho, vel, 0, direction_sets, "W")
        f_eq = calc_equilibrium(rho, vel, e, w)
        f_in = correct_distr_func(f_in, f_eq, direction_sets, "W")
        f_out = collision_step(f_in, f_eq, tau)
        f_out = bounce_back(f_in, f_out, obstacle, e_inverse)
        f_in = stream_step(f_out, e)
        f_in = correct_outflow(f_in, direction_sets, "E")
        save_state(step, f_in)

        u = np.sqrt(vel[0]*vel[0]*+vel[1]*vel[1])
        v_max = max(v_max, u.max())
        plt.imshow(u.transpose(), vmin=0, vmax=v_max, origin='lower')
        plt.savefig("Bild_%i.jpeg" % step)
    print(
        "%i Schritte in %2.f sek mit %.2f M su/s"
        % (timesteps, time.time()-t_0,
           (n_x*n_y*timesteps*1e-6/(time.time()-t_0)))
         )
    return


if __name__ == '__main__':
    main()
