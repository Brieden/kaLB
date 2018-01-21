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
    """ Function
    :Example:

    >>> import template
    >>> a = template.MainClass1()
    >>> a.function1(1,1,1)
    """
    return


def calc_macroscopic(f, e):
    """ Function"""
    rho = np.sum(f, axis=0)
    vel = np.dot(e.transpose(), f.transpose((1, 0, 2))) / rho
    return rho, vel


def calc_equilibrium(rho, vel, e, w):
    """ Function"""
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
    """ Function"""
    return f_in - (f_in - f_eq)/tau


def stream_step(f_out, e):
    """ Function"""
    f_in = np.empty_like(f_out)
    for i in range(9):
        f_in[i] = np.roll(f_out[i], shift=e[i], axis=(0, 1))
    return f_in


def correct_macroscopic(rho, vel):
    """ Function"""
    return rho, vel


def correct_distr_func(f_in, f_eq):
    """ Function"""
    return f_in


def bounce_back(f_in, f_out, obstacle, e_inverse):
    """ Function"""
    for i, j in enumerate(e_inverse):
        f_out[i, obstacle] = f_in[j, obstacle]
    return f_out


def correct_outflow(f_in):
    """ Function"""
    f_in[(6, 3, 7), -1] = f_in[(6, 3, 7), -2]
    return f_in


def save_state(step, f_in):
    """ Function"""
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
    v_max, t_0, t_n  = 0, time.time(), time.time()
    for step in range(timesteps):   # main loop

        print("Schritt: %i mit %.2f Mpoints/sek" %(step, (time.time()-t_n)/(n_x*n_y)*1e6))
        t_n = time.time()
        rho, vel = calc_macroscopic(f_in, e)
        rho, vel = correct_macroscopic(rho, vel)   # Zou-He
        f_eq = calc_equilibrium(rho, vel, e, w)
        f_in = correct_distr_func(f_in, f_eq)    # Zou-He
        f_out = collision_step(f_in, f_eq, tau)
        f_out = bounce_back(f_in, f_out, obstacle, e_inverse)
        f_in = stream_step(f_out, e)
        f_in = correct_outflow(f_in)
        save_state(step, f_in)

        u = np.sqrt(vel[0]*vel[0]*+vel[1]*vel[1])
        v_max = max(v_max, u.max())
        plt.imshow(u.transpose(), vmin=0, vmax = v_max, origin='lower')
        plt.savefig("Bild_%i.jpeg"%step)
    print("%i Schritte in %2.f sek mit %.2f Mpoints/sek" %(timesteps, time.time()-t_0,
                                (time.time()-t_0)/(n_x*n_y)*1e6/timesteps))
    return

if __name__ == '__main__':
    main()
