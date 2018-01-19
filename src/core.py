# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def initialize():
    return


def calc_macroscopic(f, e):
    rho = np.sum(f, axis=0)
    vel = np.dot(e.transpose(), f.transpose((1, 0, 2))) / rho
    return rho, vel


def calc_equilibrium(rho, vel, e, w):
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
    return f_in - (f_in - f_eq)/tau


def stream_step(f_out, e):
    f_in = np.empty_like(f_out)
    for i in range(9):
        f_in[i] = np.roll(f_out[i], shift=e[i], axis=(0, 1))
    return f_in


def correct_macroscopic(rho, vel):
    return rho, vel


def correct_distr_func(f_in, f_eq):
    return f_in


def bounce_back(f_in, f_out, obstacle, e_inverse):
    for i, j in enumerate(e_inverse):
        f_out[i, obstacle] = f_in[j, obstacle]
    return f_out


def correct_outflow(f_in):
    f_in[(6, 3, 7), -1] = f_in[(6, 3, 7), -2]
    return f_in


def save_state(step, f_in):
    return


def main():

    timesteps = 100
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

    for step in range(timesteps):   # main loop

        print(step)
        rho, vel = calc_macroscopic(f_in, e)
        rho, vel = correct_macroscopic(rho, vel)   # Zou-He
        f_eq = calc_equilibrium(rho, vel, e, w)
        f_in = correct_distr_func(f_in, f_eq)    # Zou-He
        f_out = collision_step(f_in, f_eq, tau)
        f_out = bounce_back(f_in, f_out, obstacle, e_inverse)
        f_in = stream_step(f_out, e)
        f_in = correct_outflow(f_in)
        save_state(step, f_in)

#        plt.imshow(np.sqrt(vel[0]**2+vel[1]**2).transpose(), origin='lower')
#        plt.show()
    return
