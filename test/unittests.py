# -*- coding: utf-8 -*-
"""
Unittests for core functions of the algorithm
"""
import unittest
import numpy as np
from src.d2q9_simulation import Simulation


class test_Something(unittest.TestCase):

    def setUp(self):
        self.nx, self.ny = 500, 500
        self.test_sim = Simulation()
        self.test_sim.f_in = -np.random.uniform(low=-1, high=0, size=(9, self.nx, self.ny))
        self.test_sim.rho = np.random.uniform(low=0.1, high=1.0, size=(self.nx, self.ny))
        self.test_sim.vel = np.random.normal(loc=0, scale=0.02, size=(2, self.nx, self.ny))

    def test_calc_macroscopic(self):
        control_rho = np.zeros_like(self.test_sim.rho)
        control_vel = np.zeros_like(self.test_sim.vel)
        f = self.test_sim.f_in
        for xi in range(self.nx):
            for yi in range(self.ny):
                for i in range(9):
                    control_rho[xi, yi] += f[i, xi, yi]
                    control_vel[0, xi, yi] += f[i, xi, yi] * self.test_sim.e[i,0]
                    control_vel[1, xi, yi] += f[i, xi, yi] * self.test_sim.e[i,1]
        control_vel /= control_rho
        self.test_sim.calc_macroscopic()
        self.assertTrue(np.allclose(control_rho, self.test_sim.rho))
        self.assertTrue(np.allclose(control_vel, self.test_sim.vel))

    def test_stream_step(self):
        control_f_in = np.empty_like(self.test_sim.f_in)
        f_out = self.test_sim.f_in
        for i in range(9):
            for xi in range(self.nx):
                x_next = xi + self.test_sim.e[i,0]
                if x_next < 0: x_next = self.nx - 1
                if x_next > self.nx - 1: x_next = 0
                for yi in range(self.ny):
                    y_next = yi + self.test_sim.e[i,1]
                    if y_next < 0: y_next = self.ny - 1
                    if (y_next > (self.ny - 1)):y_next = 0
                    control_f_in[i, x_next, y_next] = f_out[i, xi, yi]   
        self.test_sim.f_out = self.test_sim.f_in
        self.test_sim.stream_step()
        self.assertTrue(np.alltrue(control_f_in == self.test_sim.f_in))

    def test_calc_equilibrium(self):
        e = self.test_sim.e
        w = self.test_sim.w
        rho = self.test_sim.rho
        vel = self.test_sim.vel
        vel_sqared = vel[0] * vel[0] + vel[1] * vel[1]
        control_f_eq = np.empty_like(self.test_sim.f_in)
        for i in range(9):
            e_x_vel = e[i,0] * vel[0] + e[i,1] * vel[1]
            control_f_eq[i] = rho * w[i] * (
                1
                + (3 * e_x_vel)
                + (4.5 * e_x_vel * e_x_vel)
                - (1.5 * vel_sqared)
            )
        self.test_sim.f_eq = np.empty_like(self.test_sim.f_in)
        self.test_sim.calc_equilibrium()
        self.assertTrue(np.allclose(control_f_eq, self.test_sim.f_eq))


if __name__ == '__main__':
    unittest.main()
