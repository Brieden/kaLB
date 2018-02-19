# -*- coding: utf-8 -*-
"""
Unittests for core functions of the algorithm
"""
import unittest
import numpy as np
from src.d2q9_simulation import Simulation


class test_Simulation(unittest.TestCase):
    """
    Unittestclass for Simulation class
    """

    def setUp(self):
        """
        Create testdata

        Preparing testdata by randomly generating suitable arrays
        and mockup a dummy Simulation.
        """
        self.nx, self.ny = 500, 500
        self.test_sim = Simulation()
        self.test_sim.f_in = -np.random.uniform(1, 0, (9, self.nx, self.ny))
        self.test_sim.rho = -np.random.uniform(-1, 0, (self.nx, self.ny))
        self.test_sim.vel = np.random.normal(0, 0.02, (2, self.nx, self.ny))

    def test_calc_macroscopic(self):
        """
        Unittest for calc_macroscopic method

        Calculate control-values for density *rho*
        and velocity *vel* in an easy to understand way,
        to make sure NumPy methods used correctly.
        """

        # control calculation
        control_rho = np.zeros_like(self.test_sim.rho)
        control_vel = np.zeros_like(self.test_sim.vel)
        f = self.test_sim.f_in
        e = self.test_sim.e
        for xi in range(self.nx):
            for yi in range(self.ny):
                for i in range(9):
                    control_rho[xi, yi] += f[i, xi, yi]
                    control_vel[0, xi, yi] += f[i, xi, yi] * e[i, 0]
                    control_vel[1, xi, yi] += f[i, xi, yi] * e[i, 1]
        control_vel /= control_rho

        # do calculation in dummy Object
        self.test_sim.calc_macroscopic()

        # test
        self.assertTrue(np.allclose(control_rho, self.test_sim.rho))
        self.assertTrue(np.allclose(control_vel, self.test_sim.vel))

    def test_stream_step(self):
        """
        Unittest for stream_step method

        Calculate control-value for distribution function :math:`f_i`
        with "streamed" a.k.a. shifted values in an easy to understand way,
        to make sure NumPy *roll* methods are used correctly.
        """

        # control calculation
        control_f_in = np.empty_like(self.test_sim.f_in)
        f_out = self.test_sim.f_in
        for i in range(9):
            for xi in range(self.nx):
                x_next = xi + self.test_sim.e[i, 0]
                if x_next < 0:
                    x_next = self.nx - 1
                if x_next > self.nx - 1:
                    x_next = 0
                for yi in range(self.ny):
                    y_next = yi + self.test_sim.e[i, 1]
                    if y_next < 0:
                        y_next = self.ny - 1
                    if (y_next > (self.ny - 1)):
                        y_next = 0
                    control_f_in[i, x_next, y_next] = f_out[i, xi, yi]

        # do calculation in dummy Object
        self.test_sim.f_out = self.test_sim.f_in
        self.test_sim.stream_step()

        # test
        self.assertTrue(np.alltrue(control_f_in == self.test_sim.f_in))

    def test_calc_equilibrium(self):
        """
        Unittest for calc_equilibrium method

        Calculate control-value for equilibrium distribution function
        :math:`f_i^{eq}` in an easy to understand way,
        to make sure optimized calculation is done correctly.
        """

        # control calculation
        e = self.test_sim.e
        w = self.test_sim.w
        rho = self.test_sim.rho
        vel = self.test_sim.vel
        vel_sqared = vel[0] * vel[0] + vel[1] * vel[1]
        control_f_eq = np.empty_like(self.test_sim.f_in)
        e_dot_vel = np.dot(e, vel.transpose(1, 0, 2))
        for i in range(9):
            control_f_eq[i] = rho * w[i] * (
                1
                + (3 * e_dot_vel[i])
                + (4.5 * e_dot_vel[i] * e_dot_vel[i])
                - (1.5 * vel_sqared)
            )

        # do calculation in dummy Object
        self.test_sim.f_eq = np.empty_like(self.test_sim.f_in)
        self.test_sim.calc_equilibrium()

        # test
        self.assertTrue(np.allclose(control_f_eq, self.test_sim.f_eq))


if __name__ == '__main__':
    unittest.main()
