#####################
Testing
#####################

kaLB provides some unittests for the core methods of *Simulation* class,
as well as a systemtest to test the program as a whole.

Unit tests
==========
The unittests for some of kaLBs most important core methods are located in a separate package called 'test'.
There are three unittests for *calc_macroscopic*, *calc_equilibrium* and *stream_step*.

For every test a mockup Simulation instance is created
and initialized with randomly generated values in the expected range
for velocity, density and the distribution function.
The tests are performed on a 500x500 Grid and evaluated by NumPy's *numpy.allclose()* method.

To run the unittest you can either do it from the root of the package using setuptools::

	$ python setup.py test

Or using unittest::

	$ python -m unittest test/test_Simulation.py

System test
===========
The complete system is tested by the well-known scenario: flow in the pipe.
The analytical solution  of the velocity is a parabolic velocity profile.
The test can be started like a normal simulation via the json file "system_test.json".
On the basis of the name of the json files, kaLB recognizes the system test
and carries out the evaluation after the simulation:

>>> $ pwd
>>> /home/LBM_for_Fluid_Simulations/example
>>> $ python ../src/kaLB.py -i system_test.json

It is very easy to change the system_test.json so that the test fails.
The analysis is started default after 1000 time steps,
therefore one can change the system for example to 100x100 lattice point.
which fails because the fluid does not have enough time to build up its expected current profile.
All settings can be made via the system_test.json.
The analysis always starts after the simulation and evaluates the speed field from the hdf5 file
to the last saved time on the last row in the tube.
The scenario is constructed in such a way that a liquid flows in on the left side,
walls are defined at the top and bottom and the right side allows outflow.
