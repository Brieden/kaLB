.. kaLB documentation master file, created by
   sphinx-quickstart on Fri Jan 19 18:35:55 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:tocdepth: 3


kaLB = kaum ausgereiftes Lattice Boltzmann
==========================================

Welcome, this is the documentation for kaLB, a fluid dynamics simulation software.
It was developed as part of the course scientific software development
at the *Philipps University Marburg* in the winter term 2017/18.
It was created by:

* Sven Brieden
* Simon Schmitt

kaLB uses the lattice boltzmann method to compute a 2DQ9 fluid dynamics simulation.
The user is able to specify a simulation using a .json inputfile.
kaLB then processes this input, performs a simulation and returns a hdf5 file
containing the macroscopic values density and velocity during simulation steps.
In addition the user can use picture output and snapshots;
all of them at required output frequencies.

One can choose from four different boundary conditions
including bounce-back (a.k.a. no-slip) and Zou-He.
furthermore the user can specify arbitrary obstacles
and can even import pictures that kaLB will use as obstacles.

Besides this functionality kaLB provides an extra script
to create short videos from generated hdf5 files.

.. toctree::
   :numbered:
   :maxdepth: 2

   getting_started.rst
   manual.rst
   conventions_contribute.rst
   tests.rst
   functions_of_the_code.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
