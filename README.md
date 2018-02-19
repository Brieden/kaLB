# kaLB = kaum ausgereiftes Lattice Boltzmann

This repository contains the fluid dynamics simulation software **kaLB**.
It is part of the course *scientific software development*.

Creators:
* Sven Brieden
* Simon Schmitt

kaLB uses the *lattice boltzmann method* to compute a *2DQ9* fluid dynamics simulation.
The user is able to specify a simulation using a .json inputfile.
kaLB then processes this input, performs a simulation and returns a hdf5 file
containing the macroscopic values density and velocity during simulation steps.
In addition the user can use picture output and snapshots;
all of them at required output frequencies.

One can choose from four different boundary conditions
including *bounce-back* (a.k.a. no-slip) and *Zou-He*.
furthermore the user can specify arbitrary obstacles
and can even import pictures that kaLB will use as obstacles.

Besides this functionality kaLB provides an extra script
to create short videos from generated hdf5 files.

## Installation & Quick Start

These instructions will get you a copy of the project
up and running on your local machine:

##### Clone repository

	$ git clone git@git.physik.uni-marburg.de:Sven_Brieden/LBM_for_Fluid_Simulations.git

##### Install

Change directory to repo's root and run:

	$ python setup.py install

All dependencies should be installed automaticly

##### Run your first simulation

Change directory to examples folder and run:

	$ python ./../src/kaLB.py -i kaLB_example.json

wait until kaLB has done it's job!

##### Output & Video

Your output data will be located at *./output/*. Run:

	$ python ./../src/hdf5_to_mpeg.py -i ./output/kaLB_example_raw_data.hdf5
        
kaLB will create a short Video of your first Simulation!

This should get you going.
For a more detailed explanation on how to use kaLB refer to the documentation!

## Documentation

The documentation should provide you with all the necessary information you may need about kaLB.

To build the documentation navigate to the repos root directory and run:

    $ pip install .['doc']
    $ python setup.py build_sphinx

The documentation is now located in *./doc/build/html/index.html*
