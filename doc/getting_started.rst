Getting Started
===============

Here we show you how to install and start kaLB with the provided examples.

Virtual Environment
-------------------

We recommend to create a python virtualenvironment
and perform all simulations and development from within this environment.

with *virtualenv* just run::

	$ virtualenv <NAME-OF-YOUR-ENV>
	$ source <NAME-OF-YOUR-ENV>/bin/activate

Or if you're using conda::

	$ conda create --name <NAME-OF-YOUR-ENV> python=3
	$ conda activate <NAME-OF-YOUR-ENV>


Installation
------------

To install kaLB you can either use pip or use setuptools directly.
Navigate to kaLB's root directory and use:

Via setuptools::

	$ python setup.py install

Or with pip::

	$ pip install .

All dependencies should be installed automagically!

If you want to install kaLB in development mode use::

	$ python setup.py develop

or respectively::

	$ pip install -e .


Documentation
-------------

To generate the documentation you need to have sphinx installed.
Thankfully you can just use setuptools do install required version if you don't have it already::

	$ pip install .['doc']

To now build the documentation we recommend using setuptools again::

	$ python setup.py build_sphinx


Run a Simulation
----------------

To perform a Simulation you need to have an input *.json*-file
that holds all information about the simulation you want to perform.
For more information about inputfiles see :ref:`link-to-manual`.

Execute *kaLB.py* from the directory your inputfile is located.

kaLB comes with some predefined input files for different scenarios.

Navigate to *examples* folder and start kaLB with::

	$ python ./../src/kaLB.py --input kaLB_example.json --performance_feedback --show_obstacle

kaLB should display the matrix that represents the obstacle.
As soon as you dismiss it, the simulation will start.

During simulation an *output* directory will be created
to store pictures, snaphots and the hdf5 file.

Once the simulation has finished,
you will be presented with a performance feedback
that gives you an idea how many grid-points where calculated per second.

Congratulation! You have finished your first simulation.
Maybe you want to try another simulation right away?

Use "**Lid_Driven_Cavity.json**" for your next "**--input**" value!


Create a Video
--------------

After a Simulation you might want to see how your simulation hast developed over time.
Of course you could just look at the output-picutures,
but thats not that great, isn't it?

In this case you can use *hdf5_to_mpeg.py* script
to create a video from the generated hdf5 file.

After you have simulated the *kaLB_example.json* you could use::

	$ python ./../src/hdf5_to_mpeg.py -i ./output/kaLB_example_raw_data.hdf5

to create a video!

.. note::
	To make sure you're having enough data-points for a nice video, the *output frequency* should be adequate.


Tests
-----

kaLB provides unittests and a systemtest.

.. seealso::
	To learn more about testing and how to use it, see :ref:`link-to-testing`.

