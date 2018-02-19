Manual: How to use a calf properly
==================================

.. _link-to-inputfile:

Parameters of a simulation: The JSON file
-----------------------------------------
The Json file consists of four blocks, all of which must be present:

* simulation parameters
* boundary conditions
* obstacle parameters
* output configuration

Now we explain what these blocks contain:

.. seealso::
	Have a look inside *kaLB_example.json* to see a good inputfile!

simulation parameters
^^^^^^^^^^^^^^^^^^^^^

Here, all parameter are obligatory.

* **simulation name and simulation id:**
    Here a clear identification of the simulation should be given.
    It can be used by the user for later assignment.

* **time steps and step offset:**
    time steps > 0

	step offset >= 0

	Time steps determine the number of steps the program should do.
	Step offset is an offset for the label in the output.
	It is intended to be used, if you resume an already finished or interrupted simulation with a snapshot.
	Use the amount of steps in the previous iteration for this parameter. For a new simulation use 0.

* **lattice points x and y:**
    lattice points >= 10

    In order to get reasonably meaningful results, the grid should be at least 10x10 points.
	Keep in mind, that a computaion can get quite heavy and take a long time for large grids.
	Furthermore, big arrays can result in RAM shortage.
	**Be careful with large grids** (~	4000x4000).

* **tau:**
	tau > 1/2

	tau represents the viscosity of the fluid.
	Big values should result in a viscous fluid, small values vice versa.
	using *tau*=1 is fine and should be a good starting point if you're interested in a specific flow-scenario.
	**Be careful if tau is close to 1/2** code can easily get numerically instable.


boundary conditions
^^^^^^^^^^^^^^^^^^^

Here, all parameter are obligatory.

* **N, O, S, W:**
	specify dictionary containing a boundary-type for every cardinal direction border.
	
	valid dictionaries are:

	1. {"type": "bounce_back"}
		Introduce a thin 1-grid-point obstacle at the border.

	2. {"type": "outflow"}
		Inflow at these borders is repressed.

	3. {"type": "periodic"}
		Opposite direction has to be 'periodoc' as well!

	4. {"type": "zou-he", "v_x" : 0.04, "v_y" : 0}
		sets velocity components x and y at the border to these values in every iteration.

.. note::
	velocity components should not exceed 0.1 to get a numerically stable simulation.


obstacle parameters
^^^^^^^^^^^^^^^^^^^

obstacle parameters has to be a list.
An empty list is allowed, tho, if you don't want to have any obstacles.

3 types of obstacle types are supported:

1. **recktangle obstacle:** creating a rectangular obstacle. Make sure coordinates are meaningful.
	* **bottom_left:** bottom left corner coordinates of rectangular obstacle as list of 2
	* **top_right:** top right corner coordinates of rectangular obstacle as list of 2

2. **cylindrical obstacle:** create a circular obstacle.
	* **x-position:** x-coordinate of the center of the circle
	* **y-position:** y-coordinate of the center of the circle
	* **radius:** radius of the circle

3. **png import:**
	**file name:** name of the picture file.

	*PNG* should be in the same directory as *.json* inputfile.
	Other picture-types might work, but are not tested and therefore not officially supported!

.. seealso::
	Have a look at the documentation of their according helper functions in :ref:`link-to-utilils`.

.. note::
	you can always use the *-\\-show_obstacle* argument with your simulation
	to have a look at your obstacles before simulation starts.

output configuration
^^^^^^^^^^^^^^^^^^^^

output configuration has to be a dictionary,
but an empty dictionary is valid,
but not recommended, since there is no output in this case.

there are 3 types of output configuration:

1. **picture output configuration:** save velocity pictures at some timesteps during simulation
	* **file name:** name pre-fix for saved pictures
	* **file type:** file type for saved pictures
	* **output frequency:** number of iteration-steps between output

2. **raw data output configuration:** write density and velocity at some timesteps during simulation in hdf5 file
	* **file name:** name for saved hdf5 file
	* **output frequency:** number of iteration-steps between output

3. **snapshot:** save snapshots of distribution function at some timesteps during simulation
		* **output frequency:** number of iteration-steps between output
















