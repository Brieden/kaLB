.. _link-to-manual:

Manual: to use a calf properly
==============================


Parameters of a simulation: The JSON file
-----------------------------------------
The Json file consists of four blocks, all of which must be present:

* simulation parameters
* boundary conditions
* obstacle parameters
* output configuration

Now we explain what these blocks are made of:

simulation parameters
^^^^^^^^^^^^^^^^^^^^^
Here, no parameter should be allowed.

* **simulation name and simulation id:**
    Here a clear identification of the simulation should be given.
    It can be used by the user for later assignment.

* **time steps and step offset:**
    time steps > 0

    step offset >= 0

    Time steps determine the number of steps the program should do.
    Step offset is an offset for the label in the output.
    It makes a lot of sense to continue a simulation with a snapshot.

* **lattice points x and y:**
    lattice points >= 10
    In order to get reasonably meaningful results, the grid should be at least 10x10 points.

boundary conditions
^^^^^^^^^^^^^^^^^^^

obstacle parameters
^^^^^^^^^^^^^^^^^^^

output configuration
^^^^^^^^^^^^^^^^^^^^















