.. _link-to-manual:

Manual: calf flies
=====================

The command line parameter
--------------------------
+----+------------------------+-------------+-------------------+
|    | command line parameter | argument is | type of input     |
+====+========================+=============+===================+
|\-i |-input                  | required    | path to file      |
+----+------------------------+-------------+-------------------+
|\-o |-output                 | optional    | path to directory |
+----+------------------------+-------------+-------------------+
|\-np|-no_progessbar          | optional    | flag              |
+----+------------------------+-------------+-------------------+
|\-p |-performance_feedback   | optional    | flag              |
+----+------------------------+-------------+-------------------+
|\-so|-show_obstacle          | optional    | flag              |
+----+------------------------+-------------+-------------------+
|\-s |-snapshot               | optional    | path to file      |
+----+------------------------+-------------+-------------------+
|\-h |-help                   | optional    | flag              |
+----+------------------------+-------------+-------------------+








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












Understand the output
---------------------
As previously shown, three different outputs can be defined in the json file.::

        output/
        ├─ snapshot/
        │   ├─ snap_2000.npy   100 MB
        │   └─ snap_4000.npy   100 MB
        ├─ pic_1000.png         10 KB
        ├─ pic_2000.png         10 KB
        ├─ pic_3000.png         10 KB
        ├─ pic_4000.png         10 KB
        ├─ pic_5000.png         10 KB
        └─ raw_data.hdf5         7 GB


picture output configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This output is intended to check during the running simulaton whether the simulation runs well.
You can see at a glance whether the simulation is still running and that there are still no major numerical problems.

.. note::
    A high output rate of images slows down the simulation code very much and can fill the memory quickly.

raw data output configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This file format is very suitable for storing the data of the macroscopic value
and density for later processing in high storage frequency.
An example for the processing of the hdf5 file is the script :ref:`link-to-hdf5-to-mpeg.py`,
which creates a video from the data.
This script can be used as a basis to achieve a high-quality visualization of the simulation results.

The Hdf5 files are constructed as follows::

        kaLB_example_raw_data.hdf5/
          |
          └─raw data output/
                |
                ├─ velocity/
                │   ├─ velocity_02000
                │   ├─ velocity_04000
                │   ├─ ...
                │   ├─ velocity_88000
                │   └─ velocity_90000
                └─ density/
                    ├─ density_02000
                    ├─ density_04000
                    ├─ ...
                    ├─ density_88000
                    └─ density_90000


snapshot
^^^^^^^^
Snapshot in information technology is a full copy of a system or object.
It is a kind of backup to continue the simulation later from this point.
There are many scenarios in which you do not want to restart a simulation from the beginning.

Further advantages are that when starting with a former snapshot,
all parameters are read from the Json file and only the speed and density of the status are
taken over by the snapshot. This allows the user to continue to compute with other parameters.

To start a simulation with a snapshot of a previous simulation use the following call::

        $ python ./../src/kaLB.py --input kaLB_example.json -snapshot  output/snapshots/snap_02000.npy



Test: does the code do what it should?
--------------------------------------
kaLB provides unittests and a systemtest.

.. seealso::
    :ref:`link-to-testing`

