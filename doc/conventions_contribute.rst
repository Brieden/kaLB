Coding Conventions and Contributing
===================================

First off, thanks for taking the time to contribute!
But please read here which standards we use for this project.
We find coding standards extremely important -
one of the most important things one can implement on a collaborative project.

Code Style
----------
The Python code in this repo is meant to follow the PEP8 style
guide (a stylized version http://pep8.org). Except for the line length, we stick to the Pep8 awards.
Since we have set the line length to 100 in the  .../LBM_for_Fluid_Simulations/setup.cfg
you can test the code on Pep8 in all project directorys with the following command:

>>> pycodestyle .

Furthermore we use four spaces for a tab.

Before you push something please test:
--------------------------------------
* is every class / function / codeblock well documented
* system test
* unittest
* Pep 8 conformity
* simulation speed.

The simulation speed is currently around 3M su / s on an Intel® Core ™ i7-4500U CPU @ 1.80GHz × 4.
This value should not be far below.

The changelog file
------------------
The changelog format is based on
[Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
