#!/usr/bin/env python

from setuptools import setup
import numpy as np

def do_setup():
    setup(name = 'kaLB',
          version = '0.1',
          description='A python Lattice Boltzmann (LBM) simulation',
          packages = ['kaLB'],
          author='Simon Schmitt, Sven Brieden',
          author_email='Briedens@uni-marburg.de, simon.schmitt@physik.uni-marburg.de',
          classifiers=[
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Physics'
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Developers',
          'Operating System :: POSIX :: Linux',
          ],
          keywords = 'lbm fluid cfd lattice boltzmann computational',
          install_requires = [
              'numpy >= 1.13.3',
          ],
          extras_require = {
              'doc': ['Sphinx >= 1.6.5'],
              'visualization': ['matplotlib >= 1.3.1'],
          }
         )

if __name__ == "__main__":
    do_setup()
