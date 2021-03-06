from setuptools import setup

setup(
    name='kaLB',
    version='1.0',
    author='Simon Schmitt, Sven Brieden',
    author_email='Briedens@uni-marburg.de, simon.schmitt@physik.uni-marburg.de',
    description='A python Lattice Boltzmann (LBM) simulation',
    keywords='lbm fluid cfd lattice boltzmann computational',
    packages=['src'],
    test_suite='test',
    long_description='not sure what to put here',
    install_requires=[
        'numpy >= 1.13.3',
        'h5py >= 2.7.1',
        'matplotlib >= 2.1.1'
    ],
    extras_require={
        'doc': ['Sphinx >= 1.6.5'],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Physics'
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: POSIX :: Linux',
    ],
)
