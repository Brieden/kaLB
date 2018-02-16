from setuptools import setup


setup(
    name='kaLB',
    version='0.1',
    author='Simon Schmitt, Sven Brieden',
    author_email='Briedens@uni-marburg.de, simon.schmitt@physik.uni-marburg.de',
    description='A python Lattice Boltzmann (LBM) simulation',
    keywords='lbm fluid cfd lattice boltzmann computational',
    packages=['src'],
    long_description='not sure what to put here', # TODO verweis auf README???
    install_requires=['numpy >= 1.13.3'],
    extras_require={
        'doc': ['Sphinx >= 1.6.5'],
        'visualization': ['matplotlib >= 1.3.1'],
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
