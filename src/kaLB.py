#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
kaLB is a fluid dynamic simulation implementation using the lattice boltzmann method.

Use this programm to start a simulation,
based on settings provided by the json inputfile.

This file is just the starter file and it is intended to work together with
d2q9_simulation.py and utilities.py.

kaLB = kaum ausgereiftes Lattice Boltzmann
"""
import json
import argparse
from src.d2q9_simulation import Simulation
from src import utilities


def parse_arguments():
    """
    Parse commandline arguments.

    :return: args
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input', required=True, type=str,
        help="Specify path to input file."
    )
    parser.add_argument(
        '-o', '--output', required=False, type=str,
        default='./output/',
        help="Specify path where to save the output."
    )
    parser.add_argument(
        '-np', '--no_progessbar', required=False, action='store_true',
        help="With this command, you can hide the progressbar."
    )
    parser.add_argument(
        '-p', '--performance_feedback', required=False, action='store_true',
        help="With this command, you can activate a performance-feedback at the end."
    )
    parser.add_argument(
        '-so', '--show_obstacle', required=False, action='store_true',
        help="With this command, you can visually check your obstacle before simulation starts."
    )
    parser.add_argument(
        '-s', '--snapshot', required=False, type=str,
        help="Specify path to the existing snapshot that you want to use as initial condition."
    )

    return parser.parse_args()


def open_json(filename):
    """
    Open *.json* input-file.

    :param filename: path to input-file in *.json* format.

    :return: input_data as dictionary.
    """

    try:
        input_data = json.load(open(filename))
    except IOError as error:
        print("ERROR: could not open input-file. " + str(error))
        quit()
    return input_data


def main():
    """
    Main function
    """

    args = parse_arguments()

    # is the inputfile the system test?
    if args.input == "system_test.json":
        do_systemtest = True
    else:
        do_systemtest = False

    # load json
    input_file = args.input
    json_file = open_json(input_file)

    # start simulation
    sim = Simulation(inputfile=json_file, args=args)
    sim.run_simulation()

    # system test analysis
    if do_systemtest:
        utilities.system_test(args)


if __name__ == '__main__':
    main()
