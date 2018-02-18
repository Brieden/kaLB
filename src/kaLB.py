#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
kaLB = kaum ausgereiftes Lattice Boltzmann
"""
import json
import argparse
from src.d2q9_simulation import Simulation
from src import utilities
import os


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
        '-v', '--verbose', required=False, action='store_true',
        help="With this command, you can activate an issue of many hints."
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
    Function
    """
    args = parse_arguments()
    do_systemtest = False
    if args.input == "system_test.json":
        do_systemtest = True
    input_file = args.input
    json_file = open_json(input_file)

    sim = Simulation(inputfile=json_file, args=args)
    sim.run_simulation()
    if do_systemtest:
        utilities.system_test(args)


if __name__ == '__main__':
    main()
    print("\nEnde Gelände ")
