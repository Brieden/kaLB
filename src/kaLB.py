#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
kaLB = kaum ausgereiftes Lattice Boltzmann
"""
import json
import argparse
from d2q9_simulation import Simulation
from utilities import *


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
    return input_data


def parse_arguments():
    """
    Parse commandline arguments.

    :return: args
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input', required=False, type=str,
        default="./../examples/kaLB_example.json",
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
        '-t', '--test', required=False, action='store_true',
        help="With this command, you can start a systemtest"
    )

    return parser.parse_args()


def main():
    """
    Function
    """
    args = parse_arguments()
    if args.test:
        args.input = "./../examples/system_test.json"

    input_file = args.input
    data = open_json(input_file)

    sim = Simulation()
    sim.initialize_from_json(data, args)
    sim.run_simulation()
    if args.test:
        system_test(args)


if __name__ == '__main__':
    main()
    print("Ende Gelände ")
