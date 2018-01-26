#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
kaLB = kaum ausgereiftes Lattice Boltzmann
"""
import json
import argparse
from d2q9_simulation import Simulation


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
    return parser.parse_args()


def main():
    """
    Function
    """
    args = parse_arguments()
    input_file = args.input
    data = open_json(input_file)

    sim = Simulation()
    sim.initialize_from_json(data)
    sim.run_simulation()


if __name__ == '__main__':
    main()
