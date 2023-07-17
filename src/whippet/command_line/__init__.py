#
# Copyright (C) 2020 RFI
#
# Author: James Parkhurst
#
# This code is distributed under the GPLv3 license, a copy of
# which is included in the root directory of this package.
#
import logging
import time
from argparse import ArgumentParser
from typing import List
from whippet import simulate_and_reconstruct
from whippet import config


__all__ = ["main"]


# Get the logger
logger = logging.getLogger()


def get_description() -> str:
    """
    Get the program description

    """
    return "Simulate some images"


def get_parser(parser: ArgumentParser = None) -> ArgumentParser:
    """
    Get the parser

    """

    # Initialise the parser
    if parser is None:
        parser = ArgumentParser(prog="whippet", description=get_description())

    # Add arguments
    parser.add_argument(
        "-c",
        "--config",
        dest="config",
        type=str,
        default=None,
        help="The configuration",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Print out more information",
    )

    # Return the parser
    return parser


def main_impl(args):
    """
    Do the analysis

    """

    # Get the start time
    start_time = time.time()

    # Set the logger
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level, format="%(msg)s")

    # Simulate and reconstruct
    simulate_and_reconstruct(config.load(args.config))

    # Print output
    logger.info("Time taken: %.1f seconds" % (time.time() - start_time))


def main(args: List[str] = None):
    """
    Do the alignment

    """
    main_impl(get_parser().parse_args(args=args))
