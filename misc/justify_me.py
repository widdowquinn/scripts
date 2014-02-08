#!/usr/bin/env python
#
# justify_me.py
#
# A very short script to summarise time use from my personal lab book.
#
# My lab book is in LaTeX, and I record my work under \section{} and
# \subsection{} headers. The contents of these headers are in the format:
#
# SUBJECT: HHMM-HHMM; HHMM-HHMM;...
#
# The lab book source is stored in a hierarchical directory structure:
#
# YEAR/MM_month/YYYY-MM-DD/YYYY-MM-DD.tex
#
# This script takes a directory as input argument, searches all subdirectories
# for .tex source files, and scans them for \section{} and \subsection{}
# lines, parsing the content to calculate how much time was spent under each
# heading. It reports some summary statistics.
#
# (c) L.Pritchard 2014

###
# IMPORTS

from argparse import ArgumentParser

import logging
import logging.handlers

import os
import sys
import traceback

###
# FUNCTIONS


# Parse command-line
def parse_cmdline(args):
    """ Parse command-line arguments
    """
    parser = ArgumentParser(prog="justify_me.py")
    parser.add_argument("-o", "--outfile", dest="outfilename",
                        action="store", default=None,
                        help="Output file")
    parser.add_argument("-i", "--indir", dest="indirname",
                        action="store", default='.',
                        help="Input tab-separated plaintext table")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true",
                        help="Give verbose output")
    return parser.parse_args()


# Traverse subdirectories, collecting .tex files and processing the headers
def process_labbooks():
    """ Starting from the input directory, traverse all subdirectories,
        finding .tex files. Process each .tex file to find time spent under
        each heading, and collate.
    """
    # Traverse subdirectories and get list of lab book locations
    texfiles = []
    for root, dirs, files in os.walk(args.indirname):
        texfiles.extend([os.path.join(root, f) for f in files if 
                         os.path.splitext(f)[-1] == '.tex'])

    # Process each book, returning a dictionary of time spent, keyed by
    # section/subsection header
    times = scrape_times(texfiles)


###
# SCRIPT

if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    logger = logging.getLogger('justify_me.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info(args)

    # Make sure the input directory is a directory
    if not os.path.isdir(args.indirname):
        logger.error("Input path %s is not a directory (exiting)" %
                     args.indirname)
        sys.exit(1)

    # Do we have an output file?  No? Then use stdout
    if args.outfilename is None:
        outfhandle = sys.stdout
        logger.info("Using stdout for output")
    else:
        logger.info("Using %s for output" % args.outfilename)
        try:
            outfhandle = open(args.outfilename, 'w')
        except:
            logger.error("Could not open output file: %s (exiting)" %
                         args.outfilename)
            logger.error(''.join(
                traceback.format_exception(sys.last_type,
                                           sys.last_value,
                                           sys.last_traceback)))
            sys.exit(1)

    # Process lab books
    process_labbooks()
