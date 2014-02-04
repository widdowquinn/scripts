#!/usr/bin/env python
#
# tabular_to_wikitable.py
#
# Script to convert plain text tab-separated tables to MediaWiki-compatible
# wikitable markup, preserving headers and titles where possible.
#
# USAGE: tabular_to_mediawiki.py [-h] [-o OUTFILENAME] [-i INFILENAME] [-v]
#                                [--header HEADER] [-t TITLE] [--skip SKIP]
#                                [-s] [-c]
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -o OUTFILENAME, --outfile OUTFILENAME
#                         Output MediaWiki format table
#   -i INFILENAME, --infile INFILENAME
#                         Input tab-separated plaintext table
#   -v, --verbose         Give verbose output
#   --header HEADER       If not passed a string of comma-separated headers,
#                         assumes first line of input tab-separated file is the
#                         header line
#   -t TITLE, --title TITLE
#                         Set the title of the resulting MediaWiki table
#   --skip SKIP           Number of lines to skip at the start of the input
#                         file, before reading the header
#   -s, --sortable        Make the MediaWiki table sortable
#   -c, --collapsible     Make the MediaWiki table collapsible
#
# (c) The James Hutton Institute 2014
# Authors: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD6 9LH,
# Scotland,
# UK
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
    parser = ArgumentParser(prog="tabular_to_mediawiki.py")
    parser.add_argument("-o", "--outfile", dest="outfilename",
                        action="store", default=None,
                        help="Output MediaWiki format table")
    parser.add_argument("-i", "--infile", dest="infilename",
                        action="store", default=None,
                        help="Input tab-separated plaintext table")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true",
                        help="Give verbose output")
    parser.add_argument("--header", dest="header",
                        action="store", default=None,
                        help="If not passed a string of comma-separated " +
                        "headers, assumes first line of input tab-" +
                        "separated file is the header line")
    parser.add_argument("-t", "--title", dest="title",
                        action="store", default=None,
                        help="Set the title of the resulting MediaWiki table")
    parser.add_argument("--skip", dest="skip",
                        action="store", default=0, type=int,
                        help="Number of lines to skip at the start of the " +
                        "input file, before reading the header")
    parser.add_argument("-s", "--sortable", dest="sortable",
                        action="store_true",
                        help="Make the MediaWiki table sortable")
    parser.add_argument("-c", "--collapsible", dest="collapsible",
                        action="store_true",
                        help="Make the MediaWiki table collapsible")
    return parser.parse_args()


# Report last exception as string
def last_exception():
    """ Returns last exception as a string
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Process the input stream
def process_stream(infh, outfh):
    """ Processes the input stream, assuming tab-separated plaintext input,
        with one row per line.
        If args.header is not set, assumes that the first line contains
        column headers.
        If args.title is set, adds a title to the output.
        If args.sortable is True, makes the table sortable.
        If args.collapsible is True, makes the table collapsible.
    """
    # Read in the input stream into a list of lines
    try:
        tbldata = list(infh.readlines())
    except:
        logger.error("Could not process input (exiting)")
        logger.error(last_exception())
        sys.exit(1)
    logger.info("Read %d lines from input" % len(tbldata))

    tbldata = tbldata[int(args.skip):]
    logger.info("Skipping %d lines from input" % args.skip)

    # How many columns are we expecting?
    cols = max([len(e.split('\t')) for e in tbldata])
    logger.info("Table appears to contain %d columns" % cols)

    # Prepare the table style
    tblclass = 'class="wikitable%s"'
    classes = []
    if args.sortable:
        classes.append("sortable")
    if args.collapsible:
        classes.append("mw-collapsible")
    # Apologies for the ternary operator
    tblclass = tblclass % (str(' ' if len(classes) else '') +
                           ' '.join(classes))
    initstr = '{|%s' % tblclass

    # Do we have a title?
    if not args.title:
        titlestr = "|+"
    else:
        try:
            titlestr = "|+ %s" % args.title
        except:
            logger.error("Could not process the passed title: %s (exiting)" %
                         args.title)
            logger.error(last_exception())
            sys.exit(1)

    # Do we have a passed header?
    if args.header is not None:
        try:
            headerstr = process_header(cols)
        except:
            logger.error("Could not process the passed header string: " +
                         "%s (exiting)" % args.header)
            logger.error(last_exception())
            sys.exit(1)
    else:  # We'll pop the first line of the row list as a header
        headerstr = '! ' + \
                    ' !! '.join([e.strip() for e in
                                 tbldata.pop(0).split('\t')])

    # Construct table body
    tblbody = '|-\n' + '\n|-\n'.join(['| ' + ' || '.join([e.strip() for e in
                                                          r.split('\t')])
                                      for r in tbldata])

    # Close table off
    tblclose = "|}\n"

    # Write out MediaWiki format table
    outfh.write('\n'.join([initstr, titlestr, headerstr, tblbody, tblclose]))


# Process a header that was passed from the command line
def process_header(cols):
    """ Splits the proposed header from the command line on commas, and checks
        for the appropriate number of columns (passed as cols). If this number
        is incorrect, a warning is given.
        If the number of columns in the header is too few, the cells are
        padded.
        If the number of columens is too great, the header is accepted as-is.
    """
    if not len(args.header):  # empty string passed
        headers = []
        headercount = 0
    else:
        headers = str(args.header).split(',')
        headercount = len(headers)
    if headercount > cols:
        logger.warning("Number of column headings (%d) " % headercount +
                       "greater than columns in data (%d)." % cols)
    elif headercount < cols:
        logger.warning("Number of column headings (%d) " % headercount +
                       "less than columns in data (%d): padding." % cols)
        while headercount < cols:
            headers.append("col%d" % (headercount + 1))
            headercount += 1
    headerstr = '! ' + ' !! '.join(headers)
    return headerstr


###
# SCRIPT

if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    logger = logging.getLogger('tabular_to_mediawiki.py')
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

    # Do we have an input file?  No? Then use stdin
    if args.infilename is None:
        infhandle = sys.stdin
        logger.info("Using stdin for input")
    else:
        logger.info("Using %s for input" % args.infilename)
        try:
            infhandle = open(args.infilename, 'rU')
        except:
            logger.error("Could not open input file: %s (exiting)" %
                         args.infilename)
            logger.error(''.join(
                traceback.format_exception(sys.last_type,
                                           sys.last_value,
                                           sys.last_traceback)))
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

    # Process input stream
    process_stream(infhandle, outfhandle)
