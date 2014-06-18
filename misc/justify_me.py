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
# TODO: Turn output into JSON/graphical output
#
# (c) L.Pritchard 2014

###
# IMPORTS

from argparse import ArgumentParser
from collections import defaultdict

import logging
import logging.handlers

import os
import re
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
    return [(os.path.split(texfile)[-1], scrape_time(texfile)) for 
            texfile in texfiles]
    

# Report time spent by lab book day
def report_by_day(times, outstream):
    """ Report time spent by lab book day
    """
    for filename, tlist in sorted(times):
        outstream.write("\n%s:\n" % filename)
        total = 0
        for topic, t in sorted(tlist):
            if t:
                outstream.write("\t%30s:\t%.2fh\n" % (topic, t/60.))
                total += t
        outstream.write("Total time recorded: %.2fh\n" % (total/60.))


# Report total time recorded in lab boo
def report_total_time(times, outstream):
    """ Report time recorded across all lab books
    """
    totals = defaultdict(int)
    days = 0
    outstream.write("\nTotal time recorded:\n")
    for filename, tlist in sorted(times):
        days += 1
        for topic, t in sorted(tlist):
            if t:
                totals[topic.upper()] += t
    total = sum(totals.values())
    for topic, t in sorted(totals.items()):
        outstream.write("\t%30s:\t%.2fh\t(%.2f%%)\n" %
                        (topic, t/60., 100.*t/total))
    outstream.write("Total time recorded: %.2fh\n" % (total/60.))
    outstream.write("Total time recorded per lab book: %.2fh\n" % 
                    (total/60./days))


# Takes an iterable of .tex files and processes \section and \subsection
# headers to scrape times spent under the header
def scrape_time(filename):
    """ Loops over a .tex file and scrapes the
        time spent (in format HHMM-HHMM) from each \section and \subsection
        header.
    """
    section_re = r"((?<=\\section\{).*(?=\})|(?<=\\subsection\{).*(?=\}))"
    logger.info("Scraping %s" % filename)
    with open(filename, 'rU') as fh:
        data = fh.read()
        # Get \section{} and \subsection{} elements
        matches = [m for m in re.findall(section_re, data) if
                   len(m.strip()) and ':' in m]
        # Process matches into subject, time values
        times = [process_match(m) for m in matches]
    return times


# Convert the section/subsection headers into a topic name and time spent
def process_match(match):
    """ Takes a string regex match for a (sub)section header of format:
        TOPIC: HHMM-HHMM; HHMM-HHMM...
        and returns a tuple of (TOPIC, TIME SPENT IN MINUTES)
    """
    time_re = "[0-9]{4}-[0-9]{4}"
    topic, times = match.split(':', 1)
    topic = topic.strip()
    times = re.findall(time_re, times)
    if not len(times):
        return (topic, 0)
    return (topic, calc_time(times))


# Convert string times HHMM-HHMM into time spent
def calc_time(times):
    """ Takes a list of times in HHMM-HHMM format, and returns the difference
        between the first and second times
    """
    cumt = 0
    for t in times:
        t1, t2 = t.split('-')
        tm = (int(t2[2:]) - int(t1[2:])) % 60
        th = 60 * ((int(t2[:2]) - int(t1[:2])) % 24)
        if int(t2[2:]) < int(t1[2:]):
            th -= 60
        cumt +=  tm + th
    return cumt

    

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
    times = process_labbooks()

    # Report time spent by day
    report_by_day(times, outfhandle)

    # Report total time spent
    report_total_time(times, outfhandle)
