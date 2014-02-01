#!/usr/bin/env python
#
# restrict_long_contigs.py
#
# USAGE: restrict_long_contigs.py [options] <input_directory> \
#                                       <output_directory>
#
# Options:
#   -h, --help            show this help message and exit
#   -l MINLEN, --minlen=MINLEN
#                         Minimum length of sequence
#   -s SUFFIX, --filesuffix=SUFFIX
#                         Suffix to indicate the file was processed
#   -v, --verbose         Give verbose output
#
# Non-PSL dependencies: Biopython (www.biopython.org)
#
# A short script that takes as input a directory containing (many) FASTA files
# describing biological sequences, and writes to a new, named directory
# multiple FASTA files containing the same sequences, but restricted only to
# those sequences whose length is greater than a passed value.
#
# Example usage: You have a directory with many sets of contigs from different
#  assemblies. This script will produce a new directory of the same data where
#  the contig lengths are restricted to being greater than a specified length.
#
# Copyright (C) 2013 The James Hutton Institute
# Author: Leighton Pritchard
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

from Bio import SeqIO
from optparse import OptionParser

import logging
import logging.handlers
import os
import re
import sys

###
# GLOBALS

# File extensions that indicate FASTA content
fasta_ext = ['.fa', '.fas', '.fasta']


###
# FUNCTIONS

# Parse cmd-line
def parse_cmdline(args):
    """ Parse command-line arguments. Note that the input and output
        directories are positional arguments
    """
    usage = "usage: %prog [options] <input_directory> <output_directory>"
    parser = OptionParser(usage)
    parser.add_option("-l", "--minlen", dest="minlen",
                      action="store", default=1000,
                      help="Minimum length of sequence")
    parser.add_option("-s", "--filesuffix", dest="suffix",
                      action="store", default="_restricted",
                      help="Suffix to indicate the file was processed")
    parser.add_option("-v", "--verbose", dest="verbose",
                      action="store_true", default=False,
                      help="Give verbose output")
    return parser.parse_args()


# Get list of FASTA files from a directory
def get_fasta_filenames(indir, extensions=fasta_ext):
    """ Identifies files in the passed directory whose extensions indicate
        that they may be FASTA files. Returns the path to the file,
        including the parent directory.
    """
    filelist = [f for f in os.listdir(indir) if
                os.path.splitext(f)[-1].lower() in extensions]
    logger.info("Identified %d FASTA files in %s:" % (len(filelist),
                                                      indir))
    if not len(filelist):       # We want there to be at least one file
        logger.error("No FASTA files found in %s" % indir)
        sys.exit(1)
    return filelist


# Restrict sequence length in a named FASTA file, writing it to
# the named location
def restrict_seq_length(infile, outfile, minlen):
    """ Takes an input FASTA file as infile, and writes out a corresponding
        file to outfile, where sequences shorter than minlen are not included
    """
    logger.info("Restricting lengths of %s to >=%d;" % (infile, minlen) +
                " writing to %s" % outfile)
    SeqIO.write([s for s in SeqIO.parse(infile, 'fasta')
                 if not len(s) < minlen],
                outfile, 'fasta')


# Process FASTA files in the directory
def process_files(indir, outdir, minlen, suffix):
    """ Takes an input directory that contains FASTA files, and writes
        to the output directory corresponding files (with the suffix appended)
        that contain only sequences of length greater than minlen.
    """
    for filename in get_fasta_filenames(indir):
        filestem, ext = os.path.splitext(filename)
        infilename = os.path.join(indir, filename)
        outfilename = os.path.join(outdir, ''.join([filestem, suffix, ext]))
        restrict_seq_length(infilename, outfilename, minlen)


###
# SCRIPT

if __name__ == '__main__':

    # Parse command-line
    # options are options, arguments are the .sff files
    options, args = parse_cmdline(sys.argv)

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    logger = logging.getLogger('restrict_long_contigs.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if options.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info(options)
    logger.info(args)

    # If there are not two positional arguments, throw an error
    if len(args) != 2:
        logger.error("Not enough arguments: script requires input and " +
                     "output directory")
        sys.exit(1)
    indir, outdir = tuple(args)

    # Make sure that the input directory exists
    if not os.path.isdir(indir):
        logger.error("Input directory %s does not exist" % indir)
        sys.exit(1)

    # If output directory does not exist, create it. If it does exist,
    # issue a warning that contents may be overwritten
    if os.path.isdir(outdir):
        logger.warning("Contents of %s may be overwritten" % outdir)
    else:
        logger.warning("Output directory %s does not exist: creating it" %
                       outdir)
        os.mkdir(outdir)

    # Check that the passed suffix is a valid string: escape dodgy characters
    #try:
    #    suffix = re.escape(options.suffix)
    #except:
    #    logger.error("Could not escape suffix string: %s" % options.suffix)
    #    sys.exit(1)

    # Make sure that the minimum length is an integer, and positive
    if not int(options.minlen) > 0:
        logger.error("Minimum length must be a positive integer, got %s" %
                     options.minlen)
        sys.exit(1)

    # Restrict sequence lengths
    process_files(indir, outdir, int(options.minlen), options.suffix)
