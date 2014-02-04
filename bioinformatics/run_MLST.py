#!/usr/bin/env python
#
# run_MLST.py
#
# This script takes a PubMLST (http://pubmlst.org) profile table, and the
# corresponding MLST sequences in FASTA format, and uses the protocol described
# in:
#
# - Larsen MV, Cosentino S, Rasmussen S, Friis C, Hasman H, et al. (2012)
#   Multilocus Sequence Typing of Total Genome Sequenced Bacteria. J Clin
#   Microbiol 50: 1355-1361. doi:10.1128/JCM.06094-11
#
# For any large-scale, persistent analysis, you're probably better off using
# BIGSdb:
#
# - Jolley KA, Maiden MCJ (2010) BIGSdb: Scalable analysis of bacterial
#   genome variation at the population level. BMC Bioinformatics 11: 595.
#   doi:10.1186/1471-2105-11-595.
#
# But if you have a one-shot, one-time use, this script may well be quicker to
# run than configuring the Apache and PostgreSQL infrastructure used by BIGSdb.
#
# DEPENDENCIES:
#
# USAGE: run_MLST.py [-h] [-o OUTFILENAME] [-i INDIRNAME] [-p PROFILE] [-v]
#                    [--blast_exe BLAST_EXE]
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -o OUTFILENAME, --outfile OUTFILENAME
#                         Output MLST classification table
#   -i INDIRNAME, --indirname INDIRNAME
#                         Directory containing MLST allele sequence files
#   -p PROFILE, --profile PROFILE
#                         Tab-separated plain text table describing MLST
#                         classification scheme.
#   -v, --verbose         Give verbose output
#   --blast_exe BLAST_EXE
#                         Path to BLASTN+ executable
#
# METHOD:
#
# The Larsen et al. (2012) method is followed:
#
# - The set of MLST sequences are used as a BLASTN query against the
#   input genome. Larsen et al. (2012) implies default parameters, so default
#   parameters are used.
# - The 'correct' allele is chosen on the basis of Length Score (LS), which
#   is calculated as LS = QL - HL + G where QL is the length of the MLST
#   query allele, HL is the length of BLASTN's reported HSP, and G is the
#   number of reported gaps in that HSP. A low LS score is better.
# - The combination of identified alleles is then used to determine the ST,
#   on the basis of an MLST profile table.
#
# The MLST profile table, and corresponding allele sequences for your
# organism can (probably) be obtained from PubMLST at http://pubmlst.org.
#
# (c) The James Hutton Institute 2014
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

#================
# IMPORTS

from argparse import ArgumentParser

import logging
import logging.handlers

import os
import sys
import traceback


#================
# FUNCTIONS

# Parse command-line
def parse_cmdline(args):
    """ Parse command-line arguments
    """
    parser = ArgumentParser(prog="run_MLST.py")
    parser.add_argument("-o", "--outfile", dest="outfilename",
                        action="store", default=None, type=str,
                        help="Output MLST classification table")
    parser.add_argument("-i", "--indirname", dest="indirname",
                        action="store", default=None, type=str,
                        help="Directory containing MLST allele sequence files")
    parser.add_argument("-p", "--profile", dest="profile",
                        action="store", default=None, type=str,
                        help="Tab-separated plain text table describing " +
                        "MLST classification scheme.")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true",
                        help="Give verbose output")
    parser.add_argument("--blast_exe", dest="blast_exe",
                        action="store", default="blastn",
                        help="Path to BLASTN+ executable")
    return parser.parse_args()


# Report last exception as string
def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))




if __name__ == '__main__':

    # Parse command-line
    # options are all options - no arguments
    options, args = parse_cmdline(sys.argv)

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    # err_handler points to sys.stderr
    # err_handler_file points to a logfile, if named
    logger = logging.getLogger('run_MLST.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if options.logfile is not None:
        try:
            logstream = open(options.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging" %
                         options.logfile)
            sys.exit(1)
    if options.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    logger.info('# calculate_ani.py logfile')
    logger.info('# Run: %s' % time.asctime())

    # Report arguments, if verbose
    logger.info(options)
    logger.info(args)

    # Have we got an input directory and profile? If not, exit.
    if options.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s" % options.indirname)
    if options.profile is None:
        logger.error("No MLST profile (exiting)")
        sys.exit(1)
    logger.info("MLST profile: %s" % options.indirname)
