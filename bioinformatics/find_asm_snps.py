#!/usr/bin/env python
#
# find_asm_snps.py
#
# This script uses MUMmer's nucmer and show-snps packages to identify SNPs
# between FASTA sequences in an input directory, in a pairwise manner, 
# producing tables and graphs of the differences between those sequences.
#
# This script is intended for bacterial assemblies, where we can assume 
# haploidy and don't have to take into account SNP frequency at a site. This
# means that if we trust the validity of the assembly at a position, we can
# assume the validity of the SNP.
#
# The script performs two tasks: 
# 1) A nucmer alignment of each pair of sequences
# 2) SNP identification with show-snps
# as outlined at http://mummer.sourceforge.net/manual/#snpdetection
#
# nucmer --prefix=ref_qry ref.fasta qry.fasta
# show-snps -Clr ref_qry.delta > ref_qry.snps
#
# The -C option restricts SNP identification to uniquely-aligning regions.
#
# DEPENDENCIES
# ============
#
# o MUMmer executables in the $PATH, or available on the command line
#
# USAGE
# =====
#
# 
# (c) The James Hutton Institute 2013
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
# The MIT License
#
# Copyright (c) 2010-2014 The James Hutton Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

#============
# IMPORTS

import collections
import logging
import logging.handlers
import math
import multiprocessing
import os
import shutil
import subprocess
import sys
import time
import traceback

from optparse import OptionParser


#=============
# GLOBALS
VERSION = "0.1"

#=============
# FUNCTIONS

# Parse command-line
def parse_cmdline(args):
    """ Parse command-line arguments for the script
    """
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outdir", dest="outdirname",
                      action="store", default=None,
                      help="Output directory")
    parser.add_option("-i", "--indir", dest="indirname",
                      action="store", default=None,
                      help="Input directory name")
    parser.add_option("-v", "--verbose", dest="verbose",
                      action="store_true", default=False,
                      help="Give verbose output")
    parser.add_option("--nucmer_exe", dest="nucmer_exe",
                      action="store", default="nucmer",
                      help="Path to NUCmer executable")    
    parser.add_option("--showsnps_exe", dest="showsnps_exe",
                      action="store", default="show-snps",
                      help="Path to show-snps executable")    
    parser.add_option("-l", "--logfile", dest="logfile",
                      action="store", default=None,
                      help="Logfile location")
    parser.add_option("-f", "--force", dest="force",
                      action="store_true", default=False,
                      help="Force file overwriting")
    parser.add_option("--noclobber", dest="noclobber",
                      action="store_true", default=False,
                      help="Don't nuke existing files")
    parser.add_option("--skip_nucmer", dest="skip_nucmer",
                      action="store_true", default=False,
                      help="Skip NUCmer runs, for testing " +\
                          "(e.g. if output already present)")
    parser.add_option("-g", "--graphics", dest="graphics",
                      action="store_true", default=False,
                      help="Generate heatmap of SNP differences")
    parser.add_option("--format", dest="gformat",
                      action="store", default="pdf",
                      help="Graphics output format [pdf|png|jpg]")
    parser.add_option("--version", dest="version",
                      action="store_true", default=False,
                      help="Print script version")
    return parser.parse_args()

# Report last exception as string
def last_exception():
    """ Returns last exception as a string, for use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value, 
                                              exc_traceback))

# Run a set of command lines using multiprocessing
def multiprocessing_run(cmdlines):
    """ Distributes the passed command-line jobs using multiprocessing.

        - cmdlines is an iterable of command line strings
    """
    logger.info("Running %d jobs with multiprocessing" % len(cmdlines))
    pool = multiprocessing.Pool()
    completed = []
    if options.verbose:
        callback_fn = logger_callback
    else:
        callback_fn = completed.append
    pool_outputs = [pool.apply_async(subprocess.call,
                                     (str(cline), ),
                                     {'stderr': subprocess.PIPE,
                                      'shell': sys.platform != "win32"},
                                     callback = callback_fn) \
                        for cline in cmdlines]
    pool.close()        # Run jobs
    pool.join()         # Collect output
    logger.info("Multiprocessing jobs completed:\n%s" % completed)

# Multiprocessing callback to logger
def logger_callback(val):
    """ Basic callback for multiprocessing just to log status of each job

        - val is an integer returned by multiprocessing, describing the run
             status
    """
    logger.info("Multiprocessing run completed with status: %s" % val)

# Run NUCmer pairwise on the input files, using multiprocessing
def pairwise_nucmer(filenames):
    """ Run NUCmer to generate pairwise alignment data for each of the 
        input FASTA files.

        - filenames is a list of input FASTA filenames, from which NUCmer 
              command lines are constructed

        We loop over all FASTA files in the input directory, generating NUCmer
        command lines for each pairwise comparison, and then pass those 
        command lines to be run using multiprocessing.
    """
    logger.info("Running pairwise NUCmer comparison to generate *.delta")
    cmdlines = []
    for idx, f1 in enumerate(filenames[:-1]):
        cmdlines.extend([make_nucmer_cmd(f1, f2) \
                            for f2 in filenames[idx+1:]])
    logger.info("NUCmer command lines:\n\t%s" % '\n\t'.join(cmdlines))
    if not options.skip_nucmer:
        multiprocessing_run(cmdlines)
    else:
        logger.warning("NUCmer run skipped!")

# Construct a command-line for NUCmer
def make_nucmer_cmd(f1, f2):
    """ Construct a command-line for NUCmer pairwise comparison, and return as 
        a string

        - f1, f2 are the locations of two input FASTA files for analysis

        We use the -mum option so that we consider matches that are unique in
        both the reference and the query.
    """
    prefix = os.path.join(options.outdirname, "%s_vs_%s" % \
                              (os.path.splitext(os.path.split(f1)[-1])[0],
                               os.path.splitext(os.path.split(f2)[-1])[0]))
    cmd = "%s -mum -p %s %s %s" % (options.nucmer_exe, prefix, f1, f2)
    return cmd

# Construct a command-line for show-snps
def make_showsnps_cmd(stem):
    """ Construct a command-line for show-snps SNP-finding, and return as
        a string

        - stem is the stem for the location of the .delta file for input

        We use the -C option to report only SNPs in uniquely-aligned regions.
    """
    outfilename = stem + '.snps'
    deltafilename = stem + '.delta'
    #prefix = os.path.join(options.outdirname, "%s_vs_%s" % \
    #                          (os.path.splitext(os.path.split(f1)[-1])[0],
    #                           os.path.splitext(os.path.split(f2)[-1])[0]))
    #cmd = "%s -mum -p %s %s %s" % (options.nucmer_exe, prefix, f1, f2)
    cmd = "%s -Clr %s > %s" % (options.showsnps_exe, deltafilename, 
                               outfilename)
    return cmd

# Get the list of FASTA files from the input directory
def get_fasta_files():
    """ Return a list of FASTA files in the input directory
    """
    infiles = get_input_files(options.indirname,# '.fna')
                              '.fasta', '.fas', '.fa', '.fna')
    logger.info("Input files:\n\t%s" % '\n\t'.join(infiles))
    return infiles

# Get list of input files in a directory
def get_input_files(dir, *ext):
    """ Returns a list of files in the input directory with the passed 
        extension

        - dir is the location of the directory containing the input files
        
        - *ext is a list of arguments describing permissible file extensions
    """
    filelist = [f for f in os.listdir(dir) \
                    if os.path.splitext(f)[-1] in ext]
    return [os.path.join(dir, f) for f in filelist]
    

# Create output directory if it doesn't exist
def make_outdir():
    """ Make the output directory, if required.

        This is a little involved.  If the output directory already exists,
        we take the safe option by default, and stop with an error.  We can, 
        however, choose to force the program to go on, in which case we can 
        either clobber or not the existing directory.  The options turn out
        as the following, if the directory exists:

        DEFAULT: stop
        FORCE: continue, and remove the existing output directory
        NOCLOBBER+FORCE: continue, but do not remove the existing output
    """
    if os.path.exists(options.outdirname):
        if not options.force:
            logger.error("Output directory %s would " % options.outdirname +\
                             "overwrite existing files (exiting)")
            sys.exit(1)
        else:
            logger.info("Removing directory %s and everything below it" % \
                            options.outdirname)
            if options.noclobber:
                logger.warning("NOCLOBBER: not actually deleting directory")
            else:
                shutil.rmtree(options.outdirname)
    logger.info("Creating directory %s" % options.outdirname)
    try:
        os.makedirs(options.outdirname)   # We make the directory recursively
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if options.noclobber and options.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)

# Find SNPs using nucmer
def find_nucmer_snps():
    """ Identify SNPs as described at 
        http://mummer.sourceforge.net/manual/#snpdetection

        This is a fairly blunt way to identify SNPs between two sequences that
        you already believe to be correct, though it is possible to filter
        on the basis of other confidence measures once the SNPs are called.

        This identifies SNPs only in uniquely-aligned regions. Since the
        alignment should be symmetrical, we ought only need to consider 
        ref v qry, and not also qry v ref. But we're currently doing both.

        The SNPs are left in the MUMmer show-snps output format.
    """
    logger.info("Finding SNPs with nucmer")
    infiles = get_fasta_files()
    pairwise_nucmer(infiles)

# Run NUCmer pairwise on input files, using multiprocessing
def pairwise_nucmer(filenames):
    """ Run NUCmer to generate pairwise alignment data for each input FASTA 
        file pair

        - filenames is a list of input FASTA filenames, from which NUCmer
          command lines may be constructed.

        We loop over the list of filenames, generating a NUCmer command line
        for each pairwise comparison, and associated show-snps command line.
        We also pass those command lines to multiprocessing, to be run.
    """
    logger.info("Compiling pairwise NUCmer jobs")
    cmdlines = []
    for idx, f1 in enumerate(filenames[:-1]):
        cmdlines.extend([make_nucmer_cmd(f1, f2) \
                            for f2 in filenames[idx+1:]])
    logger.info("NUCmer command lines:\n\t%s" % '\n\t'.join(cmdlines))
    logger.info("%d" % len(cmdlines))
    if not options.skip_nucmer:
        multiprocessing_run(cmdlines)
    else:
        logger.warning("NUCmer run skipped!")
    # We use the named output files from cmdlines to define inputs to 
    # show-snps, and pass these to multiprocessing
    logger.info("Compiling pairwise show-snps jobs")
    stems = [l.split()[3] for l in cmdlines] # Get output filestem
    snpcmdlines = [make_showsnps_cmd(stem) for stem in stems]
    logger.info("show-snps command lines:\n\t%s" % '\n\t'.join(snpcmdlines))
    logger.info("%d" % len(snpcmdlines))
    if not options.skip_nucmer:
        multiprocessing_run(snpcmdlines)
    else:
        logger.warning("show-snps run skipped!")


#=============
# SCRIPT

if __name__ == '__main__':

    # Parse command-line
    # options are all options - no arguments
    options, args = parse_cmdline(sys.argv)

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    # err_handler points to sys.stderr
    # err_handler_file points to a logfile, if named
    logger = logging.getLogger('find_asm_snps.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = \
                  logging.Formatter('%(levelname)s: %(message)s')
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
    logger.info('# find_asm_snps.py logfile')
    logger.info('# Version: %s' % VERSION)
    logger.info('# Run: %s' % time.asctime())
 
    # Report arguments, if verbose
    logger.info(options)
    logger.info(args)

    # If we're asked for a version number, give it to stdout and exit normally
    if options.version:
        print >> sys.stdout, str(VERSION)
        sys.exit(0)

    # Have we got an input and output directory? If not, exit.
    if options.indirname is None:
        logger.error("No input directory name (exiting)")
        logger.error("Use find_asm_snps.py -h for help text")
        sys.exit(1)
    logger.info("Input directory: %s" % options.indirname)
    if options.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    make_outdir()
    logger.info("Output directory: %s" % options.outdirname)

    # Find SNPs
    find_nucmer_snps()
