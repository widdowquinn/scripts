#!/usr/bin/env python
#
# run_signalp.py [euk|gram+|gram-] <FASTAfile> [-o|--outfilename <output file>]
#
# Takes a FASTA format file containing protein sequences as input, and runs
# a local copy of SIGNALP on the contents, collecting the output
# generated with the -short option
# This script splits the input file into a number of smaller FASTA files
# suitable for distributed processing using, e.g. multiprocessing, with
# a maximum number of sequences per file of 1000 to avoid falling foul
# of SIGNALP's undocumented upper limit.  SIGNALP is run independently on
# each split file, and the results concatenated into the output file.
# If no output file is specified, then the output file shares a stem with
# the input file.
# The organism type is specified at the command line in the same way as for
# SIGNALP
#
# (c) The Scottish Crop Research Institute 2010
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


from Bio import SeqIO
from optparse import OptionParser
import multiprocessing                  # Needs Python 2.6
import os
import sys
import subprocess


###
# FUNCTIONS

# Parse the commandline
def parse_cmdline(args):
    """ Parse command-line args
    """
    usage = "usage: %prog [options] <organism_type> <infile>"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outfilename", dest="outfilename",
                      action="store", default=None,
                      help="Location to write collated signalp output")
    parser.add_option("--nomulti", dest="nomulti",
                      action="store_true", default=False,
                      help="Don't use multiprocessing")
    parser.add_option("-v", "--verbose", dest="verbose",
                      action="store_true", default=False,
                      help="Give verbose output")
    return parser.parse_args()


# Split a list into chunks of size n
def chunklist(l, n):
    """ Generator to split the passed list into chunks of size n
    """
    for i in xrange(0, len(l), n):
        seqs = [s[:60] for s in l[i:i+n]]
        yield seqs


# Split sequences into approximately even-sized subgroups.  There's no
# point submitting sequences longer than about 60aa to SIGNALP, so we
# also trim them, but do not modify the sequence ID or description.
def split_sequences(seqlist, nomulti):
    """ Loop over the passed list of sequences, and divide them into an
        approximately evenly-sized number of subsets of sequences with
        sequence length 60.
        The number of subsets is determined by the following constraints where
        n is the total number of sequences, c is the number of available
        processes, and s is the subset size:
        1. s <= 1000
        2. if n <= c then s = 1
        3. if c < n < 1000c then s = n/c
        4. if n >= 1000c then s = 1000
    """
    if nomulti:
        cpus = 1
    else:
        cpus = multiprocessing.cpu_count()    # Count of available CPUs

    # What size of subset do we need?
    if len(seqlist) <= cpus:
        subset_size = 1
    elif len(seqlist) < 1000 * cpus:
        subset_size = len(seqlist)/cpus
    else:
        subset_size = 1000

    # Split the sequences into subsets of the required size, without loss
    return list(chunklist(seqlist, subset_size))


# Run the passed list of command lines using multiprocessing
def mp_run(clines, poolsize):
    """ We create a new multiprocessing Pool to handle the several command
        lines we're passed
    """
    # Create processing pool and run jobs
    pool = multiprocessing.Pool(processes=poolsize)
    complete = []
    pool_output = [pool.apply_async(subprocess.call,
                                    (str(cline), ),
                                    {'stderr': subprocess.PIPE,
                                     'shell': sys.platform != "win32"})
                   for cline in clines]
    pool.close()           # Run jobs
    pool.join()


# Create a set of command-lines to run SIGNALP on the passed list of input
# files, and to return the list of output filenames
def run_signalp(filenames, organism, nomulti):
    """ Generates a list of SIGNALP command lines with the passed organism
        option and either runs them all on one processor, or distributes
        jobs, depending on the nomulti option.
    """
    clines = []                 # Holds command-lines
    outfilenames = []           # Holds output filenames

    # Generate command-lines
    for f in filenames:
        outfilename = os.path.splitext(f)[0] + '.sigp'
        outfilenames.append(outfilename)
        cline = 'signalp -t %s -short %s > %s ' % (organism, f,
                                                   outfilename)
        clines.append(cline)

    # Pass the command-lines on to be run
    if nomulti:
        for cline in clines:
            subprocess.call(cline, shell=sys.platform != 'win32')
    else:
        mp_run(clines, multiprocessing.cpu_count())

    # Clean up input files
    for f in filenames:
        os.unlink(f)

    # Return the list of output files
    return outfilenames

###
# SCRIPT

if __name__ == '__main__':

    # Parse cmdline
    options, args = parse_cmdline(sys.argv)
    organism = args[0]
    infilename = args[1]

    # If our desired tmp directory doesn't exist, create it
    if not os.path.isdir('/var/tmp/signalp'):
        os.makedirs('/var/tmp/signalp')

    # Load sequences - we turn them into a list here as we have to do it
    # sometime - and split them into a suitable set size
    subsets = split_sequences(list(SeqIO.parse(infilename, 'fasta')),
                              options.nomulti)

    # Write each sequence subset out to its own new FASTA file and keep
    # a record of them
    filestem = os.path.splitext(os.path.split(infilename)[-1])[0]
    subset_filenames = []
    for idx, ss in enumerate(subsets):
        new_filename = os.path.join('/var/tmp/signalp',
                                    '%s.%d.fas' % (filestem, idx))
        subset_filenames.append(new_filename)
        SeqIO.write(ss, new_filename, 'fasta')

    # Run SIGNALP on each of the sequence subsets
    output_filenames = run_signalp(subset_filenames, organism,
                                   options.nomulti)

    # Collate the output filenames together
    outfilename = os.path.splitext(infilename)[0]+'.sigp' if not\
        options.outfilename else options.outfilename
    # Seed output with the first output file
    cline = 'cat %s > %s' % (output_filenames[0], outfilename)
    if options.verbose:
        print cline
    subprocess.call(cline, shell=sys.platform != 'win32')
    # Process the rest of the files to remove headers
    for ofile in output_filenames[1:]:
        cline = "cat %s | awk 'FNR>2{print}' >> %s" % (ofile, outfilename)
        if options.verbose:
            print cline
        subprocess.call(cline, shell=sys.platform != 'win32')

    # Clean up /var/tmp/signalp
    for ofile in output_filenames:
        os.unlink(ofile)
