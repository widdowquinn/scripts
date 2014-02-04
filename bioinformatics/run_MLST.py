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
# - Biopython
# - Python 2.6+ (for multiprocessing)
# - Pandas
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
# - The 'correct' allele is chosen on the basis of: (1) Length Score (LS),
#   which is calculated as LS = QL - HL + G where QL is the length of the MLST
#   query allele, HL is the length of BLASTN's reported HSP, and G is the
#   number of reported gaps in that HSP - a low LS score is better; and (2)
#   percentage identity of the match. This is the method proposed in Larsen
#   et al. (2012)
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

import csv
import multiprocessing
import os
import shutil
import subprocess
import sys
import time
import traceback

from Bio import SeqIO

import pandas as pd

#================
# FUNCTIONS

# Parse command-line
def parse_cmdline(args):
    """ Parse command-line arguments
    """
    parser = ArgumentParser(prog="run_MLST.py")
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default=None, type=str,
                        help="Output MLST classification data directory")
    parser.add_argument("-i", "--indirname", dest="indirname",
                        action="store", default=None, type=str,
                        help="Directory containing MLST allele sequence files")
    parser.add_argument("-g", "--genomedir", dest="genomedir",
                        action="store", default=None, type=str,
                        help="Directory containing genome sequence files")
    parser.add_argument("-p", "--profile", dest="profile",
                        action="store", default=None, type=str,
                        help="Tab-separated plain text table describing " +
                        "MLST classification scheme.")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None, type=str,
                        help="Logfile location")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true",
                        help="Give verbose output")
    parser.add_argument("-f", "--force", dest="force",
                        action="store_true",
                        help="Force overwriting of output directory")
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


# Process the allele sequences from the input directory
def load_alleles():
    """ Identify the FASTA files in the input directory, and the genes/
        allele numbers involved.

        Return a dictionary of (filename, number of alleles) tuples, 
        keyed by gene name.
    """
    # We keep all our data in a dictionary:
    allele_dict = {}
    # Get the list of FASTA files to process
    logger.info("Processing input directory %s" % args.indirname)
    try:
        fastalist = get_input_files(args.indirname, 
                                    '.fasta', '.fas', '.fna', '.fa')
    except:
        logger.error("Could not identify FASTA files in directory %s " %
                     args.indirname + "(exiting)")
        logger.error(last_exception())
        sys.exit(1)
    # Loop over the FASTA sequences, and identify the genes and alleles.
    # We assume that the standard PubMLST identifier format is being used:
    # >gene_N
    # where gene is the gene name, and N is the allele number
    for filename in fastalist:
        try:
            alleles = list(SeqIO.parse(filename, 'fasta'))
            logger.info("Read %d alleles from %s" % (len(alleles), filename))
            gene_name = alleles[0].id.split('_')[0]
            allele_dict[gene_name] = (filename, len(alleles))
        except:
            logger.error("Could not process FASTA file %s (exiting)" % 
                         filename)
            logger.error(last_exception())
            sys.exit(1)
    return allele_dict


# Process genome sequences from the input directory
def load_genomes():
    """ Identify the genome FASTA sequences, and associate with a filename
    """
    # Keep data in a dictionary
    genomes = {}
    # Get the list of FASTA files to process
    try:
        fastalist = get_input_files(args.genomedir, 
                                    '.fasta', '.fas', '.fna', '.fa')
    except:
        logger.error("Could not identify FASTA files in directory %s " %
                     args.indirname + "(exiting)")
        logger.error(last_exception())
        sys.exit(1)
    # Loop over genome sequences and associate an identifier with a filename
    for filename in fastalist:
        try:
            genome = list(SeqIO.parse(filename, 'fasta'))
            logger.info("Read %d sequences from %s" % (len(genome), filename))
            isolate = genome[0].id.replace("|", "_")
            genomes[isolate] = filename
        except:
            logger.error("Could not process genome file %s (exiting)" % 
                         filename)
            logger.error(last_exception())
            sys.exit(1)
    return genomes


# Run pairwise BLASTN in all combinations of allele sequence and genome
# sequence input
def run_blast(alleles, genomes):
    """ We pass a series of BLAST-2-sequences commands to multiprocessing,
        writing to the prescribed output directory.
    """
    logger.info("Running BLASTN-2-sequences to identify alleles in each " +
                "genome")
    cmdlines = []
    blastoutfiles = {}
    for allele, data in alleles.items():
        for isolate, filename in genomes.items():
            outfilename = os.path.join(args.outdirname, "%s_vs_%s.tab" %
                                       (allele, isolate))
            cmdlines.append(make_blast_cmd(allele, isolate, data[0], 
                                           filename, outfilename))
            blastoutfiles[(isolate, allele)] = outfilename
    logger.info("Generated %d command-lines" % len(cmdlines))
    multiprocessing_run(cmdlines)
    return blastoutfiles
            

# Assign alleles to isolates
def assign_alleles(filedict, df):
    """ Takes a dataframe (df) and a dictionary of filenames keyed by (row,
        column) of that dataframe, and reads the contents of the file -
        a BLASTN comparison of alleles to the isolate, assigning the 'best'
        allele number to each genome on the basis of the LS score and
        percentage identity, as in Larsen et al. (2012).
    """
    for k, v in filedict.items():
        best_allele = find_best_allele(v)
        df.loc[k[0], k[1]] = best_allele
    return df


# Find the 'best' allele from a BLASTN output file
def find_best_allele(filename):
    """ Identifies the 'best' allele according to Larsen et al. (2012) from
        the BLASTN output file
    """
    with open(filename, 'r') as fh:
        csvreader = csv.DictReader(fh, delimiter='\t',
                                   fieldnames=['qseqid', 'sseqid', 'pident',
                                               'qlen', 'length', 'gaps', 
                                               'mismatch'])
        # Hold 'best' allele by LS and percentage identity
        best_ls = (None, None, None)
        best_pident = (None, None, None)
        # Identify 'best' allele
        for row in csvreader:
            allele = row['qseqid'].split('_')[-1]
            ls = int(row['qlen']) - int(row['length']) + int(row['gaps'])
            pident = float(row['pident'])
            if pident > best_pident[1]:
                best_pident = (ls, pident, allele)
            if best_ls is None:
                best_ls = (ls, pident, allele)
            elif ls <= best_ls and pident > best_ls[1]:
                best_ls = (ls, pident, allele)                
        if best_ls[-1] != best_pident[-1]:
            logger.warning("%s: Same allele not identified by LS and %%ID" %
                           filename)
            logger.warning("Best LS: %s" % str(best_ls))
            logger.warning("Best %%ID: %s" % str(best_pident))
        return best_ls[-1]
        
                

# Construct a BLASTN-2-sequences command line (default settings)
def make_blast_cmd(qid, sid, qfilename, sfilename, outfilename):
    """ Construct a default BLASTN-2-sequences command line from the passed
        query identifier and filename (qid, qfilename), and subject 
        identifier and filename (sid, sfilename), placing the result in 
        the file outfilename.
    """
    cmdline = "%s -query %s -subject %s -out %s" %\
              (args.blast_exe, qfilename, sfilename, outfilename)
    # Custom format for output
    cmdline += " -outfmt '6 qseqid sseqid pident qlen length gaps mismatch'"
    return cmdline


# Run a set of command lines using multiprocessing
def multiprocessing_run(cmdlines):
    """ Distributes the passed command-line jobs using multiprocessing.

        - cmdlines is an iterable of command line strings
    """
    logger.info("Running %d jobs with multiprocessing" % len(cmdlines))
    pool = multiprocessing.Pool()
    completed = []
    if args.verbose:
        callback_fn = logger_callback
    else:
        callback_fn = completed.append
    pool_outputs = [pool.apply_async(subprocess.call,
                                     (str(cline), ),
                                     {'stderr': subprocess.PIPE,
                                      'shell': sys.platform != "win32"},
                                     callback=callback_fn)
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


# Get a list of files with a particular (set of) extension(s) from the passed
# directory
def get_input_files(dir, *ext):
    """ Returns a list of files in the input directory with the passed
        extension

        - dir is the location of the directory containing the input files

        - *ext is a list of arguments describing permissible file extensions
    """
    filelist = [f for f in os.listdir(dir)
                if os.path.splitext(f)[-1] in ext]
    return [os.path.join(dir, f) for f in filelist]


# Create output directory if it doesn't exist
def make_outdir():
    """ Make the output directory, if required.

        This is a little involved.  If the output directory already exists,
        we take the safe option by default, and stop with an error.  We can,
        however, choose to force the program to go on, in which case we 
        clobber the existing directory.

        DEFAULT: stop
        FORCE: continue, and remove the existing output directory
    """
    if os.path.exists(args.outdirname):
        if not args.force:
            logger.error("Output directory %s would " % args.outdirname +
                         "overwrite existing files (exiting)")
            sys.exit(1)
        else:
            logger.info("Removing directory %s and everything below it" %
                        args.outdirname)
            shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s" % args.outdirname)
    try:
        os.makedirs(args.outdirname)   # We make the directory recursively
    except OSError:
        logger.error(last_exception)
        sys.exit(1)


if __name__ == '__main__':

    # Parse command-line
    # options are all options - no arguments
    args = parse_cmdline(sys.argv)

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    # err_handler points to sys.stderr
    # err_handler_file points to a logfile, if named
    logger = logging.getLogger('run_MLST.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging" %
                         args.logfile)
            sys.exit(1)
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    logger.info('# run_MLST.py logfile')
    logger.info('# Run: %s' % time.asctime())

    # Report arguments, if verbose
    logger.info(args)

    # Have we got an input/output directory and profile? If not, exit.
    if args.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s" % args.indirname)

    if args.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    logger.info("Output directory: %s" % args.outdirname)

    if args.profile is None:
        logger.error("No MLST profile (exiting)")
        sys.exit(1)
    logger.info("MLST profile: %s" % args.indirname)

    # Process the allele sequences
    alleles = load_alleles()
    logger.info("Found sequences for %d MLST genes:" % len(alleles))
    for k, v in alleles.items():
        logger.info("%s: %s %d alleles" % (k, os.path.split(v[0])[-1], v[1]))

    # Process genome sequences
    genomes = load_genomes()
    logger.info("Found sequences for %d genomes:" % len(genomes))
    for k, v in genomes.items():
        logger.info("%s: %s" % (k, os.path.split(v)[-1]))

    # Make the output directory (if needed)
    make_outdir()

    # BLASTN allele sequences against the input genome
    blastoutfiles = run_blast(alleles, genomes)

    # Create a dataframe with genomes as rows, and alleles as columns.
    # We'll populate these as we process our BLAST output
    df = pd.DataFrame(index=genomes.keys(), columns=alleles.keys())

    # Assign alleles to isolates
    df = assign_alleles(blastoutfiles, df)
    print df
