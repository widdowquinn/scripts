#!/usr/bin/env python
#
# genbank_get_genomes_by_taxon.py
#
# A script that takes an NCBI taxonomy identifier (or string, though this is
# not reliable for taxonomy tree subgraphs...) and downloads all genomes it 
# can find from NCBI in the corresponding taxon subgraph with the passed
# argument as root.
#
# (c) TheJames Hutton Institute 2015
# Author: Leighton Pritchard

import logging
import logging.handlers
import os
import shutil
import sys
import time
import traceback

from argparse import ArgumentParser
from collections import defaultdict

from Bio import Entrez, SeqIO

# Parse command-line
def parse_cmdline(args):
    """Parse command-line arguments"""
    parser = ArgumentParser(prog="genbacnk_get_genomes_by_taxon.py")
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default=None,
                        help="Output directory")
    parser.add_argument("-t", "--taxon", dest="taxon",
                        action="store", default=None,
                        help="NCBI taxonomy ID")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Give verbose output")
    parser.add_argument("-f", "--force", dest="force",
                        action="store_true", default=False,
                        help="Force file overwriting")
    parser.add_argument("--noclobber", dest="noclobber",
                        action="store_true", default=False,
                        help="Don't nuke existing files")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    parser.add_argument("--format", dest="format",
                        action="store", default="gbk,fasta",
                        help="Output file format [gbk|fasta]")
    parser.add_argument("--email", dest="email",
                        action="store", default=None,
                        help="Email associated with NCBI queries")
    return parser.parse_args()


# Set contact email for NCBI
def set_ncbi_email():
    """Set contact email for NCBI."""
    Entrez.email = args.email
    logger.info("Set NCBI contact email to %s" % args.email)
    Entrez.tool = "genbank_get_genomes_by_taxon.py"
    

# Create output directory if it doesn't exist
def make_outdir():
    """Make the output directory, if required.
    This is a little involved.  If the output directory already exists,
    we take the safe option by default, and stop with an error.  We can,
    however, choose to force the program to go on, in which case we can
    either clobber the existing directory, or not.  The options turn out
    as the following, if the directory exists:
    DEFAULT: stop and report the collision
    FORCE: continue, and remove the existing output directory
    NOCLOBBER+FORCE: continue, but do not remove the existing output
    """
    if os.path.exists(args.outdirname):
        if not args.force:
            logger.error("Output directory %s would " % args.outdirname +
                         "overwrite existing files (exiting)")
            sys.exit(1)
        else:
            logger.info("--force output directory use")
            if args.noclobber:
                logger.warning("--noclobber: existing output directory kept")
            else:
                logger.info("Removing directory %s and everything below it" %
                            args.outdirname)
                shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s" % args.outdirname)
    try:
        os.makedirs(args.outdirname)   # We make the directory recursively
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if args.noclobber and args.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)


# Get assembly UIDs for the root taxon
def get_asm_uids(taxon_uid):
    """Returns a set of NCBI UIDs associated with the passed taxon.

    This query at NCBI returns all assemblies for the taxon subtree
    rooted at the passed taxon_uid.
    """
    asm_ids = set()  # Holds unique assembly UIDs
    query = "txid%s[Organism:exp]" % taxon_uid
    logger.info("ESearch for %s" % query)
    
    # Perform initial search with usehistory
    handle = Entrez.esearch(db="assembly", term=query, format="xml",
                            usehistory="y")
    record = Entrez.read(handle)
    result_count = int(record['Count'])
    logger.info("Entrez ESearch returns %d assembly IDs" % result_count)
    
    # Recover all child nodes
    batch_size = 250
    for start in range(0, result_count, batch_size):
        tries, success = 0, False
        while not success and tries < 20:
            try:
                batch_handle = Entrez.efetch(db="assembly", retmode="xml",
                                             retstart=start,
                                             retmax=batch_size,
                                             webenv=record["WebEnv"],
                                             query_key=record["QueryKey"])
                batch_record = Entrez.read(batch_handle)
                asm_ids = asm_ids.union(set(batch_record))
                success = True
            except:
                tries += 1
                logger.warning("Entrez batch query failed (#%d)" % tries)
        if not success:
            logger.error("Too many download attempt failures (exiting)")
            sys.exit(1)
    logger.info("Identified %d unique assemblies" % len(asm_ids))
    return asm_ids

# Get contig UIDs for a specified assembly UID
def get_contig_uids(asm_uid):
    """Returns a set of NCBI UIDs associated with the passed assembly.

    The UIDs returns are for assembly_nuccore_insdc sequences - the 
    assembly contigs."""
    logger.info("Finding contig UIDs for assembly %s" % asm_uid)
    contig_ids = set()  # Holds unique contig UIDs
    links = Entrez.read(Entrez.elink(dbfrom="assembly", db="nucleotide",
                                     retmode="gb", from_uid=asm_uid))
    contigs = [l for l in links[0]['LinkSetDb'] \
               if l['LinkName'] == 'assembly_nuccore_insdc'][0]
    contig_uids = set([e['Id'] for e in contigs['Link']])
    logger.info("Identified %d contig UIDs" % len(contig_uids))
    return contig_uids

# Run as script
if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

    # Set up logging
    logger = logging.getLogger('genbank_get_genomes_by_taxon.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Was a logfile specified? If so, use it
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

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info("genbank_get_genomes_by_taxon.py: %s" % time.asctime())
    logger.info("command-line: %s" % ' '.join(sys.argv))
    logger.info(args)

    # Have we got an output directory? If not, exit.
    if args.email is None:
        logger.error("No email contact address provided (exiting)")
        sys.exit(1)
    set_ncbi_email()

    # Have we got an output directory? If not, exit.
    if args.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    make_outdir()
    logger.info("Output directory: %s" % args.outdirname)

    # We might have more than one taxon in a comma-separated list
    taxon_ids = args.taxon.split(',')
    logger.info("Passed taxon IDs: %s" % ', '.join(taxon_ids))    

    # Get all NCBI assemblies for each taxon UID
    asm_dict = defaultdict(set)
    for tid in taxon_ids:
        asm_dict[tid] = get_asm_uids(tid)
    for tid, asm_uids in asm_dict.items():
        logger.info("Taxon %s: %d assemblies" % (tid, len(asm_uids)))

    # Get links to the nucleotide database for each assembly UID
    contig_dict = defaultdict(set)
    for tid, asm_uids in asm_dict.items():
        for asm_uid in asm_uids:
            contig_dict[asm_uid] = get_contig_uids(asm_uid)
    for asm_uid, contig_uids in contig_dict.items():
        logger.info("Assembly %s: %d contigs" % (asm_uid, len(contig_uids)))
