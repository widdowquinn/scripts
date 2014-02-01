#!/usr/bin/python
#
# get_NCBI_cds_from_protein.py
#
# Given input of protein sequences with suitably-formatted identifier strings
# (e.g. those in the NCBI nr database, having >gi|[0-9]|.* as their
# identifier), this script uses Entrez tools to find corresponding coding
# sequences where possible.
#
# If the --sto option is used, sequence identifiers of the format
# <seqid>|/<start>-<end> >gi|[0-9].*|/[0-9]*-[0-9]* will be interpreted as
# Stockholm sequence headers, as described at
# http://sonnhammer.sbc.su.se/Stockholm.html
#
# Unless the --keepcache option is passed, all data obtained from Entrez
# retained in local cache files, to aid debugging and limit network traffic,
# will be deleted when the script completes. By default, caches will have a
# filestem reflecting the time that the script was run.  Cache filestems can
# be specified by the -c or --cachestem option; if the cache already exists, it
# will be reused, and not overwritten.  This is useful for debugging.
#
# An email address must be provided for Entrez, using the -e or --email option.
#
# Queries will be batched in groups of 500 using EPost by default.  Batch size
# can be specified with -b or --batchsize.
#
# Sequences can be provided via stdin, and output sequences are written to
# stdout, by default.  Specifying a filename with -o or --outfilename writes
# output to that file.
#
# The script runs as follows:
#
# * Collect input sequences with valid sequence identifiers
# * Check identifiers against an ELink cache of protein_nuccore searches.
#   Where there is an entry in the cache, use this.  For sequences with no
#   cache entry, compile identifiers in batches and submit queries in a
#   protein_nuccore search with EPost.  Add the recovered results to the
#   ELink cache.
# * Collect identifiers from the Elink results, and check against a cache of
#   GenBank headers.  For identifiers with no cache entry, compile identifiers
#   in batches and submit an EFetch for the GenBank headers using EPost,
#   caching the results.
# * For each of the query protein sequences, choose the shortest potential
#   nucleotide coding sequence from the GenBank header cache, and check its
#   identifier against the full GenBank record cache.  For identifiers with
#   no cache entry, compile them in batches and submit an EFetch for the
#   complete GenBank record, and cache the results.
# * For each of the query protein sequences, get the corresponding CDS from the
#   cached full GenBank sequence on the basis of a matching GI:[0-9]* value
#   in a db_xref qualifier (and/or a matching protein_id qualifier if there's
#   a suitable accession in the query sequence header.
# * Extract the CDS nucleotide sequence from the parent sequence, and write
#   to the output stream, with appropriate headers.
#
# While the script runs, the important query and CDS data are kept in the
# results dictionary, keyed by query ID, with value a dictionary containing
# the sequence identifier, query AA sequence, query GI ID, and eventually the
# matching CDS feature, its sequence and, if the Stockholm header format is
# set, the sequence that codes for the region specified in the query header.
#
# (c) The James Hutton Institute 2012
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

from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from itertools import chain
from optparse import OptionParser

import logging
import logging.handlers

import StringIO
import os
import re
import sys
import time
import traceback


###
# FUNCTIONS

# Parse command-line
def parse_cmdline(args):
    """ Parse command-line arguments
    """
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outfile", dest="outfilename",
                      action="store", default=None,
                      help="Output FASTA sequence filename")
    parser.add_option("-i", "--infile", dest="infilename",
                      action="store", default=None,
                      help="Input FASTA sequence file")
    parser.add_option("-v", "--verbose", dest="verbose",
                      action="store_true", default=False,
                      help="Give verbose output")
    parser.add_option("-e", "--email", dest="email",
                      action="store", default=None,
                      help="Entrez email")
    parser.add_option("-c", "--cachestem", dest="cachestem",
                      action="store", default=time.strftime("%y%m%d%H%m%S"),
                      help="Suffix for cache filestem")
    parser.add_option("-b", "--batchsize", dest="batchsize",
                      action="store", default=500, type="int",
                      help="Batchsize for EPost submissions")
    parser.add_option("-r", "--retries", dest="retries",
                      action="store", default=10, type="int",
                      help="Maximum number of Entrez retries")
    parser.add_option("-l", "--limit", dest="limit",
                      action="store", default=None, type="int",
                      help="Limit number of sequences processed (for testing)")
    parser.add_option("--keepcache", dest="keepcache",
                      action="store_true", default=False,
                      help="Keep cache files for debugging purposes")
    parser.add_option("--sto", dest="stockholm",
                      action="store_true", default=False,
                      help="Parse .*/start-end identifiers as " +
                      "sequence regions")
    return parser.parse_args()


# Report last exception as string
def last_exception():
    """ Returns last exception as a string
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Create or verify existence of a cache file
def create_or_verify_cache(filename, dirname):
    """ Checks to see whether the passed file exists in the passed directory.
        If not, create the file.  Return path to the file.
    """
    cachename = os.path.join(dirname, filename)
    if os.path.isfile(cachename):   # File exists
        logger.info("Cache %s exists" % cachename)
    else:                           # File doesn't exist, so create it
        try:
            fhandle = open(cachename, 'w')
            fhandle.close()
            logger.info("Created cache file %s" % cachename)
        except:
            logger.error("Could not create cache file %s" % cachename)
            logger.error(last_exception())
            sys.exit(1)
    return cachename


# Establish caches
def create_caches(cachestem):
    """ Check whether caches with the passed cachename exist and, if not,
        creates them.
    """
    # Directory for caches - we're kind of assuming *nix - ensure it exists
    cachedir = os.path.expanduser("~/.lpscripts")
    if not os.path.exists(cachedir):
        try:
            logger.info("Creating directory for cache files in %s" % cachedir)
            os.mkdir(cachedir)
        except:
            logger.error("Could not create cache directory (exiting)")
            logger.error(last_exception())
            sys.exit(1)
    # Create each required cache: elink_cache_, gb_cache_, gbfull_cache_,
    # acc_cache_
    elink_cachename = \
        create_or_verify_cache("elink_cache_%s" % cachestem, cachedir)
    gb_cachename = create_or_verify_cache("gb_cache_%s" % cachestem, cachedir)
    gbfull_cachename = \
        create_or_verify_cache("gbfull_cache_%s" % cachestem, cachedir)
    acc_cachename = \
        create_or_verify_cache("acc_cache_%s" % cachestem, cachedir)
    # Return the cache filenames
    return elink_cachename, gb_cachename, gbfull_cachename, acc_cachename


# Process sequences from input stream
def process_sequences(fhandle, limit):
    """ We load in sequences from the passed stream of FASTA format sequences,
        but retain only those that have a suitably-formatted sequence ID:

        >gi\|[0-9]*\|

        That is, having a GI accession as the first accession in the identifier
    """
    logger.info("Parsing sequences")
    if limit:
        logger.info("Limiting to %d sequences" % limit)
    regex = re.compile('gi\|[0-9]*\|.*')  # Sequence identifiers must match
    resdict = {}                          # Dictionary to hold data
    attempts = 0
    for s in SeqIO.parse(fhandle, 'fasta'):
        attempts += 1
        if re.match(regex, s.id) is None:
            logger.warning("Sequence ID %s is unexpected format" % s.id)
        else:
            logger.info("Retaining sequence %s" % s.id)
            resdict[s.id] = {'gid': s.id.split('|')[1],
                             'sid': s.id, 'aa': s.seq}
        if limit and len(resdict) >= limit:
            break    # For testing
    logger.info("Loaded %d sequences (from %d attempts)" %
                (len(resdict), attempts))
    if not attempts:
        logger.error("No suitable sequences found (exiting)")
        logger.error(last_exception())
        sys.exit(1)
    return resdict


# Load elink_cache and return dictionary keyed by query ID, with values
# a list of match accessions
def load_elink_cache(filename):
    """ Load the contents of the elink cache, returning a dictionary keyed
        by query GI IDs, with values a list of match GI IDs.
    """
    try:
        data = [l.strip().split() for l in open(filename, 'rU') if
                len(l.strip())]
    except:
        logger.error("Could not parse ELink cache %s (exiting)" % filename)
        logger.error(last_exception)
        sys.exit(1)
    cachedict = {}
    for line in data:
        if len(line) > 1:
            cachedict[line[0]] = line[1:]
        else:
            cachedict[line[0]] = None
    return cachedict


# Load accession_cache and return dictionary keyed by GI ID, with values
# a list of GenBank accessions
def load_accession_cache(filename):
    """ Load the contents of the accession cache, returning a dictionary keyed
        by query GI IDs, with value the GenBank accession
    """
    try:
        data = [l.strip().split() for l in open(filename, 'rU') if
                len(l.strip())]
    except:
        logger.error("Could not parse accession cache %s (exiting)" % filename)
        logger.error(last_exception)
        sys.exit(1)
    cachedict = {}
    for line in data:
        if len(line) > 1:
            cachedict[line[0]] = line[1]
        else:
            cachedict[line[0]] = None
    return cachedict


# Load accession_cache and return dictionary keyed by GenBank ID, with values
# a dictionary with keys 'len' and 'gi' relating to the length and GI ID of
# the accession
def load_gb_cache(filename):
    """ Load the contents of the GenBank header summary cache, returning a
        dictionary keyed by GenBank accession, with value a dictionary
        containing the length and GI ID of the record.
    """
    try:
        data = [l.strip().split() for l in open(filename, 'rU') if
                len(l.strip())]
    except:
        logger.error("Could not parse GenBank summary cache %s (exiting)" %
                     filename)
        logger.error(last_exception)
        sys.exit(1)
    cachedict = {}
    for line in data:
        if len(line) != 3:
            logger.error("Problem with GenBank summary cache " +
                         "%s (%s) (exiting)" % (filename, line))
            logger.error(last_exception)
            sys.exit(1)
        cachedict[line[0]] = {'len': int(line[1]), 'gi': line[2]}
    return cachedict


# Extend the ELink cache by a query result
def extend_elink_cache(matches, filename, maxretries):
    """ Adds the passed protein_nuccore match to the ELink cache.  We assume
        that there was one query ID, and there's only one result in the list.

        If the query is not already in the cache, the data are added outright.
        If the query is in the cache, no new data is added

        The match format is
        <query id>\t<match id 1>\t<match id 2>...
    """
    match = matches[0]
    qid = match['IdList'][0]
    #logger.info("Parsing ELink cache %s" % filename)
    elink_match_ids = set(load_elink_cache(elink_cachename).keys())
    # If the query is already in the cache, return
    if qid in elink_match_ids:
        logger.warning("Query ID %s already in ELink cache" %
                       match['IdList'][0])
        return
    # Query not in cache, so find the appropriate accessions, and add a new
    # line to the end of cache
    try:
        match_ids = [i['Id'] for i in match['LinkSetDb'][0]['Link']]
    except IndexError:
        logger.warning("No ELink matches for %s" % qid)
        return
    entry_str = [qid] + match_ids
    try:
        cache_handle = open(filename, 'a')
        cache_handle.write("%s\n" % '\t'.join(entry_str))
        cache_handle.close()
        logger.info("Cached ELink matches for %s" % qid)
    except:
        logger.error("Could not cache ELink matches for %s " % qid +
                     "(exiting)")
        logger.error(last_exception())
        sys.exit(1)


# Extend the accession cache by an association between a GI ID and an accession
def extend_accession_cache(gi_id, accession, filename):
    """ Adds the passed GI ID and associated GenBank accession match to the
        accession cache.  We assume that there was one query ID, and there's
        only one GenBank accession for that GI ID.

        If the query is not already in the cache, the data are added outright.
        If the query is in the cache, no new data is added

        The match format is
        <GI id>\t<accession>\n
    """
    cached_gi_ids = set(load_accession_cache(filename).keys())
    # If the query is already in the cache, return
    if gi_id in cached_gi_ids:
        logger.warning("GI ID %s already in Accession cache" % gi_id)
        return
    # Query not in cache so add a new line to the end of cache
    try:
        cache_handle = open(filename, 'a')
        cache_handle.write("%s\t%s\n" % (gi_id, accession))
        cache_handle.close()
        logger.info("Cached ELink matches for %s" % gi_id)
    except:
        logger.error("Could not cache accession for %s " % gi_id +
                     "(exiting)")
        logger.error(last_exception())
        sys.exit(1)


# Use Entrez EUtils to find protein_nuccore links for each of the passed
# sequences
def cache_elinks(resdict, elink_cachename, maxretries):
    """ Loop over the passed results dictionary, checking the GI accessions of
        each against the passed elink cache data.  Those accessions without
        entries in the cache submitted to Entrez to obtain
        protein_nuccore matches.  These are written to the elink_cache.

        We cannot use batching with EPost here, as the returned link results
        are not associated with the corresponding query sequences.

        Also, EPost queries don't quite match EFetch direct queries:
        In [79]: history = Entrez.read(Entrez.epost("protein", id="15242228"))
        In [85]: Entrez.read(Entrez.elink(dbfrom="protein",
                 linkname="protein_nuccore",query_key=history['QueryKey'],
                 webenv=history['WebEnv']))
        Out[85]: [{u'LinkSetDb': [{u'DbTo': 'pubmed', u'Link': [{u'Id':
                  '240256493'}, {u'Id': '145358304'}], u'LinkName':
                  'protein_nuccore'}], u'DbFrom': 'protein', u'IdList':
                  ['15242228'], u'LinkSetDbHistory': [], u'ERROR': []}]
        In [86]: Entrez.read(Entrez.elink(dbfrom="protein",
                 linkname="protein_nuccore",id="15242228"))
        Out[86]: [{u'LinkSetDb': [{u'DbTo': 'nuccore', u'Link': [{u'Id':
                  '240256493'}, {u'Id': '145358304'}], u'LinkName':
                  'protein_nuccore'}], u'DbFrom': 'protein', u'IdList':
                  ['15242228'], u'LinkSetDbHistory': [], u'ERROR': []}]
    """
    # Check sequence IDs against the ELink cache
    logger.info("Checking ELink matches in cache")
    elink_match_ids = set(load_elink_cache(elink_cachename).keys())
    seqids = set([resdict[s]['gid'] for s in resdict])  # All query GI IDs
    # List the query GI IDs not in cache
    notincache = list(seqids.difference(elink_match_ids))
    # Use EPost to find protein_nuccore matches for those sequences not in
    # the cache
    if not len(notincache):
        logger.info("All input sequences found in ELink cache")
        return
    logger.info("%d queries not found in ELink cache" % len(notincache))
    logger.info("Finding ELink matches in protein_nuccore")
    # Run ELink submissions for each query.  We can't use EPost because of the
    # need to associate results with a single query
    for query_id in notincache:
        matches = elink_fetch_with_retries(query_id, "protein",
                                           "protein_nuccore", maxretries)
        extend_elink_cache(matches, elink_cachename, maxretries)


# Fetch ELink matches to passed query GI ID
def elink_fetch_with_retries(qid, dbfrom, linkname, maxretries):
    """ Runs an ELink search with the passed parameters, with a maximum
        number of retries.
        Returns the result parsed by Entrez.read.
    """
    tries = 0
    while tries < maxretries:
        try:
            matches = Entrez.read(Entrez.elink(dbfrom=dbfrom,
                                               linkname=linkname,
                                               id=qid))
            return matches
        except:
            tries += 1
            logger.error("EFetch submission failed (try %d/%d)" %
                         (tries, maxretries))
            logger.error(last_exception())
    # Won't get here unless there's failure
    logger.error("Could not complete EFetch submission (exiting)")
    sys.exit(1)


# Run an EFetch on a single ID with passed parameters
def efetch_with_retries(qid, db, rettype, retmode, maxretries):
    """ Run an Entrez EFetch for a single query ID, with passed parameters.
        Returns a handle to a completely buffered string, after a read()
        operation, sanity check, and retokenising.
    """
    tries = 0
    while tries < maxretries:
        try:
            data = Entrez.efetch(db=db, rettype=rettype,
                                 retmode=retmode, id=qid).read()
            if rettype in ['gb', 'gbwithparts'] and retmode == 'text':
                assert data.startswith('LOCUS')
            logger.info("Successful EFetch for %s query ID" % qid)
            # Return data string as stream
            return StringIO.StringIO(data)
        except:
            tries += 1
            logger.error("Query ID %s EFetch failed (%d/%d)" %
                         (qid, tries, maxretries))
            logger.error(last_exception())
    if success is False:
        logger.error("Query ID %s EFetch failed (exiting)" % qid)
        sys.exit(1)


# Cache GenBank accessions for each ELink ID
def cache_accessions(resdict, elink_cachename, acc_cachename, maxretries):
    """ Associate each ELink ID with the appropriate accession.
        We have to query by individual GI ID, as we want to associate each
        ID with its own accession, and an EPost will hide this.
    """
    t0 = time.time()
    logger.info("Caching GenBank IDs for ELink matches")
    # Get ELink match IDs
    result_gids = [v['gid'] for k, v in resdict.items()]
    cached_elink_data = load_elink_cache(elink_cachename)
    result_elink_keys = \
        set(result_gids).intersection(set(cached_elink_data.keys()))
    elink_ids = set(chain.from_iterable([cached_elink_data[gid] for gid in
                                         result_elink_keys]))
    #elink_ids = set(chain.from_iterable([v for v in \
    #                        load_elink_cache(elink_cachename).values()]))
    logger.info("Found %d unique ELink IDs" % len(elink_ids))
    logger.info("%.2fs" % (time.time() - t0))
    # First we check the accession cache and, if the GI ID is not in there,
    # we need to run an EFetch for the accession.
    accession_index = load_accession_cache(acc_cachename).keys()
    notincache = [e for e in elink_ids if e not in accession_index]
    if not len(notincache):
        logger.info("All ELink IDs in accession cache")
        return
    logger.info("%d ELink IDs not in accession cache" % len(notincache))
    logger.info("Associating %d ELink IDs with accessions" % len(notincache))
    for elink_id in notincache:
        acc = efetch_with_retries(elink_id, 'nucleotide', 'acc', 'text',
                                  maxretries).read().strip()
        extend_accession_cache(elink_id, acc, acc_cachename)


# Run EPost to get history for a batch of query IDs
def epost_history_with_retries(qids, db, maxretries):
    """ Run an Entrez EPost on a list of queryids, returning the generated
        history, as parsed by Entrez read
    """
    tries = 0
    while tries < maxretries:
        try:
            history = Entrez.read(Entrez.epost(db, id=','.join(qids)))
            logger.info("Successful EPost for %d queries" % len(qids))
            # Return history
            return history
        except:
            tries += 1
            logger.error("EPost failed for %d queries (%d/%d)" %
                         (len(qids), tries, maxretries))
            logger.error(last_exception())
    logger.error("EPost failed for %d queries (exiting)" % len(qids))
    sys.exit(1)


# Run EFetch to get results for a batch of query IDs
def efetch_history_with_retries(history, db, rettype, retmode, maxretries):
    """ Run an Entrez EFetch on a passed history, returning the result
        data as a tokenised string.
    """
    tries = 0
    while tries < maxretries:
        try:
            data = Entrez.efetch(db=db, rettype=rettype, retmode=retmode,
                                 webenv=history['WebEnv'],
                                 query_key=history['QueryKey']).read()
            if rettype in ['gb', 'gbwithparts'] and retmode == 'text':
                assert data.startswith('LOCUS')
            logger.info("Successful EFetch for history %s" %
                        history['QueryKey'])
            # Return histor
            return StringIO.StringIO(data)
        except:
            tries += 1
            logger.error("EFetch failed for history %s (retrying)" %
                         history['QueryKey'])
            logger.error(last_exception())
    logger.error("EFetch failed for history %s failed (exiting)" %
                 history['QueryKey'])
    sys.exit(1)


# Cache GenBank headers for each unique nucleotide GI ID in the ELink cache
def cache_gb_headers(resdict, elink_cachename, acc_cachename, gb_cachename,
                     batchsize, maxretries):
    """ Identify all unique nucleotide GI IDs in the ELink cache, and
        check whether they are already in the GenBank header cache.  Batch
        those that are not as EPost queries, then run EFetch searches for
        the 'gb' rettype.
    """
    t0 = time.time()
    logger.info("Caching GenBank headers for ELink matches")
    # We get a set of GenBank accessions in the cache, and the set of
    # GenBank accessions from the Elink and accession caches; we generate
    # a list of GenBank accessions that are not currently in the cache, and
    # that correspond to our query sequences.
    # We get a list of result GI IDs, and find their intersection with the
    # keys for entries in the ELink cache
    result_gids = [v['gid'] for k, v in resdict.items()]
    cached_elink_data = load_elink_cache(elink_cachename)
    result_elink_keys = \
        set(result_gids).intersection(set(cached_elink_data.keys()))
    # Then we flatten the list of values for this set of keys, to get a
    # list of potentially coding nucleotide sequences for our queries
    elink_ids = set(chain.from_iterable([cached_elink_data[gid] for gid in
                                         result_elink_keys]))
    # We get a set of GenBank IDs for already cached GenBank summaries,
    # and the GI ID:GenBank accession tally data from cache
    cached_gb_ids = set(load_gb_cache(gb_cachename).keys())
    cached_acc_data = load_accession_cache(acc_cachename)
    # Identify which records we have cached accession information for,
    # but not a cached summary - this is a superset of the sequences we want
    gb_notincache = set(cached_acc_data.values()).difference(cached_gb_ids)
    # Restrict this set of records with no summary in the cache only to those
    # records that are associated with one of our queries
    elink_notincache = [k for k, v in cached_acc_data.items() if
                        v in gb_notincache and k in elink_ids]
    logger.info("%.2fs" % (time.time() - t0))
    if not len(elink_notincache):  # Nothing left to do
        logger.info("All recovered ELink IDs in GenBank cache")
        return
    else:
        logger.info("%d unique ELink IDs not in cache" %
                    len(elink_notincache))
    # Batch the unique, uncached IDs and run EPost on each batch
    histories = []
    for batch in [elink_notincache[i:i+batchsize] for i in
                  range(0, len(elink_notincache), batchsize)]:
        histories.append(epost_history_with_retries(batch, 'nucleotide',
                                                    maxretries))
    # Run EFetch searches with each history
    for history in histories:
        logger.info("Fetching batch result %s/%d" %
                    (history['QueryKey'], len(histories)))
        # CACHING FIX?
        # Instead of only retrieving 'gb' type here, retrieve full record,
        #   but modify summary cache extension to write full record to
        #   the full gbcache, if it's not already in there.  We would index
        #   the full cache once (above) to get the keys, and pass keys
        #   with the extend GB summary function
        gbhandle = efetch_history_with_retries(history, 'nucleotide',
                                               'gb', 'text', maxretries)
        extend_gb_summary_cache(gb_cachename, gbhandle)
    logger.info("%.2fs" % (time.time() - t0))


# Cache full GenBank records for the shortest coding sequence record
# corresponding to each ELink GI ID
def cache_gb_full(resdict, elink_cachename, acc_cachename, gb_cachename,
                  gbfull_cachename, batchsize, maxretries):
    """ For each query GI ID in the ELink cache, we find the shortest GenBank
        record from the GenBank header cache, and download the full sequence.
    """
    t0 = time.time()
    logger.info("Caching full GenBank headers")
    # Index/load data from caches
    logger.info("Indexing caches")
    elink_data = load_elink_cache(elink_cachename)
    elink_accessions = load_accession_cache(acc_cachename)
    #gb_index = SeqIO.index(gb_cachename, 'gb')
    gb_index = load_gb_cache(gb_cachename)
    logger.info("%.2fs" % (time.time()-t0))
    # A set of GI numbers corresponding to the cached full GenBank records
    logger.info("Identifying ELink IDs for cached records")
    gbfull_keyset = set(SeqIO.index(gbfull_cachename, 'gb').keys())
    accset = set([a for e, a in
                  elink_accessions.items()]).intersection(gbfull_keyset)
    gbfull_index = set([e for e, a in elink_accessions.items() if
                        a in accset])
    logger.info("Found %d accessions in full GenBank cache" %
                len(gbfull_index))
    # Loop over query IDs, and find the shortest of the possible GenBank
    # coding sequence records
    logger.info("%.2fs" % (time.time()-t0))
    logger.info("Finding shortest GenBank Records")
    shortest_genbank_ids = set()
    query_ids = [v['gid'] for k, v in resdict.items() if v['gid'] in
                 elink_data]
    for qid, match_ids in \
            [(qid, match_ids) for qid, match_ids in elink_data.items()
             if qid in query_ids]:
        match_accessions = [elink_accessions[i] for i in match_ids]
        records = sorted([(gb_index[acc]['len'], gb_index[acc]['gi'])
                          for acc in match_accessions])
        # Find GI number of shortest record, and add to list
        shortest_genbank_ids.add(records[0][1])
    logger.info("%.2fs" % (time.time()-t0))
    logger.info("Found %d unique shortest GenBank records" %
                len(shortest_genbank_ids))
    # We need to restrict the disclosed GenBank IDs to those not already in
    # the cache
    notincache = list(shortest_genbank_ids.difference(gbfull_index))
    if not len(notincache):
        logger.info("All shortest full GenBank records already in cache")
        return
    logger.info("Caching %d shortest full GenBank records" % len(notincache))
    # Run batch EPost fetch of the full GenBank records, adding these to the
    # full GenBank cache
    histories = []
    for batch in [notincache[i:i+batchsize] for i in
                  range(0, len(notincache), batchsize)]:
        histories.append(epost_history_with_retries(batch, 'nucleotide',
                                                    maxretries))
    # Run EFetch searches with each history
    logger.info("%.2fs" % (time.time()-t0))
    for history in histories:
        logger.info("Fetching batch result %s/%d" %
                    (history['QueryKey'], len(histories)))
        gbdata = efetch_history_with_retries(history, 'nucleotide',
                                             'gbwithparts', 'text',
                                             maxretries).read()
        extend_gb_cache(gbfull_cachename, gbdata)
    logger.info("%.2fs" % (time.time()-t0))


# Add passed GenBank summary data to the named GenBank cache
def extend_gb_summary_cache(filename, handle):
    """ Add summary of the passed GenBank handle to the cache in the passed
        filename.
    """
    # CACHING FIX?
    # Here, we check each record anyway, but if we have a list of keys from the
    # full GenBank cache index, we can also write the full record to the
    # end of the full GenBank cache.  This would save a future lookup for that
    # record.
    try:
        cache_handle = open(filename, 'a')
        #cache_handle.write("%s\n" % data)
        # Rather than write the whole record, just keep the data we'll use
        # later.
        count = 0
        for record in SeqIO.parse(handle, 'gb'):
            cache_handle.write("%s\t%s\t%s\n" % (record.id, len(record),
                                                 record.annotations['gi']))
            count += 1
        cache_handle.close()
        logger.info("Extended cache %s with %d record summaries" %
                    (filename, count))
    except:
        logger.error("Could not extend cache %s" % filename)
        logger.error(last_exception())
        sys.exit(1)


# Add passed GenBank data to the named GenBank cache
def extend_gb_cache(filename, data):
    """ Add the passed data (as a string) to the cache in the passed
        filename.
    """
    try:
        cache_handle = open(filename, 'a')
        cache_handle.write("%s\n" % data)
        cache_handle.close()
        logger.info("Extended cache %s with %d bytes of data" %
                    (filename, len(data)))
    except:
        logger.error("Could not extend cache %s" % filename)
        logger.error(last_exception())
        sys.exit(1)


# Given the caches of ELink data, corresponding GenBank accessions, and full
# GenBank records, extract the CDS sequence for each query
def extract_cds(resdict, elink_cachename, acc_cachename, gbfull_cachename,
                stockholm):
    """ Using the data in the ELink, accession and full GenBank record caches,
        extract a CDS nucleotide sequence for each query.
    """
    logger.info("Extracting CDS for each query")
    if stockholm:
        logger.info("Expecting headers in Stockholm format")
    # Load/index caches
    logger.info("Indexing caches")
    elink_data = load_elink_cache(elink_cachename)
    elink_accessions = load_accession_cache(acc_cachename)
    gb_index = SeqIO.index(gbfull_cachename, 'gb')
    # Loop over query sequences, and find the appropriate CDS feature in the
    # correct GenBank record
    logger.info("Locating CDS features")
    for qhead, vals in resdict.items():
        logger.info("Finding CDS feature for %s" % qhead)
        qid = vals['gid']
        try:
            acc = [elink_accessions[l] for l in elink_data[qid] if
                   elink_accessions[l] in gb_index]
        except KeyError:
            logger.error("No GenBank records in accession cache for %s" % acc)
            continue
        if len(acc) > 1:
            logger.error("More than one GenBank record in accession " +
                         "cache for %s (%s)" % (qhead, acc) + " (exiting)")
            sys.exit(1)
        gbdata = gb_index[acc[0]]
        # We obtain the coding feature as a SeqFeature, and the coding sequence
        # as a Seq.  We need both because the feature can contain information
        # necessary to recover the actual coding sequence, such as an offset
        # qualifier.
        coding_feature = find_feature_by_gi(qid, gbdata, 'CDS')
        if not coding_feature:
            logger.warning("Could not find feature %s in %s" %
                           (qid, acc))
            continue
        else:
            if 'locus_tag' in coding_feature.qualifiers:
                logger.info("Found coding feature, locus_tag: %s" %
                            coding_feature.qualifiers['locus_tag'][0])
                results[qhead]['match_id'] =
                coding_feature.qualifiers['locus_tag'][0]
            elif 'gene' in coding_feature.qualifiers:
                logger.info("Found coding feature, gene: %s" %
                            coding_feature.qualifiers['gene'][0])
                results[qhead]['match_id'] = \
                    coding_feature.qualifiers['gene'][0]
            else:
                logger.info("Found coding feature, protein_id: %s" %
                            coding_feature.qualifiers['protein_id'][0])
                results[qhead]['match_id'] = \
                    coding_feature.qualifiers['protein_id'][0]
        # Add the complete CDS feature and sequence to the results dictionary
        results[qhead]['cds_feature'] = coding_feature
        results[qhead]['cds'] = coding_feature.extract(gbdata.seq)
        # Test whether the recovered CDS matches the passed query.  Report
        # whether it does, and keep the test_coding result in the results
        # dictionary
        query_cds = test_coding(results[qhead], stockholm)
        if query_cds:
            logger.info("CDS regenerates protein sequence for %s" % qhead)
        else:
            logger.warning("Possible issue with CDS translation for %s" %
                           qhead)
        results[qhead]['query_cds'] = query_cds
    return results


# Use a GI ID to extract a feature from a SeqRecord on the basis of a
# db_xref qualifier
def find_feature_by_gi(gid, record, ftype):
    """ Loops over the ftype features in the passed SeqRecord, checking
        db_xref qualifiers for a match to the passed gid.
        Returns the first feature identified, or None if no feature found.
    """
    for feature in [f for f in record.features if f.type == ftype]:
        try:
            if 'GI:%s' % gid in feature.qualifiers['db_xref']:
                return feature
        except KeyError:
            continue
    return None


# Test whether the recovered nt sequence actually codes for the aa sequence
# Not general - specific to our result dictionary structure
def test_coding(result, stockholm):
    """ We translate the ['cds'] sequence for the passed result, and compare
        that to the ['aa'] sequence, returning the appropriate nt sequence
        for a match, and False otherwise.
        We account for ID/start-end notation, only translating the specified
        region, if the ID line is valid Stockholm format.
    """
    # Translate sequence - making sure to translate account for a codon_start
    # offset, if there is one.
    if 'codon_start' in result['cds_feature'].qualifiers:
        if result['cds_feature'].qualifiers['codon_start'][0] == '1':
            cds = result['cds'].translate()
        else:
            offset = int(result['cds_feature'].qualifiers['codon_start'][0])-1
            cds = result['cds'][offset:].translate()
    else:
        cds = result['cds'].translate()
    # Do we need to account for Stockholm format
    if not stockholm:
        if str(cds) == str(result['aa']):
            return result['cds']
        else:
            logger.warning("CDS translation issue:\n%s\n%s\n" %
                           (result['aa'], cds))
            return None
    # Are we dealing with a sequence fragment or a whole sequence
    elif re.match('gi\|.*/[0-9]*-[0-9]*', result['sid']):  # fragment
        logger.info("Stockholm sequence fragment %s" % result['sid'])
        fstart, fend = tuple([int(e) for e in
                              result['sid'].split('/')[-1].split('-')])
        aa_frag = cds[fstart-1:fend]
        if str(aa_frag) == str(result['aa']):
            return result['cds'][(fstart-1)*3:(fend*3)+1]
        else:
            logger.warning("CDS translation issue:\n%s\n%s\n" %
                           (result['aa'], aa_frag))
    else:
        logger.warning("Stockholm format specified, header not Stockholm: %s" %
                       result['sid'])
        return None

###
# SCRIPT

if __name__ == '__main__':

    # Parse command-line
    # options are all options - no arguments
    options, args = parse_cmdline(sys.argv)

    # Set email address for remote Entrez
    Entrez.email = options.email

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    logger = logging.getLogger('CDS_from_AA_at_NCBI.py')
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

    # Do we have an input file?  No? Then use stdin
    if options.infilename is None:
        infhandle = sys.stdin
        logger.info("Using stdin for input")
    else:
        logger.info("Using %s for input" % options.infilename)
        try:
            infhandle = open(options.infilename, 'rU')
        except:
            logger.error("Could not open input file: %s (exiting)" %
                         options.infilename)
            logger.error(''.join(traceback.format_exception(
                sys.last_type,
                sys.last_value,
                sys.last_traceback)))
            sys.exit(1)

    # Do we have an output file?  No? Then use stdout
    if options.outfilename is None:
        outfhandle = sys.stdout
        logger.info("Using stdout for output")
    else:
        logger.info("Using %s for output" % options.outfilename)
        try:
            outfhandle = open(options.outfilename, 'w')
        except:
            logger.error("Could not open output file: %s (exiting)" %
                         options.outfilename)
            logger.error(''.join(traceback.format_exception(
                sys.last_type,
                sys.last_value,
                sys.last_traceback)))
            sys.exit(1)

    # Set up caches for ELink data, GenBank headers, GenBank full records,
    # and GenBank accessions/GI numbers
    elink_cachename, gb_cachename, gbfull_cachename, acc_cachename = \
        create_caches(options.cachestem)

    # Parse sequences in the input stream, keeping those with correctly-
    # formatted identifiers.  We obtain a dictionary keyed by the sequence
    # identifier, with value a dictionary, keyed with:
    # * gid: the GI accession number
    # * sid: the sequence identifier
    # * aa:  the amino acid sequence
    results = process_sequences(infhandle, options.limit)

    # Collect protein_nuccore matches from the cache/Entrez
    cache_elinks(results, elink_cachename, options.retries)

    # Collect GenBank accession data for the local cache
    cache_accessions(results, elink_cachename, acc_cachename, options.retries)

    # Collect GenBank headers for local cache
    cache_gb_headers(results, elink_cachename, acc_cachename, gb_cachename,
                     options.batchsize, options.retries)

    # Collect full GenBank records for shortest coding sequences in local cache
    cache_gb_full(results, elink_cachename, acc_cachename, gb_cachename,
                  gbfull_cachename, options.batchsize, options.retries)

    # Extract coding sequence from full GenBank records
    results = extract_cds(results, elink_cachename, acc_cachename,
                          gbfull_cachename, options.stockholm)

    # Write coding sequences to file
    if options.outfilename:
        logger.info("Writing CDS nucleotide sequences to %s" %
                    options.outfilename)
    # Some queries didn't have ELink matches, so we lose these...
    results_with_sequences = [r for r in results.values() if 'query_cds' in r]
    # ... and others have 'None' as a sequence, due to some sort of
    # translation issue - we lose these, too.
    SeqIO.write([SeqRecord(id=r['sid'], seq=r['query_cds'],
                           description=
                           'Coding sequence for %s derived from %s' %
                           (r['sid'], r['match_id']))
                 for r in results_with_sequences if r['query_cds']],
                outfhandle, 'fasta')
    outfhandle.close()
