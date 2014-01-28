#README.md (scripts/bioinformatics)

## Overview
This repository contains a set of scripts for use in bioinformatics work. The dependencies (and possibly languages) may vary for each script. Please refer to the individual **README** details for each script.

## Scripts
The current set of scripts includes:

* [`calculate_ani.py`](#calculate_ani): calculates whole genome similarity measures ANIb, ANIm, and TETRA, with tabular text and graphical output.
* [`draw_gd_all_core.py`](#draw_gd_all_core)
* [`find_asm_snps.py`](#find_asm_snps)
* [`get_NCBI_cds_from_protein.py`](#get_NCBI_cds_from_protein) Given input of protein sequences with suitably-formatted identifiers, retrieves the corresponding coding sequence from NCBI, using Entrez.
* [`restrict_long_contigs.py`](#restrict_long_contigs) From a directory of FASTA files, generates a new directory of corresponding FASTA files where all sequences shorter than a specified length have been removed.
* [`run_signalp.py`](#run_signalp): splits large files and parallelises for input to `SignalP`.
* [`run_tmhmm.py`](#run_tmhmm): splits large files and parallelises for input to `TMHMM`.

## Script READMEs

### <a name="calculate_ani">`calculate_ani.py`</a>

This script calculates Average Nucleotide Identity (ANI) according to one of a number of alternative methods described in, e.g.

* Richter M, Rossello-Mora R (2009) Shifting the genomic gold standard for the prokaryotic species definition. Proc Natl Acad Sci USA 106: 19126-19131. doi:10.1073/pnas.0906412106. (ANI1020, ANIm, ANIb)
* Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, et al. (2007) DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. Int J Syst Evol Micr 57: 81-91. doi:10.1099/ijs.0.64483-0.

ANI is proposed to be the appropriate *in silico* substitute for DNA-DNA 
hybridisation (DDH), and so useful for delineating species boundaries. A 
typical percentage threshold for species boundary in the literature is 95% 
ANI (e.g. Richter et al. 2009).

All ANI methods follow the basic algorithm:

- Align the genome of organism 1 against that of organism 2, and identify the matching regions
- Calculate the percentage nucleotide identity of the matching regions, as an average for all matching regions

Methods differ on: (1) what alignment algorithm is used, and the choice of parameters (this affects the aligned region boundaries); (2) what the input is for alignment (typically either fragments of fixed size, or the most complete assembly available).

* **ANIm**: uses MUMmer (NUCmer) to align the input sequences.
* **ANIb**: uses BLASTN to align 1000nt fragments of the input sequences
* **TETRA**: calculates tetranucleotide frequencies of each input sequence

This script takes as input a directory containing a set of correctly-formatted FASTA multiple sequence files. All sequences for a single organism should be contained in only one sequence file. The names of these files are used for identification, so it would be advisable to name 
them sensibly.

Output is written to a named directory. The output files differ depending on the chosen ANI method.

* **ANIm**: MUMmer/NUCmer .delta files, describing the sequence alignment; tab-separated format plain text tables describing total alignment lengths, and total alignment percentage identity
* **ANIb**: FASTA sequences describing 1000nt fragments of each input sequence; BLAST nucleotide databases - one for each set of fragments; and BLASTN output files (tab-separated tabular format plain text) - one for each pairwise comparison of input sequences. There are potentially a lot of intermediate files.
* **TETRA**: Tab-separated text file describing the Z-scores for each tetranucleotide in each input sequence.

In addition, all methods produce a table of output percentage identity (ANIm and ANIb) or correlation (TETRA), between each sequence.

If graphical output is chosen, the output directory will also contain PDF files representing the similarity between sequences as a heatmap with row and column dendrograms.

#### Usage

```
calculate_ani.py [options]

Options:
   -h, --help            show this help message and exit
   -o OUTDIRNAME, --outdir=OUTDIRNAME
                         Output directory
   -i INDIRNAME, --indir=INDIRNAME
                         Input directory name
   -v, --verbose         Give verbose output
   -f, --force           Force file overwriting
   -s, --fragsize        Sequence fragment size for ANIb
   --skip_nucmer         Skip NUCmer runs, for testing (e.g. if output already
                         present)
   --skip_blast          Skip BLAST runs, for testing (e.g. if output already
                         present)
   --noclobber           Don't nuke existing files
   -g, --graphics        Generate heatmap of ANI
   -m METHOD, --method=METHOD
                         ANI method
   --maxmatch            Override MUMmer settings and allow all matches in 
                         NUCmer
   --nucmer_exe=NUCMER_EXE
                         Path to NUCmer executable
   --blast_exe=BLAST_EXE
                         Path to BLASTN+ executable
   --makeblastdb_exe=MAKEBLASTDB_EXE
                         Path to BLAST+ makeblastdb executable
```

#### Dependencies

(mandatory)

* **Biopython** <http://www.biopython.org>
* **BLAST+** executable in the `$PATH`, or available on the command line (**ANIb**) <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>
* **MUMmer** executables in the $PATH, or available on the command line (**ANIm**) <http://mummer.sourceforge.net/>

(optional, only required for graphical output)

* **R** with shared libraries installed on the system <http://cran.r-project.org/>
* **Rpy2** <http://rpy.sourceforge.net/rpy2.html>

### <a name="draw_gd_all_core">`draw_gd_all_Core.py`</a>

This script was made available to provide demonstration GenomeDiagram usage code, as part of a [Twitter conversation](https://twitter.com/widdowquinn/status/411165587297951744/photo/1). The script renders the image in the photograph, which is the same as that at the foot of [this blog post](http://armchairbiology.blogspot.co.uk/2012/09/the-colours-man-colours.html).

**NOTE: The script will not run without modification unless you have both the required data, and the as yet unreleased `iadhore.py` module.**

#### Usage

```
draw_gd_all_core.py
```

#### Dependencies

* **Biopython** <http://www.biopython.org>
* **Reportlab** <http://www.reportlab.org> 
* **ColorSpiral** <https://github.com/widdowquinn/ColorSpiral> (this has since been incorporated into Biopython, but this script relies on the standalone library)
* **iAdhore.py** (not yet publicly available)
* **relevant data** (not yet publicly available)


### <a name="get_NCBI_cds_from_protein">`get_NCBI_cds_from_protein.py`</a>

Given input of protein sequences with suitably-formatted identifier strings (e.g. those in the NCBI nr database, having `>gi|[0-9]|.*` as their identifier), this script uses Entrez tools to find corresponding coding sequences where possible.

If the `--sto` option is used, sequence identifiers of the format `<seqid>|/<start>-<end> >gi|[0-9].*|/[0-9]*-[0-9]*` will be interpreted as Stockholm sequence headers, as described at <http://sonnhammer.sbc.su.se/Stockholm.html>

Unless the `--keepcache` option is passed, all data obtained from Entrez retained in local cache files, to aid debugging and limit network traffic, will be deleted when the script completes. By default, caches will have a  filestem reflecting the time that the script was run.  Cache filestems can be specified by the `-c` or `--cachestem` option; if the cache already exists, it  will be reused, and not overwritten.  This is useful for debugging.

An email address must be provided for Entrez, using the `-e` or `--email` option.

Queries will be batched in groups of 500 using EPost by default.  Batch size can be specified with `-b` or `--batchsize`.

Sequences can be provided via stdin, and output sequences are written to  stdout, by default.  Specifying a filename with `-o` or `--outfilename` writes output to that file.

The script runs as follows:

* Collect input sequences with valid sequence identifiers
* Check identifiers against an ELink cache of `protein_nuccore` searches. Where there is an entry in the cache, use this.  For sequences with no cache entry, compile identifiers in batches and submit queries in a `protein_nuccore `search with EPost. Add the recovered results to the ELink cache.
* Collect identifiers from the Elink results, and check against a cache of   GenBank headers.  For identifiers with no cache entry, compile identifiers in batches and submit an EFetch for the GenBank headers using EPost, caching the results.
* For each of the query protein sequences, choose the shortest potential nucleotide coding sequence from the GenBank header cache, and check its    identifier against the full GenBank record cache.  For identifiers with no cache entry, compile them in batches and submit an EFetch for the complete GenBank record, and cache the results.
* For each of the query protein sequences, get the corresponding CDS from the cached full GenBank sequence on the basis of a matching `GI:[0-9]*` value in a `db_xref` qualifier (and/or a matching `protein_id` qualifier if there's a suitable accession in the query sequence header.
* Extract the CDS nucleotide sequence from the parent sequence, and write to the output stream, with appropriate headers.

While the script runs, the important query and CDS data are kept in the results dictionary, keyed by query ID, with value a dictionary containing the sequence identifier, query AA sequence, query GI ID, and eventually the matching CDS feature, its sequence and, if the Stockholm header format is  set, the sequence that codes for the region specified in the query header.

#### Usage

```
Usage: get_NCBI_cds_from_protein.py [options]

Options:
  -h, --help            show this help message and exit
  -o OUTFILENAME, --outfile=OUTFILENAME
                        Output FASTA sequence filename
  -i INFILENAME, --infile=INFILENAME
                        Input FASTA sequence file
  -v, --verbose         Give verbose output
  -e EMAIL, --email=EMAIL
                        Entrez email
  -c CACHESTEM, --cachestem=CACHESTEM
                        Suffix for cache filestem
  -b BATCHSIZE, --batchsize=BATCHSIZE
                        Batchsize for EPost submissions
  -r RETRIES, --retries=RETRIES
                        Maximum number of Entrez retries
  -l LIMIT, --limit=LIMIT
                        Limit number of sequences processed (for testing)
  --keepcache           Keep cache files for debugging purposes
  --sto                 Parse .*/start-end identifiers as sequence regions
```

#### Dependencies

* **Biopython** <http://www.biopython.org>


### <a name="restrict_long_contigs">`restrict_long_contigs.py`</a>

A short script that takes as input a directory containing (many) FASTA files describing biological sequences, and writes to a new, named directory multiple FASTA files containing the same sequences, but restricted only to those sequences whose length is greater than a passed value.

Example usage: You have a directory with many sets of contigs from an assembly. This script will produce a new directory of the same data where the contig lengths are restricted to being greater than a specified length.

#### Usage

```
Usage: restrict_long_contigs.py [options] <input_directory> <output_directory>

Options:
  -h, --help            show this help message and exit
  -l MINLEN, --minlen=MINLEN
                        Minimum length of sequence
  -s SUFFIX, --filesuffix=SUFFIX
                        Suffix to indicate the file was processed
  -v, --verbose         Give verbose output
```

#### Dependencies

* **Biopython** <http://www.biopython.org>

### <a name="run_signalp">`run_signalp.py`</a>

This script takes a FASTA format file containing protein sequences as input, and runs a local copy of `signalp` (in the `$PATH`) on the contents, collecting the output generated with the `-short` option of `signalp`. 

The script splits the input sequence file into a number of smaller FASTA files suitable for distributed processing using, e.g. with Python's `multiprocessing` module. It prepares intermediate files with a maximum of 1000 sequences per file, to avoid falling foul of the undocumented upper input limit of `signalp`.  

The script runs `signalp` independently on  each split file, and the results are concatenated into a single output file. If no output file is specified, then the output file shares a stem with the input file. The organism type is specified at the command line in the same way as for `signalp`.

#### Usage ####

```
[python] run_signalp.py [euk|gram+|gram-] <FASTAfile> [-o|--outfilename <output file>]
```

#### Dependencies

* **signalp** <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp>
* **Python 2.6+** <http://www.python.org> (2.6+ required for `multiprocessing`)


### <a name="run_tmhmm">`run_tmhmm.py`</a>

This script takes a FASTA format file containing protein sequences as input, and runs a local copy of `tmhmm` (in the `$PATH`) on the contents, collecting the output generated with the `-short` option of `tmhmm`. 

The script splits the input sequence file into a number of smaller FASTA files suitable for distributed processing using, e.g. with Python's `multiprocessing` module. It prepares intermediate files with a maximum of 1000 sequences per file, to avoid falling foul of the undocumented upper input limit of `tmhmm`.  

The script runs `tmhmm` independently on  each split file, and the results are concatenated into a single output file. If no output file is specified, then the output file shares a stem with the input file.

#### Usage ####

```
run_tmhmm.py <FASTAfile> [-o|--outfilename <output file>]
```

#### Dependencies

* **tmhmm** <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm>
* **Python 2.6+** <http://www.python.org> (2.6+ required for `multiprocessing`)


## Licensing

Unless otherwise indicated in the script, all code is subject to the following agreement:

(c) The James Hutton Institute 2013
Author: Leighton Pritchard

Contact:
`leighton.pritchard@hutton.ac.uk`

Address: 
> Leighton Pritchard,
> Information and Computational Sciences,
> James Hutton Institute,
> Errol Road,
> Invergowrie,
> Dundee,
> DD6 9LH,
> Scotland,
> UK

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**; without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
