#README.md (scripts/bioinformatics)

## Overview
This repository contains a set of scripts for use in bioinformatics work. The dependencies (and possibly languages) may vary for each script. Please refer to the individual **README** details for each script.

## Scripts
The current set of scripts includes:

* [`calculate_ani.py`](#calculate_ani): calculates whole genome similarity measures ANIb, ANIm, and TETRA, with tabular text and graphical output.
* [`draw_gd_all_core.py`](#draw_gd_all_core)
* [`find_asm_snps.py`](#find_asm_snps)
* [`get_NCBI_cds_from_protein.py`](#get_NCBI_cds_from_protein)
* [`restrict_long_contigs.py`](#restrict_long_contigs)
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
* **Python** <http://www.python.org> (2.6+ required for `multiprocessing`)


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
* **Python** <http://www.python.org> (2.6+ required for `multiprocessing`)


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
