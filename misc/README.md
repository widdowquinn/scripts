#README.md (scripts/misc)

## Overview
This repository contains a set of scripts for a bunch of disparate purposes. They may or may not be useful to you - I hope they might be, though. The dependencies (and possibly languages) may vary for each script. Please refer to the individual **README** details for each script.

## Scripts
The current set of scripts includes:

* [`tabular_to_wikitable.py`](#tabular_to_wikitable): converts a plain text tab-separated table to a MediaWiki wikitable

## Script READMEs

### <a name="tabular_to_wikitable">`tabular_to_wikitable.py`</a>

This script converts a plain text tab-separated table to MediaWiki wikitable markup. It will read/write from/to STDIN and STDOUT unless files are indicated on the command-line. A title, if added on the command-line, will be added to the table. The first line of the table will be used for column headers, unless they are passed on the command-line.

#### Usage

```
usage: tabular_to_mediawiki.py [-h] [-o OUTFILENAME] [-i INFILENAME] [-v]
                               [--header HEADER] [-t TITLE] [-s] [-c]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILENAME, --outfile OUTFILENAME
                        Output MediaWiki format table
  -i INFILENAME, --infile INFILENAME
                        Input tab-separated plaintext table
  -v, --verbose         Give verbose output
  --header HEADER       If not passed a string of comma-separated headers,
                        assumes first line of input tab-separated file is the
                        header line
  -t TITLE, --title TITLE
                        Set the title of the resulting MediaWiki table
  -s, --sortable        Make the MediaWiki table sortable
  -c, --collapsible     Make the MediaWiki table collapsible
```

Simple conversion. You just get a wikitable, and the script assumes your first line is column headers:

```
$ cat testdata/table.tab | ./tabular_to_wikitable.py
{|class="wikitable"
|+
! header1 !! header2 !! header3
|-
| egg || chips || spam
|-
| spam || spam || spam
|-
| Norwegian || Blue || parrot
|}
$ ./tabular_to_wikitable.py -i testdata/table.tab -o testdata/wikitable.txt
$ more testdata/wikitable.txt 
{|class="wikitable"
|+
! header1 !! header2 !! header3
|-
| egg || chips || spam
|-
| spam || spam || spam
|-
| Norwegian || Blue || parrot
|}
```

Add fancy classes - sortable and collapsible only, for now:

```
$ cat testdata/table.tab | ./tabular_to_wikitable.py --sortable --collapsible
{|class="wikitable sortable mw-collapsible"
|+
! header1 !! header2 !! header3
|-
| egg || chips || spam
|-
| spam || spam || spam
|-
| Norwegian || Blue || parrot
|}
$ cat testdata/table.tab | ./tabular_to_wikitable.py -s -c -v
INFO: Namespace(collapsible=True, header=None, infilename=None, outfilename=None, sortable=True, title=None, verbose=True)
INFO: Using stdin for input
INFO: Using stdout for output
INFO: Read 4 lines from input
INFO: Table appears to contain 3 columns
{|class="wikitable sortable mw-collapsible"
|+
! header1 !! header2 !! header3
|-
| egg || chips || spam
|-
| spam || spam || spam
|-
| Norwegian || Blue || parrot
|}
```

Handling too few/too many headers. If there are too few headers, the headers are extended with *colNN* placeholders. If there are too many, the data are **not** padded. The 'correct' number of headers is equal to the greatest number of columns in any row of the input table:

```
$ cat testdata/table.tab | ./tabular_to_wikitable.py -s -c --header Fear,Surprise
WARNING: Number of column headings (2) less than columns in data (3): padding.
{|class="wikitable sortable mw-collapsible"
|+
! Fear !! Surprise !! col3
|-
| header1 || header2 || header3
|-
| egg || chips || spam
|-
| spam || spam || spam
|-
| Norwegian || Blue || parrot
|}
$ cat testdata/table.tab | ./tabular_to_wikitable.py -s -c --header "Fear,Surprise,Ruthless Efficiency,An Almost Fanatical Devotion To The Pope"
WARNING: Number of column headings (4) greater than columns in data (3).
{|class="wikitable sortable mw-collapsible"
|+
! Fear !! Surprise !! Ruthless Efficiency !! An Almost Fanatical Devotion To The Pope
|-
| header1 || header2 || header3
|-
| egg || chips || spam
|-
| spam || spam || spam
|-
| Norwegian || Blue || parrot
|}
```

Add a title:

```
$ cat testdata/table.tab | ./tabular_to_wikitable.py -s -c --header "Fear,Surprise,Ruthless Efficiency" --title "Cardinal Ximenez's To-Do List"
{|class="wikitable sortable mw-collapsible"
|+ Cardinal Ximenez's To-Do List
! Fear !! Surprise !! Ruthless Efficiency
|-
| header1 || header2 || header3
|-
| egg || chips || spam
|-
| spam || spam || spam
|-
| Norwegian || Blue || parrot
|}
```