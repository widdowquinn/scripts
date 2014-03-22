#README.md (scripts/misc)

## Overview
This repository contains a set of scripts for a bunch of disparate purposes. They may or may not be useful to you - I hope they might be, though. The dependencies (and possibly languages) may vary for each script. Please refer to the individual **README** details for each script.

## Scripts
The current set of scripts includes:

* [`excel_to_tab.scpt`](#excel_to_tab): converts an Excel workbook to a subdirectory of tab-separated plaintext files, one per worksheet
* [`rename_to_hash`](#rename_to_hash): renames passed files to theit MD5 hash
* [`tabular_to_wikitable.py`](#tabular_to_wikitable): converts a plain text tab-separated table to a MediaWiki wikitable

## Script READMEs

### <a name="excel_to_tab">`excel_to_tab.scpt`</a>

This script takes as input an Excel workbook containing one or more worksheets. It creates a new directory with the same name as the workbook, with the appended string `_extracted`. This directory contains a set of tab-separated plaintext files, one per worksheet. Each file has the same name as the corresponding worksheet, with the extension `.tab`.

#### Installation

This script is written in AppleScript, and is only expected to work on OSX.

Place the `excel_to_tab.scpt` into your `~/Library/Scripts` directory (create this directory if it does not exist).

Open AppleScript Editor (in `/Applications/Utilities`) and open the General Preferences. Check the `Show Script menu in menu bar` setting, and close AppleScript Editor. You should now see the script symbol in the top menu bar.

![The location of the script symbol, second from right](images/excel_to_tab1.png?raw=True =200x)

![The script menu, showing excel_to_tab](images/excel_to_tab2.png?raw=True =200x)

#### Usage

Click on the script symbol in the menu bar, and select the `excel_to_tab` option (it will be in the lower section). This will open a file selection dialog box. Select the appropriate Excel file, and click `Choose`. The script will generate the output directory in the same location as the Excel file.




### <a name="rename_to_hash">`rename_to_hash`</a>

This script takes a file or list of files as input, and renames them to their MD5 hash, preserving the file extension. This is useful in a limited set of circumstances, admittedly.  I use it for avoiding filename collisions, where the filename itself is unimportant.

#### Usage

```
Usage: rename_to_hash [options] ARGS 
Options:
   -h, --help    Display this message.
   -n            Dry-run; only show what would be done.
```

As an example of usage:

```
$ mkdir temp && cd temp && touch file.{a..e} && mkdir dir_{a..e} && touch f.{a..e}
$ ls
dir_a/  dir_c/  dir_e/  f.b     f.d     file.a  file.c  file.e
dir_b/  dir_d/  f.a     f.c     f.e     file.b  file.d
$ rename_to_hash -n *
dir_a is not a regular file (skipping)
dir_b is not a regular file (skipping)
dir_c is not a regular file (skipping)
dir_d is not a regular file (skipping)
dir_e is not a regular file (skipping)
mv -v f.a d41d8cd98f00b204e9800998ecf8427e.a
mv -v f.b d41d8cd98f00b204e9800998ecf8427e.b
mv -v f.c d41d8cd98f00b204e9800998ecf8427e.c
mv -v f.d d41d8cd98f00b204e9800998ecf8427e.d
mv -v f.e d41d8cd98f00b204e9800998ecf8427e.e
mv -v file.a d41d8cd98f00b204e9800998ecf8427e.a
mv -v file.b d41d8cd98f00b204e9800998ecf8427e.b
mv -v file.c d41d8cd98f00b204e9800998ecf8427e.c
mv -v file.d d41d8cd98f00b204e9800998ecf8427e.d
mv -v file.e d41d8cd98f00b204e9800998ecf8427e.e
$ rename_to_hash *
dir_a is not a regular file (skipping)
dir_b is not a regular file (skipping)
dir_c is not a regular file (skipping)
dir_d is not a regular file (skipping)
dir_e is not a regular file (skipping)
f.a -> d41d8cd98f00b204e9800998ecf8427e.a
f.b -> d41d8cd98f00b204e9800998ecf8427e.b
f.c -> d41d8cd98f00b204e9800998ecf8427e.c
f.d -> d41d8cd98f00b204e9800998ecf8427e.d
f.e -> d41d8cd98f00b204e9800998ecf8427e.e
d41d8cd98f00b204e9800998ecf8427e.a already exists, skipping file.a
d41d8cd98f00b204e9800998ecf8427e.b already exists, skipping file.b
d41d8cd98f00b204e9800998ecf8427e.c already exists, skipping file.c
d41d8cd98f00b204e9800998ecf8427e.d already exists, skipping file.d
d41d8cd98f00b204e9800998ecf8427e.e already exists, skipping file.e
$ ls
d41d8cd98f00b204e9800998ecf8427e.a  dir_a/                              file.a
d41d8cd98f00b204e9800998ecf8427e.b  dir_b/                              file.b
d41d8cd98f00b204e9800998ecf8427e.c  dir_c/                              file.c
d41d8cd98f00b204e9800998ecf8427e.d  dir_d/                              file.d
d41d8cd98f00b204e9800998ecf8427e.e  dir_e/                              file.e
```


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