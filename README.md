# RILogo

Authors:
Peter Menzel <pmenzel@gmail.com>
Stefan Seemann <seemann@rth.dk>

Copyright 2012,2018 Peter Menzel, see file LICENCE


## Description

RILogo is a command line program to create an RNA-RNA interaction logo from a
pair of RNA sequences or alignments and outputs an SVG file.
If alignments are given as input, then RILogo displays them as sequence logos.
Intramolecular base pairs are displayed by arcs connecting the corresponding
columns of the sequence logo.  Intermolecular base pairs are denoted by
connecting lines between both logos.

RILogo calculates the mutual information for base pairs and displays it through
a colour gradient in the connecting arcs and lines.  Additionally an 'M'
character is placed on top of the other 4 letters in the sequences logos.

The source code of RILogo is available at http://github.com/pmenzel/RILogo

A web server is available: http://rth.dk/resources/rilogo

The paper for RILogo is published as:  
P. Menzel, S.E. Seemann, J. Gorodkin  
[http://www.ncbi.nlm.nih.gov/pubmed/22826541](RILogo: Visualising RNA-RNA interactions)  
Bioinfomatics, 2012.


RILogo is distributed as open source software under the GNU Lesser General
Public Licence 3, see the file LICENCE.


## Installation

RILogo is written in C++ for Linux. It does not depend on
additional libraries.  To compile from source, simply type `make` and the
program will be compiled into the binary file `RILogo`.


## Usage

RILogo expects either one or two input files that contain the sequences as
arguments, or can read a single input file from STDIN.

```
RILogo [options] file1.fa [file2.fa]  >output.svg
```

or

```
RILogo [options] <file.fa >output.svg
```


## Input format

The input is either a single file or two files in FASTA alignment format.
In the former case, the sequences of the two interacting RNAs need to be
separated by the `&` character.  All lines have to be of the same length. A
secondary structure annotation has to be given in a special line with the
name `structure`.  Base pairs are denoted in the dot-bracket notation, see
examples below.

Interactions between the two RNAs can be denoted in the structure line only by
uppercase letters on the left side and corresponding lowercase letters on the
right side of the separator `&`.  Additionally uppercase and lowercase letters
can be used together on one side exclusively to denote pseudoknots.

### Example 1:
Simple internal structure, denoted by brackets, and interaction, denoted by `AA` and `aa`.
```
>seq1
AACGTAACGTAAACGAA&AACGTAAACGAAACGAA
>structure
..(((..BBB..)))..&..(((..bbb..)))..
```

### Example 2:
Internal structure on the left side containing a pseudoknot, which is denoted by `AA` and `aa`.
The interaction is denoted by `BB` and `bb`.
```
>seq1
AGCTAGCTAAGCTAAAGCAAAGCAA&AACGAGCCGTAA
>structure
.AAA.(((..BBB..aaa..)))..&.(((bbb)))..
```

### Example 3:
This example shows the interaction between the bacterial small RNA OxyS and
its binding site in the *fhlA* mRNA. See the files `fhlA-OxyS.fa` for the
alignment, `fhlA-OxyS.tree` for the calculated phylogenetic tree, and
`fhlA-OxyS.treedist` for the tree distances. The output of RILogo with default options using the command
```
RILogo -t fhlA-OxyS.treedist fhlA-OxyS.fa > fhlA-OxyS.svg
```
is in the file `fhlA-OxyS.svg`.


## Command line parameters
```
-m NAME       Specify name of mutual information measure.
              Options are 'mi' and 'miwp'. Default is 'miwp'.
-t FILENAME   Switch to treeMI or treeMI^WP (depending on parameter -m) measure and
              specify the name of the file with the average tree distances.
-w            Use N_d instead of N in the weighting of observed and expected base pair frequencies.
-d            Debug mode
-g            Debug mode for SVG output
-v            Verbose mode, prints calculated MI values to STDERR
```
The script `nw_avg_dist.pl` can be used to calculate the pairwise distances
from a phylogenetic tree in Newick format (uses BioPerl), which can be used
with option `-t` in RILogo.


## Output

RILogo outputs a single vector graphics file in SVG format, which is written to
STDOUT.

The SVG file can easily be converted to other file formats, for example
with one of these commands, either using Inkscape or ImageMagick:
```
inkscape -f in.svg -A out.pdf
inkscape --export-png=out.png --export-width=800 in.svg
convert in.svg out.png
```

## Configuration

RILogo can read a configuration file to customise the output.
See `default.cfg` for the customisable parameters and their
default values in RILogo.


