#!/bin/bash

RILOGO='.'
FASTTREE="/home/ptr/software/bin"

FASTA1=$1
NAME1=${FASTA1##*/}
NAME1=${NAME1%.*}

FASTA2=$2
NAME2=${FASTA2##*/}
NAME2=${NAME2%.*}

# calculate phylogenetic tree 
$FASTTREE/FastTree -gtr -nosupport -nt $FASTA1 > ${NAME1}.tmp
# get average distance for each sequence to all other sequences
$RILOGO/nw_avg_dist.pl ${NAME1}.tmp > ${NAME1}.treedist
# run RILogo
$RILOGO/RILogo $FASTA1 $FASTA2 -t ${NAME1}.treedist  > ${NAME1}_${NAME2}.svg
