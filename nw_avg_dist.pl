#!/usr/bin/perl
# by Peter Menzel
# reads a newick tree file (first arg) and prints for each leaf the average distance to all other leafs
#
use strict;
use Bio::TreeIO;
sub trim($) { my $string = shift; $string =~ s/^\s+//; $string =~ s/\s+$//; return $string; }

if(@ARGV == 0) { die "Give file name (Newick) as first argument"; }

my $treeio = Bio::TreeIO->new('-format' => 'newick', '-file'   => $ARGV[0]);
if( my $tree = $treeio->next_tree ) {
	my @taxa = $tree->get_leaf_nodes;
	foreach my $currLeaf (@taxa) {
			print trim($currLeaf->id),"\t";
			my $distsum;
			foreach my $testLeaf (@taxa) {
				$distsum += $tree->distance(-nodes => [$currLeaf,$testLeaf]);
			}
			printf"%.4f\n", $distsum / @taxa; 
	}

}
else { die "Error reading tree"; }
