#!/usr/bin/perl
## a simple module for reading in a standard fasta file
## in the ">" line, only the part in front of the first space is used as seq name
## allow sequences to be split into multiple lines
## return a hash variable
## hebin
## 6 sep 2011
##

#package fastaIN;
#use Exporter;
#@ISA = ('Exporter');
#@EXPORT = ('fastaIN');

sub fastaIN {
    my $usage = "fastaIN('infile.fa')\n";
    my $infile = shift or die $usage;
    open(IN, $infile) or 
	die "Cannot find or open $infile\n";
    my @seqs;
    $/ = ">";
    my @tmp = <IN>;
    $/ = "\n";
    close IN;
    shift @tmp; ## skip the first ">"
    foreach my $line (@tmp) {
	$line =~ s/>//;
	my @array = split(/\n/,$line); # [0] seq name [1..end] seq (may be split into multiple lines)
	my $seqName = shift @array;
	my $fasta = join("",@array);
	push(@seqs,($seqName,$fasta));
    }
    return @seqs;
}

1;
