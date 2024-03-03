#!/usr/bin/perl

$infile = $ARGV[0]; # input sam file
$out_prefix = $ARGV[1]; # demultplexed .fa prefix
$out_type = $ARGV[2]; # generate fa or fq files

open IN, "<", $infile or die; #sam file

$ct = 0;
while (my $line = <IN>){
	if ($line =~ /[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+)\t[^\t]+\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)+\t([^\t]+)/){
		my $sid = $1;
		my $pos = $2;
		my $cigar = $3;
		my $qseq = $4;
		my $qual = $5;

		
		if ($out_type eq "fa"){
			my $qid = ">" . $sid . "__" . $ct . "__" . $pos . "__" . $cigar . "\n";
			$hash{$sid}{out} .= $qid . $qseq . "\n";
		} elsif ($out_type eq "fq"){
			my $qid = "@" . $sid . "__" . $ct . "__" . $pos . "__" . $cigar . "\n";
			$hash{$sid}{out} .= $qid . $qseq . "\n+\n" . $qual . "\n";
		}
		
		$hash{$sid}{ct}++; 
		$ct++;
	}
}

foreach my $sid (sort keys %hash){
	print $out_prefix . "\t" . $sid . "\t" . $hash{$sid}{ct} . "\n";
	open OUT, ">", $out_prefix . "__" . $sid . "." . $out_type or die;
	print OUT $hash{$sid}{out};
	close OUT;
}
