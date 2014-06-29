#!/bin/perl
use Bio::AlignIO;
my @files = <*.msf>;
foreach $file (@files) {
    my $alnio = Bio::AlignIO->new( -file=>$file);
    $file =~ s{\.[^.]+$}{};
    open OUTFILE, '>'.$file.'.seq';
    my $started = 0;
    while (my $aln = $alnio->next_aln) {
        foreach $seq ($aln->each_seq) {
            $myseq = $seq->seq;
            $myseq =~ s/\.//g;
	    printf OUTFILE "# ";
	    printf OUTFILE $seq->id();
	    printf OUTFILE "\n";
	    printf OUTFILE $myseq;
	    printf OUTFILE "\n";
        }
    }
    close OUTFILE;
}
