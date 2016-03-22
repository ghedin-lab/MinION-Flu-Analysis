#!/usr/bin/perl -w

use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq;

my $file = Bio::SeqIO->new(-file=>"$ARGV[0]");

while(my $seq_obj = $file->next_seq){

        print $seq_obj->id,"\t",length($seq_obj->seq),"\n";

}
