#! /usr/bin/perl -w
use strict;
use warnings;

use IO::File;
use FileHandle;


if( ! $ARGV[ 0 ] )
{
		print "no input\n\n";

}
elsif( $ARGV[ 0 ] eq "--detect_kmer_motifs" )
{
	&detect_kmer_motifs ($ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4], $ARGV[5]);
}
elsif( $ARGV[ 0 ] eq "--annotate_pCRE" )
{
	&annotate_pCRE ($ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4]);
}


##################################################################
## detect motifs with k-mer length in observed and control sequences
##################################################################
sub detect_kmer_motifs
{
	my ($kmer, $dir, $infile_control, $infile_test, $outfile) = @_;
	$infile_control = "$dir$infile_control";
	$infile_test = "$dir$infile_test";
	$outfile = "$dir$outfile";

	# get the 6-mer strings
	print "get the $kmer"."-mer strings\n";
	my @merKlist = ();
	my @bases4 = ('A', 'T', 'C', 'G');
	if ($kmer == 6)
	{
		foreach my $i1 (0..$#bases4)
		{
			foreach my $i2 (0..$#bases4)
			{
				foreach my $i3 (0..$#bases4)
				{
					foreach my $i4 (0..$#bases4)
					{
						foreach my $i5 (0..$#bases4)
						{
							foreach my $i6 (0..$#bases4)
							{
								my $string = qq{$bases4[$i1]$bases4[$i2]$bases4[$i3]$bases4[$i4]$bases4[$i5]$bases4[$i6]};
								push @merKlist, $string;
							}
						}
					}
				}
			}
		}
	} # 6-mer
	elsif ($kmer == 7)
	{
		foreach my $i1 (0..$#bases4)
		{
			foreach my $i2 (0..$#bases4)
			{
				foreach my $i3 (0..$#bases4)
				{
					foreach my $i4 (0..$#bases4)
					{
						foreach my $i5 (0..$#bases4)
						{
							foreach my $i6 (0..$#bases4)
							{
								foreach my $i7 (0..$#bases4)
								{
									my $string = qq{$bases4[$i1]$bases4[$i2]$bases4[$i3]$bases4[$i4]$bases4[$i5]$bases4[$i6]$bases4[$i7]};
									push @merKlist, $string;
								}							
							}
						}
					}
				}
			}
		}
	} # 7-mer
	elsif ($kmer == 8)
	{
		foreach my $i1 (0..$#bases4)
		{
			foreach my $i2 (0..$#bases4)
			{
				foreach my $i3 (0..$#bases4)
				{
					foreach my $i4 (0..$#bases4)
					{
						foreach my $i5 (0..$#bases4)
						{
							foreach my $i6 (0..$#bases4)
							{
								foreach my $i7 (0..$#bases4)
								{
									foreach my $i8 (0..$#bases4)
									{
										my $string = qq{$bases4[$i1]$bases4[$i2]$bases4[$i3]$bases4[$i4]$bases4[$i5]$bases4[$i6]$bases4[$i7]$bases4[$i8]};
										push @merKlist, $string;
									}
								}							
							}
						}
					}
				}
			}
		}
	} # 8-mer
	elsif ($kmer == 9)
	{
		foreach my $i1 (0..$#bases4)
		{
			foreach my $i2 (0..$#bases4)
			{
				foreach my $i3 (0..$#bases4)
				{
					foreach my $i4 (0..$#bases4)
					{
						foreach my $i5 (0..$#bases4)
						{
							foreach my $i6 (0..$#bases4)
							{
								foreach my $i7 (0..$#bases4)
								{
									foreach my $i8 (0..$#bases4)
									{
										foreach my $i9 (0..$#bases4)
										{	
											my $string = qq{$bases4[$i1]$bases4[$i2]$bases4[$i3]$bases4[$i4]$bases4[$i5]$bases4[$i6]$bases4[$i7]$bases4[$i8]$bases4[$i9]};
											push @merKlist, $string;
										}										
									}
								}							
							}
						}
					}
				}
			}
		}
	} # 9-mer
	elsif ($kmer == 10)
	{
		foreach my $i1 (0..$#bases4)
		{
			foreach my $i2 (0..$#bases4)
			{
				foreach my $i3 (0..$#bases4)
				{
					foreach my $i4 (0..$#bases4)
					{
						foreach my $i5 (0..$#bases4)
						{
							foreach my $i6 (0..$#bases4)
							{
								foreach my $i7 (0..$#bases4)
								{
									foreach my $i8 (0..$#bases4)
									{
										foreach my $i9 (0..$#bases4)
										{	
											foreach my $i10 (0..$#bases4)
											{	
												my $string = qq{$bases4[$i1]$bases4[$i2]$bases4[$i3]$bases4[$i4]$bases4[$i5]$bases4[$i6]$bases4[$i7]$bases4[$i8]$bases4[$i9]$bases4[$i10]};
												push @merKlist, $string;
											}
										}
									}
								}							
							}
						}
					}
				}
			}
		}
	} # 10-mer

	my $Nstrings = scalar(@merKlist);
	print "$Nstrings $kmer -mer strings\n";
	
	# read test fasta file
	print "read test fasta file\n";
	my $IN_test = new IO::File ("$infile_test") || die "Fail to open file: $infile_test !\n";
    $/ = ">";
	my $line = <$IN_test>;
	my %seqs_test = ();
	while (<$IN_test>)
    {
		$_ =~ s/>$//;
	#	print $_;
		my @tmparr = split(/\n/, $_);
		my $title = $tmparr[0];
		shift @tmparr;
		my $seq = join("", @tmparr);
		$seq =~ s/^\s+//;
		$seq =~ s/\s+$//;
		die "$title already in test seqs\n" if exists $seqs_test{$title};
		$seqs_test{$title} = $seq;
	}
	$/ = "\n";
	close $IN_test;

	# read control fasta file
	print "read control fasta file\n";
	my $IN_control = new IO::File ("$infile_control") || die "Fail to open file: $infile_control !\n";
    $/ = ">";
	$line = <$IN_control>;
	my %seqs_control = ();
	while (<$IN_control>)
    {
		$_ =~ s/>$//;
	#	print $_;
		my @tmparr = split(/\n/, $_);
		my $title = $tmparr[0];
		shift @tmparr;
		my $seq = join("", @tmparr);
		$seq =~ s/^\s+//;
		$seq =~ s/\s+$//;
		die "$title already in test seqs\n" if exists $seqs_control{$title};
		$seqs_control{$title} = $seq;
	}
	$/ = "\n";
	close $IN_control;

	# calculate gene counts for each 6-mer string
	print "calculate gene counts for each 6-mer string and save\n";
    open OUT, ">$outfile" || die "Fail to open file: $outfile !\n";
    print OUT "id\tkmer\tNmer_t\tNmer_c\tNother_t\tNother_c\n";
	my $Ntotal_t = scalar(keys %seqs_test);
 	my $Ntotal_c = scalar(keys %seqs_control);
	my $i = 0;
   	foreach my $this_kmer (@merKlist)
    {
		$i ++;

    	# in test
		my $count_t = 0;
		foreach my $t (keys %seqs_test)
		{
			my $seq = $seqs_test{$t};
			if ($seq =~ m/$this_kmer/)
			{
				$count_t ++;
			}
		}
		# in control
		my $count_c = 0;
		foreach my $t (keys %seqs_control)
		{
			my $seq = $seqs_control{$t};
			if ($seq =~ m/$this_kmer/)
			{
				$count_c ++;
			}
		}
		my $miss_t = $Ntotal_t - $count_t;
		my $miss_c = $Ntotal_c - $count_c;

    	print OUT qq{$i\t$this_kmer\t$count_t\t$count_c\t$miss_t\t$miss_c\n};
    }
	close OUT;
}

##################################################################
## 1. transform dna seq of motifs to meme format
## 2. run tomtom to annotate k-mer motifs
##################################################################
sub annotate_pCRE
{
	my ($dir, $tfname, $infile_motifslist, $dir_motifdb) = @_;
	my $infile = "$dir$tfname.p";
	my $outfile = "$dir$tfname.anno";

	# read motif names
	my $INnames = new IO::File ("input/reg_motifsDB.list") || die "Fail to open file: reg_motifsDB.list !\n";
	my %names = ();
	while(<$INnames>)
	{
		chomp;
		my ($l, $id_db, $name) = split(/\s+/, $_);
		$names{$id_db} = $name;
	}
	close $INnames;
	print scalar(keys %names)." id-name pairs\n";

	# read motif file
	print "read motif file\n  $infile\n  $outfile\n";
	my $IN = new IO::File ("$infile") || die "Fail to open file: $infile !\n";
	my $OUT = new IO::File ("$outfile", "w") || die "Fail to open file: $outfile !\n";
	my $line = <$IN>;
	my $cmd;
	print $OUT qq{miid\tkmer\tmotif\tNmer_t\tNmer_c\tNother_t\tNother_c\toddratio\tp\tadjp\t}.
		qq{Target_ID\tOptimal_offset\tanno_p\tanno_E\tanno_q\tOverlap\tTarget_consensus\tOrientation\ttargetName\n};
	while (<$IN>)
    {
		# id      kmer    Nmer_t  Nmer_c  Nother_t        Nother_c        odds.ratio      
		# p.value.fish    p.adjust.fish
		chomp;
	#	print $_;
		my ($id, $motif, $Nmer_t, $Nmer_c, $Nother_t, $Nother_c,
			$oddratio, $p, $adjp) = split(/\t/, $_);
		my $k = length($motif);
		print "($id, $motif, $Nmer_t, $Nmer_c, $Nother_t, $Nother_c, $oddratio, $p, $adjp)\n";
		# 1  transform dna seq of motifs to meme format: iupac2meme -dna AACAAA > mo1.meme
		my $iupac2meme_path = qx{which iupac2meme};
		chomp($iupac2meme_path);
		if (!-e $iupac2meme_path) {
    		die "Error: iupac2meme not found in PATH.\n";
		}

		$cmd = qq{$iupac2meme_path -dna $motif > $dir/$tfname/$id.meme};
		`$cmd`;

		# 2 run tomtom to annotate k-mer motifs
		my $tomtom_path = qx{which tomtom};
		chomp($tomtom_path);
		if (!-e $tomtom_path) {
    		die "Error: tomtom not found in PATH.\n";
		}
		$cmd = qq{$tomtom_path $dir/$tfname/$id.meme  \\
	$dir_motifdb/JASPAR/JASPAR2022_CORE_plants_non-redundant_v2.meme  \\
	$dir_motifdb/CIS-BP_2.00/Arabidopsis_thaliana.meme  \\
	$dir_motifdb/CIS-BP_2.00/Arabidopsis_lyrata.meme  \\
	$dir_motifdb/CIS-BP_2.00/Capsella_rubella.meme  \\
	$dir_motifdb/CIS-BP_2.00/Brassica_napus.meme  \\
	$dir_motifdb/CIS-BP_2.00/Brassica_oleracea.meme  \\
	$dir_motifdb/CIS-BP_2.00/Brassica_rapa.meme  \\
	$dir_motifdb/CIS-BP_2.00/Eutrema_salsugineum.meme  \\
	$dir_motifdb/CIS-BP_2.00/Carica_papaya.meme  \\
	$dir_motifdb/CIS-BP_2.00/Zea_mays.meme  \\
	$dir_motifdb/CIS-BP_2.00/Oryza_sativa.meme  \\
	$dir_motifdb/ARABD/ArabidopsisDAPv1.meme  \\
	$dir_motifdb/ARABD/ArabidopsisPBM_20140210.meme  \\
	-oc outdir_$tfname };
#		print "$cmd\n";
		`$cmd`;
	
		$cmd = qq{mv outdir_$tfname/tomtom.tsv $dir/$tfname/$id.tomtom.tsv};
		`$cmd`;
		$cmd = qq{mv outdir_$tfname/tomtom.html $dir/$tfname/$id.tomtom.html};
		`$cmd`;

		# 3 combine motifs info and annotation, append motif name
		my $INtomtom = new IO::File ("$dir/$tfname/$id.tomtom.tsv") || die "Fail to open file: $dir/$tfname/$id.tomtom.tsv !\n";
		$line = <$INtomtom>;
		while(<$INtomtom>)
		{
			# Query_ID        Target_ID       Optimal_offset  
			# p-value E-value q-value Overlap Query_consensus Target_consensus  Orientation
			chomp;
			next if $_ =~ m/^#/;
			next if $_ =~ m/^\s*$/;
			my ($query, $target, $offset, $anno_p, $anno_E, $anno_q, 
				$overlap, $query_con, $target_con, $orientation) = split(/\t/, $_);
			my $name = 'NA';
			if(exists $names{$target})
			{
				$name = $names{$target};
			}
			print $OUT qq{$id\t$k\t$motif\t$Nmer_t\t$Nmer_c\t$Nother_t\t$Nother_c\t$oddratio\t$p\t$adjp\t}.
				qq{$target\t$offset\t$anno_p\t$anno_E\t$anno_q\t$overlap\t$target_con\t$orientation\t$name\n};
		}
		close $INtomtom;
	}
	close $IN;
	close $OUT;

}

